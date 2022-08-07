#!/usr/bin/python3
import os
import umap
import pickle
import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from re import sub, compile
from tqdm import tqdm
from numpy import linspace
from shutil import copy
from typing import Tuple, List
from tempfile import mkdtemp
from functools import wraps
from detect_delimiter import detect
from gzip import open as gzip_open
from platform import system as platform_check
from multiprocessing import Pool, cpu_count


platform = platform_check()
pat_header1 = compile(rb'^##')
pat_header2 = compile(rb'^#')


# plink_dir = 'D:/Users/Sheng/Documents/Github/pageant/bin'  # Set the plink directory
# save_dir = 'D:/Users/Sheng/Documents/Github/pageant/output/population_QC'  # Set the save directory
# output_dir = 'D:/Users/Sheng/Documents/Github/pageant/output'  # Set the output directory
# ps_ref = 'D:/Users/Sheng/Documents/Github/pageant/personal_genome/hapmap3.vcf.gz'
# prune = True
# temp_dir = 'D:/Users/Sheng/Documents/Github/pageant/temp'
# auto_sep = '\t'
# IID_name = 'IID'
# POP_name = 'population'


def auto_sep(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            fun_res = func(*args, **kwargs)
        except AssertionError as e:
            fun_res = func(*args, **kwargs, sep=detect(e.args[0]))
        except Exception as e:
            raise e
        return fun_res

    return wrapper


class PLINKERROR(Exception):
    def __init__(self, *args):
        super(PLINKERROR, self).__init__(*args)


class RUNPLINKERROR(Exception):
    def __init__(self, *args):
        super(RUNPLINKERROR, self).__init__(*args)


def select_list(ob_list: list, index) -> list:
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def detect_plink(plink: str = 'plink' or 'plink2'):
    a = subprocess.Popen(f'which {plink}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, sterr = a.communicate()
    a.wait()
    return True if stdout else False


def cal_sha(file: str):
    """
    calculate the sha256 of a file
    :param file: file path
    :return: the sha256 value of the file
    """
    from hashlib import sha256
    with open(file, 'rb') as f:
        return sha256(f.read()).hexdigest()


def verify_data(data_file: str, source_files: str or list, mode: str) -> bool:
    if mode == 'single':
        sha_cal = cal_sha
        with open(data_file, encoding='UTF-8') as f:
            f.read(2)
            if f.read(4) == 'DONE':
                if f.read(61) == sha_cal(source_files)[3:]:
                    return True
        return False
    else:
        return False


def open_gz(file: str):
    if file.split('.')[-1] == 'gz':
        return gzip_open
    else:
        return open


def process_line(line: bytes) -> Tuple[np.ndarray, str]:
    line = line.strip().split(b'\t')
    geno_dat = line[6:]
    return np.array([int(geno.decode()) if geno != b'NA' else None for geno in geno_dat], dtype='f'), line[1].decode()


def add_res(geno_res: np.ndarray, iid: str, data: list, iid_list: list, pbar: tqdm) -> None:
    data.append(geno_res)
    iid_list.append(iid)
    pbar.update(1)


@auto_sep
def get_columns_num(file: str, need_columns: list, sep: str = '\t') -> List[int]:
    nopen = open_gz(file)
    with nopen(file, 'rb') as f:
        header = f.readline()
    header = header.strip().split(sep.encode())
    assert len(header) > 1, sep.encode().join(header).decode()
    return [header.index(i.encode()) + 1 for i in need_columns]


@auto_sep
def select_columns(file: str, need_col: list, out_file: str, method: str = 'number', header: bool = True,
                   skipspace: bool = False, sep: str = '\t') -> None:
    first = True
    if method == 'name':
        need_col = get_columns_num(file, need_col)
    elif method != 'number':
        raise AttributeError("Unsupported method")
    nopen = open_gz(file)
    with nopen(file, 'rb') as fr:
        with open(out_file, 'wb') as fw:
            for line in fr:
                if skipspace:
                    line = line.lstrip(b' ')
                    line = sub(f' *{sep} *'.encode(), sep.encode(), line)
                if not first or header:
                    if line.strip():
                        fw.write(sep.encode().join(select_list(line.strip().split(sep.encode()), need_col)) + b'\n')
                if first:
                    assert len(line.strip().split(sep.encode())) > 1, \
                        sep.encode().join(line.strip().split(sep.encode())).decode()
                    first = False


def run_plink_cmd(cmd: str, plink='plink2', delete_log=True) -> None:
    """
    Creat a subprocess which not display in GUI.
    :param delete_log:
    :param plink: plink name
    :param cmd: cmd string
    :return: A finished subprocess with no result
    """
    if '--out' in cmd:
        out_file = cmd.split('--out ')[-1].split(' ')[0]
    else:
        out_file = os.path.join('.', plink)
    if detect_plink(plink):
        plink_file = plink
    else:
        sys_suffix = '.exe' if platform == 'Windows' else ''
        plink_file = os.path.join(plink_dir, plink + sys_suffix)
        assert os.path.isfile(plink_file), f"Cannot find {plink}{sys_suffix} in plink directory"
    cmd = plink_file + ' ' + cmd
    a = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, sterr = a.communicate()
    a.wait()
    if b'Error' in sterr:
        raise PLINKERROR("PlinkError:\n" + sterr.decode().strip())
    if delete_log:
        try:
            os.remove(out_file + '.log')
        except FileNotFoundError:
            raise RUNPLINKERROR("RunPlinkError:\n" + sterr.decode().strip())


def file_recognize(file: str) -> Tuple[str, str]:
    file_extension = os.path.basename(file).replace('.gz', '').split('.')[-1]
    if file_extension == 'vcf':
        file_prefix = 'vcf'
        return file_prefix, file
    elif file_extension == 'bed':
        file_prefix = 'bfile'
    elif file_extension == 'pgen':
        file_prefix = 'pfile'
        if os.path.isfile('.'.join(file.split('.')[:-1]) + '.pvar.zst'):
            return file_prefix, '.'.join(file.split('.')[:-1]) + ' vzs'
    elif file_extension == 'txt':
        file_prefix = '23file'
        return file_prefix, file
    else:
        raise TypeError('Unrecognized file extension')
    file_name = '.'.join(file.split('.')[:-1])
    return file_prefix, file_name


def ps_verify(ps_file: str, work_dir: str) -> bool:
    ps_save = os.path.join(work_dir, 'ps_data')
    if os.path.isfile(ps_save):
        if verify_data(ps_save, ps_file, 'single'):
            return True
    return False


def get_vcf_header(vcf_list: list, binary=False) -> str or bytes:
    for vcf_file in vcf_list:
        nopen = open_gz(vcf_file)
        if nopen:
            with nopen(vcf_file, 'rb') as f:
                for line in f:
                    if pat_header2.search(line):
                        if not pat_header1.search(line):
                            return line if binary else line.decode()


def fill_missing_alt(vcf: str, reference_alt: str, output: str) -> str:
    ref_data = pd.read_csv(reference_alt, sep='\t', dtype={0: str})
    header = get_vcf_header([vcf])
    raw_data = pd.read_csv(vcf, sep='\t', dtype={0: str}, comment='#', header=None, names=header.strip().split('\t'))
    raw_data = raw_data.loc[raw_data['REF'] != 'N', :]
    fill_data = pd.merge(raw_data, ref_data, how='left', on='ID', suffixes=("", "_ref"))
    fill_data['ALT'] = fill_data.apply(
        lambda a: a['ALT'] if a['ALT'] != '.' else a['ALT_ref'] if a['REF'] == a['REF_ref'] else a['REF_ref'], axis=1)
    fill_data.iloc[:, 0:10].to_csv(output, sep='\t', index=False)
    return output


def get_ps_need(ps_file: str, work_dir: str, temp_dir: str, prune: bool) -> None:
    with open(os.path.join(work_dir, 'ps_data'), 'w', encoding='UTF-8') as f:
        f.write(f'###{cal_sha(ps_file)}\n')
        file_prefix, file_name = file_recognize(ps_file)
        if prune:
            # prune data
            run_plink_cmd(f'--{file_prefix} {file_name} --maf 0.01 --geno 0.01 --hwe 1e-50 --indep-pairwise 50 10 0.2 '
                          f'--out {os.path.join(temp_dir, "prune_pop")}', plink='plink2')
        # produce PCA outcome and frequency data
        get_pca_need(file_prefix, file_name, temp_dir, prune, work_dir)
        get_umap_need(file_prefix, file_name, temp_dir, prune, work_dir)
        f.seek(2)
        f.write('DONE')


def get_pca_need(file_prefix: str, file_name: str, temp_dir: str, prune: bool, save_dir: str) -> None:
    convert_data(file_prefix, file_name, temp_dir, prune, 'pca')
    run_plink_cmd(f'--pfile {os.path.join(temp_dir, "pop")} vzs --freq counts --pca 2 allele-wts '
                  f'--out {os.path.join(temp_dir, "ref_pca")}', plink='plink2')
    if prune:
        copy(os.path.join(temp_dir, 'prune_pop.prune.in'), save_dir)
    copy(os.path.join(temp_dir, 'ref_pca.eigenvec.allele'), save_dir)
    copy(os.path.join(temp_dir, 'ref_pca.eigenvec'), save_dir)
    copy(os.path.join(temp_dir, 'ref_pca.eigenval'), save_dir)
    copy(os.path.join(temp_dir, 'ref_pca.acount'), save_dir)


def convert_data(file_prefix: str, file_name: str, temp_dir: str, prune: bool, mode: str) -> None:
    if prune:
        run_plink_cmd(f'--{file_prefix} {file_name} --maf 0.01 --geno 0.01 --hwe 1e-50 --indep-pairwise 50 10 0.2 '
                      f'--out {os.path.join(temp_dir, "prune_pop")}', plink='plink2')
    plink_cmd_prune = f'--extract {os.path.join(temp_dir, "prune_pop.prune.in")} ' if prune else ''
    plink_cmd_recode = '--make-pgen vzs ' if mode != 'umap' else '--recode A '
    run_plink_cmd(
        f'--{file_prefix} {file_name} '
        f'--maf 0.01 '
        f'--geno 0.01 '
        f'--hwe 1e-50 '
        f'--autosome '
        f'{plink_cmd_prune}'
        f'--mac 1 '
        f'{plink_cmd_recode}'
        f'--out {os.path.join(temp_dir, "pop")}', plink='plink2')


def get_umap_need(file_prefix: str, file_name: str, temp_dir: str, prune: bool, save_dir: str) -> None:
    convert_data(file_prefix, file_name, temp_dir, prune, 'umap')
    # Umap
    with open(f'{os.path.join(temp_dir, "pop.raw")}', 'rb') as f:
        total = 1
        dat = []
        freq_mean = []
        iid_list = []
        with tqdm(desc='Getting genotype data...', unit=' samples', unit_scale=True) as pbar:
            with Pool(processes=cpu_count()) as pool:
                rs_id = pd.DataFrame([rs.split('_') for rs in f.readline().decode().strip().split('\t')[6:]])
                pool.apply_async(process_line, args=(f.readline(),),
                                 callback=lambda a: add_res(*a, dat, iid_list, pbar))
                for line in f:
                    total += 1
                    pbar.format_meter(0, total, 0)
                    pool.apply_async(process_line, args=(line,), callback=lambda a: add_res(*a, dat, iid_list, pbar))
                pool.close()
                pool.join()
    dat = np.vstack(tuple(dat))
    for col in range(dat.shape[1]):
        freq_mean.append(np.nanmean(dat[:, col]))
        dat[:, col] = np.nan_to_num(dat[:, col], nan=freq_mean[col], copy=False)
    rs_id.to_csv(os.path.join(save_dir, "snp_info"), sep='\t', index=False, header=False)
    rs_id.insert(2, 'freq_mean', freq_mean)
    rs_id.to_csv(os.path.join(save_dir, "umap_snp"), sep='\t', index=False, header=False)
    ref_data = pd.DataFrame([iid.split('_') for iid in iid_list] if iid_list[0].count('_') == 1 else iid_list)
    ref_data.columns = ['FID', 'IID'] if ref_data.shape[1] > 1 else ['IID']
    ref_data.to_csv(os.path.join(save_dir, "sample_info"), index=False)
    umap_fit = umap.UMAP(random_state=42)
    umap_fit.fit(dat)
    with gzip_open(os.path.join(save_dir, 'umap_fit'), 'wb', compresslevel=6) as f:
        f.write(pickle.dumps(umap_fit))


def pca_sample(sample_vcf: str, work_dir: str, temp_dir: str) -> Tuple[float, float]:
    # extract data
    plink_prune_cmd = f'--extract {os.path.join(work_dir, "prune_pop.prune.in")} ' if prune else ''
    run_plink_cmd(f'--vcf {sample_vcf} '
                  f'--autosome '
                  f'{plink_prune_cmd}'
                  f'--recode vcf '
                  f'--out {os.path.join(temp_dir, "sample_prune_in")}',
                  plink='plink2')
    # Fill missing Alt
    fill_vcf = fill_missing_alt(os.path.join(temp_dir, "sample_prune_in.vcf"), os.path.join(work_dir, "ref_pca.acount"),
                                os.path.join(temp_dir, "fill.vcf"))
    # Calculate sample's PC values
    run_plink_cmd(f'--vcf {fill_vcf} --read-freq {os.path.join(work_dir, "ref_pca.acount")} '
                  f'--max-alleles 2 '
                  f'--score {os.path.join(work_dir, "ref_pca.eigenvec.allele")} 2 5 header-read no-mean-imputation '
                  f'variance-standardize '
                  f'--score-col-nums 6-7 '
                  f'--out {os.path.join(temp_dir, "sample_pca")}', plink='plink2')
    sample_pca = pd.read_csv(os.path.join(temp_dir, 'sample_pca.sscore'), sep='\t')
    pca_vec = pd.read_csv(os.path.join(work_dir, 'ref_pca.eigenval'), sep='\t', header=None)[0]
    return sample_pca['PC1_AVG'] / (-(pca_vec[0] ** 0.5) / 2), sample_pca['PC2_AVG'] / (-(pca_vec[1] ** 0.5) / 2)


def get_metadata(metadata: str or None, ref_data: pd.DataFrame) -> pd.DataFrame:
    if metadata:
        with open(metadata) as f:
            population_sep = detect(f.readline(10))
        population_data = pd.read_csv(metadata, sep=population_sep)
        return pd.merge(ref_data, population_data, left_on='IID', right_on=IID_name)
    else:
        ref_data.insert(0, POP_name, 'Reference')
        return ref_data


def pca_plot(sample_vcf: str or None, metadata: str or None, work_dir: str, temp_dir: str, img_dir: str) -> None:
    pca_data = pd.read_csv(os.path.join(work_dir, 'ref_pca.eigenvec'), sep='\t')
    if '#IID' in pca_data.columns:
        pca_data['IID'] = pca_data['#IID'].apply(lambda a: a.split('_')[-1])
    ref_pca_data = get_metadata(metadata, pca_data)
    sample_pc = pca_sample(sample_vcf, work_dir, temp_dir) if sample_vcf else None
    draw_plot(ref_pca_data, sample_pc, POP_name, img_dir, 'PC')


def draw_plot(ref_data: pd.DataFrame, sample: Tuple[float, float] or None,
              pop_column: str, img_dir: str, mode: str = 'PC' or 'UMAP') -> None:
    grouped = ref_data.groupby(pop_column)
    colors = iter(cm.rainbow(linspace(1, 0, len(grouped))))
    plt.figure(dpi=400)
    for key, group in grouped:
        plt.scatter(x=group[f'{mode}1'], y=group[f'{mode}2'], label=key, s=66, alpha=0.35, marker='.',
                    color=next(colors))
    if sample:
        plt.scatter(x=sample[0], y=sample[1],
                    label='Me', s=66, marker='x', c='black')
    plt.legend(bbox_to_anchor=(1.05, 0.95), loc='upper left', borderaxespad=0, ncol=len(grouped) // 18 + 1)
    plt.xlabel(f'{mode}1', fontdict={'weight': 'bold'})
    plt.ylabel(f'{mode}2', fontdict={'weight': 'bold'})
    plt.savefig(os.path.join(img_dir, f'QC_{"PCA" if mode == "PC" else mode}.png'), bbox_inches='tight')
    plt.close()


def sample_umap(sample_vcf: str, umap_fit: umap.UMAP, work_dir: str, temp_dir: str) -> Tuple[float, float]:
    plink_prune_cmd = f'--extract {os.path.join(work_dir, "prune_pop.prune.in")} ' if prune else ''
    run_plink_cmd(
        f'--vcf {sample_vcf} '
        f'--autosome '
        f'{plink_prune_cmd}'
        f'--export A '
        f'--export-allele {os.path.join(work_dir, "snp_info")} '
        f'--out {os.path.join(temp_dir, "sample_umap")}', plink='plink2')
    with open(os.path.join(temp_dir, 'sample_umap.raw'), 'rb') as f:
        snps = f.readline().strip().split(b'\t')[6:]
        geno = f.readline().strip().split(b'\t')[6:]
    sample_dat = [snp.decode().split('_') + [geno[idx].decode()] for idx, snp in enumerate(snps)]
    sample_dat = [sample for sample in sample_dat if sample[2].isnumeric()]
    rs_id = pd.read_csv(os.path.join(work_dir, "umap_snp"), header=None, sep='\t')
    sample = pd.merge(rs_id, pd.DataFrame(sample_dat), how='left', left_on=[0], right_on=[0])
    sample.insert(5, 'value',
                  sample.apply(lambda a: a['2_x'] if pd.isna(a['1_y']) else
                  a['2_y'] if a['1_x'] == a['1_y'] else 2 - a['2_y'], axis=1))
    own_umap = umap_fit.transform(np.array([sample['value']]))
    return own_umap[0, 0], own_umap[0, 1]


def umap_plot(sample_vcf: str or None, metadata: str or None, work_dir: str, temp_dir: str, img_dir: str) -> None:
    with gzip_open(os.path.join(work_dir, 'umap_fit')) as f:
        umap_fit = pickle.load(f)
    umap_data = pd.read_csv(os.path.join(work_dir, "sample_info"))
    umap_data = pd.concat([umap_data, pd.DataFrame(umap_fit.embedding_, columns=['UMAP1', 'UMAP2'])], axis=1)
    sample_value = sample_umap(sample_vcf, umap_fit, work_dir, temp_dir) if sample_vcf else None
    ref_umap_data = get_metadata(metadata, umap_data)
    draw_plot(ref_umap_data, sample_value, POP_name, img_dir, 'UMAP')


def ps_plot(sample_vcf: str or None, metadata: str or None, work_dir: str, temp_dir: str, img_dir: str) -> None:
    pca_plot(sample_vcf, metadata, work_dir, temp_dir, img_dir)
    umap_plot(sample_vcf, metadata, work_dir, temp_dir, img_dir)


def ps_analyse(ps_ref: str, sample_vcf: str or None, metadata: str or None,
               temp_dir: str, save_dir: str, img_dir: str or None) -> str:
    if not img_dir:
        img_dir = save_dir
    if not ps_verify(ps_ref, save_dir):
        get_ps_need(ps_ref, save_dir, temp_dir, prune)
    ps_plot(sample_vcf, metadata, save_dir, temp_dir, img_dir)
    return 'Yes'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and UMAP analysis for genotype data')
    parser.add_argument('-s', '--sample', default=None,
                        help='The sample genotype data (VCF (*.vcf *.vcf.gz); 23andme (*.txt); '
                             'PLINK 1 binary (*.bed)); PLINK 2 binary (*.pgen)')
    parser.add_argument('-p', '--population', required=True, dest='population',
                        help='The population genotype data (vcf (*.vcf *.vcf.gz); '
                             'PLINK 1 binary (*.bed); PLINK 2 binary (*.pgen))')
    parser.add_argument('-m', '--metadata', default=None, help='The metadata for population file')
    parser.add_argument('-o', '--output', default='./output', help='The output directory')
    parser.add_argument('-t', '--temp', default=None, help='The temporary directory')
    parser.add_argument('--plink-dir', default='./bin', dest='plink',
                        help='The directory for storing PLINK executable files')
    parser.add_argument('--image', default=None, help='The directory where the images save')
    parser.add_argument('--sep', default='\t', help='The delimiter for metadata')
    parser.add_argument('--pop-col', default='population', dest='pop',
                        help="The name of the column that describe population infomation")
    parser.add_argument('--id-col', default='IID', dest='iid',
                        help="The name of the column that describe individual id")
    parser.add_argument('--prune', dest='prune', action='store_true')
    parser.add_argument('--no-prune', dest='prune', action='store_false')
    parser.set_defaults(prune=True)
    args = parser.parse_args()
    prune = args.prune
    # auto_sep = args.sep
    IID_name = args.iid
    POP_name = args.pop
    plink_dir = args.plink
    if not args.temp:
        args.temp = mkdtemp(suffix='pageant')
    ps_analyse(args.population, args.sample, args.metadata, args.temp, args.output, args.image)
    print('Finish!')
