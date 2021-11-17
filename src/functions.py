import os.path

import umap
from functools import partial
from gzip import open as gzip_open
from multiprocessing import Pool, cpu_count
from pickle import dump, load, dumps
from re import compile, sub
from tempfile import mkdtemp
from urllib.request import urlretrieve
from zipfile import ZipFile
import matplotlib.cm as cm
import numpy as np
from tqdm import tqdm
from objects import *

pattern = compile(rb'(?<=[\t])rs[0-9]*(?=[\t;])')
pat_header1 = compile(rb'^##')
pat_header2 = compile(rb'^#')
pattern_clinvar = compile(r'(?<=CLNSIG=)\w+(?=[,;])')
pattern_clinvar_rs = compile(r'(?<=RS=)\d+')
split_str = compile(r'\w{1,28}')
pattern_load_num = compile(r'\d+(?= variants processed\.)')
sample_file_pat = compile(r'\.samples$')
functions_config = configparser.ConfigParser()


# Sub_functions
def load_config_fun(main_config: configparser.ConfigParser):
    global functions_config
    functions_config = main_config
    load_config_obj(main_config)


def open_gz(file: str):
    if file.split('.')[-1] == 'gz':
        return gzip_open
    else:
        return open


def cal_sha(file: str):
    """
    calculate the sha256 of a file
    :param file: file path
    :return: the sha256 value of the file
    """
    from hashlib import sha256
    with open(file, 'rb') as f:
        return sha256(f.read()).hexdigest()


def read_file_flow(obj: str or list, vcf_only: bool, re: str, include=True) -> Iterable[IO]:
    """
    read
    :param obj:
    :param vcf_only:
    :param re:
    :param include:
    :return:
    """
    if type(obj) == list:
        for file in obj:
            yield from read_file_flow(file, vcf_only, re, include)
    else:
        if os.path.isfile(obj):
            if not ((search('vcf', os.path.basename(obj)) is None) and vcf_only):
                if (search(re, os.path.basename(obj)) is None) ^ include:
                    with open(obj, 'rb') as f:
                        yield f.read()
        else:
            file_list = [os.path.join(obj, file) for file in os.listdir(obj)]
            yield from read_file_flow(file_list, vcf_only, re, include)


def cal_multi_sha(file_list: str or list, vcf_only=True, re=r'', include=True):
    from hashlib import sha256
    value = sha256(b'')
    for file_flow in read_file_flow(file_list, vcf_only, re, include):
        value.update(file_flow)
    return value.hexdigest()


def extract_snp(vcf_file: str, temp_dir: str) -> None:
    """
    Extract snp using plink2
    :param vcf_file: Vcf file needed to extract
    :param temp_dir: Directory where the file is stored
    :return: Extracted vcf file
    """
    nopen = open_gz(vcf_file)
    need_lines = []
    need_snp_list = snp_list_rw(temp_dir, 'r')
    logging.info(f"Step 1: Get {vcf_file} index")
    idi_file = get_ref_index(vcf_file)
    logging.info(f"Step 2: Get {vcf_file} need snps' rows' number")
    with open(idi_file) as f:
        for line in f:
            temp = line.strip().split('\t')
            if temp[0] in need_snp_list:
                need_lines.append(int(temp[1]))
                need_snp_list.remove(temp[0])
    logging.info(f"Step 3: Extract {vcf_file} need snps")
    with open(os.path.join(temp_dir, os.path.basename(vcf_file)), 'wb') as fw:
        with nopen(vcf_file, 'rb') as f:
            i_line = 1
            for line in f:
                if not need_lines:
                    break
                if i_line == need_lines[0]:
                    fw.write(line)
                    need_lines.pop(0)
                i_line += 1
    logging.info(f"{vcf_file} finish extract.")


def merge_files(file_dir: str, out: str, valid=None, header=None) -> None:
    with open(out, 'wb') as fw:
        if valid:
            fw.write(b'##' + valid.encode() + b'\n')
        if header:
            fw.write(header)
        for file in os.listdir(file_dir):
            with open(os.path.join(file_dir, file), 'rb') as fr:
                for line in fr:
                    if line:
                        fw.write(line)
        if valid:
            fw.seek(2)
            fw.write(b'DONE')


def get_vcf_header(vcf_list: list, binary=False) -> str or bytes:
    for vcf_file in vcf_list:
        nopen = open_gz(vcf_file)
        if nopen:
            with nopen(vcf_file, 'rb') as f:
                for line in f:
                    if pat_header2.search(line):
                        if not pat_header1.search(line):
                            return line if binary else line.decode()


def get_gt_iid(vcf_list: list) -> pd.DataFrame:
    header = get_vcf_header(vcf_list)
    return pd.DataFrame(columns=[iid for iid in header.strip().split('\t')[9:]])


@use_time('Get reference GT data')
def get_ref_gt(vcf_list: list, data_dir: List[str], output: str, filter_vcf: str = 'filter.vcf') -> pd.DataFrame:
    if os.path.isfile(filter_vcf):
        if not verify_data(filter_vcf, vcf_list, 'multi'):
            get_ref_vcf(vcf_list, get_snp_list(data_dir))
    else:
        get_ref_vcf(vcf_list, get_snp_list(data_dir), output)
    gt_df = get_gt_iid(vcf_list)
    with open(os.path.join(output, 'population_QC', filter_vcf), 'rb') as f:
        f.readline()
        for line in f:
            temp = line.strip().split(b'\t')
            try:
                gt_df.loc[temp[2].decode()] = [get_gt(gt.decode(), temp[3].decode(), temp[4].decode(), warning=False)
                                               for gt in temp[9:]]
            except IndexError:
                raise
    return gt_df


def snp_list_rw(tempdir: str, method: str = 'w', need_snp_list: None or set or list = None):
    with open(os.path.join(tempdir, 'need_snp_list'), method + 'b' if method != 'raw' else 'w') as f:
        if method == 'w':
            dump(need_snp_list, f)
        elif method == 'r':
            return load(f)
        elif method == 'raw':
            for snp in need_snp_list:
                f.write(snp + '\n')


def rm_dir(rmdir: str) -> None:
    # Danger function, only can delete one level directory
    if os.listdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
    os.rmdir(rmdir)


def copy_files(files: list, from_dir: str, to_dir: str) -> None:
    for file in files:
        copy(os.path.join(from_dir, file), to_dir)


@use_time('Get reference subset vcf')
def get_ref_vcf(vcf_files: list, need_snp_list: set or list, output: str, out_name: str = 'filter.vcf',
                write_valid: bool = True, write_header: bool = False) -> None:
    temp_dir = mkdtemp(suffix='vcf')
    snp_list_rw(tempdir=temp_dir, need_snp_list=need_snp_list)  # temp
    func = partial(extract_snp, temp_dir=temp_dir)
    pro_num = len(vcf_files) if len(vcf_files) < os.cpu_count() else os.cpu_count()
    pool = Pool(processes=pro_num)
    pool.map(func, vcf_files)
    pool.close()
    pool.join()
    os.remove(os.path.join(temp_dir, 'need_snp_list'))
    valid = cal_multi_sha(vcf_files) if write_valid else None
    header = get_vcf_header(vcf_files, True) if write_header else None
    merge_files(temp_dir, os.path.join(output, 'population_QC', out_name), valid=valid, header=header)
    rm_dir(temp_dir)


def get_ref_freq(data_dir: List[str], vcf_files: list, output: str) -> None:
    if not os.path.isfile(os.path.join(output, 'population_QC', 'prs.ref.afreq')):
        get_ref_vcf(vcf_files, get_snp_list(data_dir, cal_ind=False, prs_ind=True), output, out_name='prs.vcf',
                    write_valid=False, write_header=True)
        run_plink_cmd(f"--vcf {os.path.join(output, 'population_QC', 'prs.vcf')} --freq "
                      f"--out {os.path.join(output, 'population_QC', 'prs.ref')}",
                      plink='plink2')
        os.remove(os.path.join(output, 'population_QC', 'prs.vcf'))


def verify_data(data_file: str, source_files: str or list, mode: str) -> bool:
    def_config = convert_dict(load_config()['file'], False)
    sha_cal = cal_sha if mode == 'single' else \
        partial(cal_multi_sha, vcf_only=False,
                re=def_config['need_suf'].replace('.', '\\.') + "|" +
                   def_config['gwas_suf'].replace('.', '\\.') + '|setting|\\.vcf',
                include=True)
    with open(data_file) as f:
        f.read(2)
        if f.read(4) == 'DONE':
            if f.read(61) == sha_cal(source_files)[3:]:
                return True
    return False


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


def same_status(files: list) -> bool:
    row = None
    for file in files:
        with open(file, 'rb') as f:
            if row:
                if row != f.readline():
                    return False
            else:
                row = f.readline()
    return True


def trans_gz(gz_file: str, out_dir: str) -> str:
    if gz_file.endswith('.gz'):
        with gzip_open(gz_file) as fr:
            with open(os.path.join(out_dir, 'prs_data'), 'w') as fw:
                for line in fr:
                    fw.write(line.decode())
        return os.path.join(out_dir, 'prs_data')
    else:
        return gz_file


def get_ref_index(vcf_file: str) -> str:
    id_index = vcf_file + '.idi'
    if os.path.isfile(id_index):
        if verify_data(id_index, vcf_file, 'single'):
            return id_index
    with open(id_index, 'w') as fw:
        nopen = open_gz(vcf_file)
        fw.write(f'###{cal_sha(vcf_file)}\n')
        with nopen(vcf_file, 'rb') as fr:
            global pattern
            i = 0
            for line in fr:
                i += 1
                if pattern.search(line):
                    fw.write(f'{pattern.search(line).group().decode()}\t{i}\n')
        fw.seek(2)
        fw.write('DONE')
    return id_index


def xlim_adjust(mean: float, x_range: tuple) -> Tuple[float, float]:
    left, right = (mean - x_range[0], x_range[1] - mean)
    return (x_range[0], mean + left) if left > right else (mean - right, x_range[1])


def get_type_list() -> Dict[str, str]:
    qual_dir = functions_config['file']['qual_data']
    quan_dir = functions_config['file']['quan_data']
    type_dict = dict()
    if qual_dir:
        type_dict.update({os.path.join(qual_dir, directory): 'qual' for directory in os.listdir(qual_dir)
                          if directory != '.DS_Store'})
    if quan_dir:
        type_dict.update({os.path.join(quan_dir, directory): 'quan' for directory in os.listdir(quan_dir)
                          if directory != '.DS_Store'})
    return type_dict


def get_subtype_list(type_list: Dict[str, str], need_type: str = 'quan' or 'qual') -> Dict[str, str]:
    return {directory: need_type for directory in type_list if type_list[directory] == need_type}


def get_gt(vcf_gt: str, ref: str, alt: str, snp: None or str = None, warning: bool = True) -> Set[str]:
    if '/' in vcf_gt:
        gt_n = vcf_gt.split('/')
    elif '|' in vcf_gt:
        gt_n = vcf_gt.split('|')
    else:
        gt_n = [vcf_gt]
    while '.' in gt_n:
        gt_n.remove('.')
    if '*' in gt_n:
        if warning:
            logging.warning(f'Find "*" in "{snp}", with which it may be difficult to compare the genotypes.')
    gt_set = set(alt.split(',')[int(i) - 1] if int(i) else ref for i in gt_n)
    return gt_set


@auto_encoding
def load_txt(description_text: str, encoding: str = 'UTF-8') -> Iterable[Tuple[str, str]]:
    with open(description_text, encoding=encoding) as f:
        for line in f:
            line = line.strip('\n\r')
            if line:
                try:
                    key, value = line.split(':')[0], ":".join(line.split(':')[1:])
                except ValueError:
                    raise Exception(f"Data format is bad in {description_text}")
                else:
                    value = value.strip(' ')
                    value = value.strip('"')
                    yield key, value


def get_cal_file() -> Tuple[bool, bool or None, str or None]:
    """
    Get the calculation file from work directory
    :return: Tuple[Whether exist, whether from GWAS, the name of the calculation file]
    """
    def_config = convert_dict(functions_config['file'], False)
    files = os.listdir()
    simple_files = list(filter(lambda a: a.endswith(def_config['need_suf']), files))
    gwas_files = list(filter(lambda a: a.endswith(def_config['gwas_suf']), files))
    gwas_gz_files = list(filter(lambda a: a.endswith(def_config['gwas_suf'] + '.gz'), files))
    if simple_files:
        return True, False, simple_files[0]
    elif gwas_files:
        return True, True, gwas_files[0]
    elif gwas_gz_files:
        return True, True, gwas_gz_files[0]
    else:
        return False, None, None


def get_snp_list(type_list: List[str], cal_ind: bool = True, prs_ind: bool = False) -> Set[str]:
    def_config = convert_dict(functions_config['name'], False)
    snp_list = set([])
    try:
        for data_dir in type_list:
            os.chdir(data_dir)
            for ind_dir in [item for item in sorted(os.listdir('.')) if os.path.isdir(item)]:
                os.chdir(ind_dir)
                exist, file_type, file = get_cal_file()
                if exist:
                    if cal_ind and not file_type:
                        snp_list = add_snp_list(file, snp_list, header=True)
                    elif prs_ind and file_type:
                        if file:
                            nopen = open_gz(file)
                            with nopen(file, 'rb') as f:
                                snp_list = load_prs_txt(f, snp_list, def_config)
                else:
                    logging.warning(f'Code({ind_dir}) has no corresponding algorithm!')
                os.chdir('..')
            os.chdir(raw_dir)
    finally:
        os.chdir(raw_dir)
    return snp_list


@auto_sep
def get_columns_num(file: str, need_columns: list, sep: str = '\t', start=1) -> List[int]:
    """
    Get corresponding columns' index from header
    :param file: Input file
    :param need_columns: The names of needed columns
    :param sep: The separator of file
    :param start: index start number
    :return: A list for corresponding columns' index
    """
    open_ = open_gz(file)
    with open_(file, 'rb') as f:
        header = f.readline()
    header = header.strip().split(sep.encode())
    assert len(header) > 1, sep.encode().join(header).decode()
    return [header.index(i.encode()) + start for i in need_columns]


def get_prs_res(result_file: str) -> float:
    with open(result_file, 'r') as f:
        temp = f.readlines()
    os.remove(result_file)
    return exp(float(temp[-1].strip().split('\t')[-1]))


def cal_trait(trait: str or None, gt: set, ref: str, inter: str, no_inter: str or None) -> str or None:
    lan_lib = convert_dict(functions_config['lan_lib'])
    if trait is None:
        if inter in lan_lib['if']:
            return 'Normal' if equal_gt(gt, trans_gt(ref)) else 'Normal' \
                if equal_gt(gt, trans_gt(gene_filp(ref))) else None
        else:
            return None
    else:
        if trait == 'Normal':
            if inter in lan_lib['if']:
                return 'Normal' if equal_gt(gt, trans_gt(ref)) else 'Normal' \
                    if equal_gt(gt, trans_gt(gene_filp(ref))) else None
            elif inter in lan_lib['and']:
                outcome = 'And_yes' if equal_gt(gt, trans_gt(ref)) else 'And_yes' \
                    if equal_gt(gt, trans_gt(gene_filp(ref))) else 'And_no'
                return outcome
            else:
                outcome = inter if equal_gt(gt, trans_gt(ref)) else inter \
                    if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
                return 'Normal' if outcome in lan_lib['none'] else outcome
        elif trait == 'And_yes':
            outcome = inter if equal_gt(gt, trans_gt(ref)) else inter \
                if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
            return 'Normal' if outcome in lan_lib['none'] else outcome
        elif trait == 'And_no':
            outcome = no_inter
            return 'Normal' if outcome in lan_lib['none'] else outcome
        else:
            return trait


# QC
# MAF
def qc_maf(failed_snps: set, file_path: str, report_dir: str, suffix: str = '') -> None:
    maf_threshold = eval(functions_config['parameter']['qc_maf'])
    data = pd.read_csv(os.path.join(file_path, 'qc.frq'), ' ', skipinitialspace=True)
    failed_snps.update(data.loc[data.MAF < maf_threshold, 'SNP'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.MAF, bins=50)
    plt.text(assign_pos(maf_threshold, plt.axis(), right_offset=0.02, left_offset=0.24), plt.axis()[3] * 0.184,
             f'MAF threshold: {maf_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(maf_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Minor allele frequency")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_MAF{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# SNPs_MISS
def qc_vmiss(failed_snps: set, file_path: str, report_dir, suffix: str = '') -> None:
    snp_miss_threshold = eval(functions_config['parameter']['qc_vmiss'])
    data = pd.read_csv(os.path.join(file_path, 'qc.lmiss'), ' ', skipinitialspace=True)
    failed_snps.update(data.loc[data.F_MISS > snp_miss_threshold, 'SNP'])
    plt.figure(dpi=300)
    if data.F_MISS.max() / snp_miss_threshold < 1e-1:
        plt.hist(data.F_MISS, range=[0, snp_miss_threshold], bins=25)
    else:
        plt.hist(data.F_MISS, bins=25)
    plt.text(assign_pos(snp_miss_threshold, plt.axis(), right_offset=0.03, left_offset=0.40), plt.axis()[3] * 0.184,
             f'Missing genotype rate\nper samplethreshold: {snp_miss_threshold}',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(snp_miss_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Missing genotype rate per sample")
    plt.ylabel("Number of samples")
    plt.savefig(os.path.join(report_dir, f'QC_var_miss{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# FMISS
def qc_smiss(failed_samples: set, file_path: str, report_dir, suffix: str = '') -> None:
    individual_miss_threshold = eval(functions_config['parameter']['qc_smiss'])
    data = pd.read_csv(os.path.join(file_path, 'qc.imiss'), ' ', skipinitialspace=True)
    failed_samples.update(data.loc[data.F_MISS > individual_miss_threshold, 'IID'])
    plt.figure(dpi=300)
    if data.F_MISS.max() / individual_miss_threshold < 1e-1:
        plt.hist(data.F_MISS, range=[0, individual_miss_threshold], bins=25)
    else:
        plt.hist(data.F_MISS, bins=25)
    plt.text(assign_pos(individual_miss_threshold, plt.axis(), right_offset=0.03, left_offset=0.40),
             plt.axis()[3] * 0.184,
             f'Missing genotype rate\nper variant threshold: {individual_miss_threshold}',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(individual_miss_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Missing genotype rate per variant")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_ind_miss{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# HARDY
def qc_hardy(failed_snps: set, file_path: str, report_dir, suffix: str = '') -> None:
    hardy_threshold = eval(functions_config['parameter']['qc_hardy'])
    data = pd.read_csv(os.path.join(file_path, 'qc.hwe'), ' ', skipinitialspace=True)
    failed_snps.update(data.loc[data.P < 10 ** -hardy_threshold, 'SNP'])
    data.loc[data.P == 0, 'P'] = 1e-295
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(-log10(data.P), bins=50)
    plt.text(assign_pos(hardy_threshold, plt.axis(), right_offset=0.02, left_offset=0.285), plt.axis()[3] * 0.816,
             f'HWE exact test p-value\nthreshold: 1e-{hardy_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(hardy_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("-log$_{10}$(HWE exact test p-value)")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_hardy{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# HET
def qc_het(failed_samples: set, file_path: str, report_dir, suffix: str = '') -> None:
    het_threshold = eval(functions_config['parameter']['qc_het'])
    data = pd.read_csv(os.path.join(file_path, 'qc.het'), ' ', skipinitialspace=True)
    het_mean, het_sd = data.F.mean(), data.F.std()
    failed_samples.update(data.loc[(data.F - het_mean > het_sd * het_threshold) |
                                   (data.F - het_mean < -het_sd * het_threshold), 'IID'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.F, bins=50)
    plt.text(assign_pos(het_mean - het_sd * het_threshold, plt.axis(), right_offset=0.02, left_offset=0.285),
             plt.axis()[3] * 0.816,
             f'F threshold:\n mean ±{het_threshold} sd',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(het_mean + het_sd * het_threshold, c='r', lw=1.5, ls='--')
    plt.axvline(het_mean - het_sd * het_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("F value")
    plt.ylabel("Number of individuals")
    plt.xlim(*xlim_adjust(het_mean, plt.xlim()))
    plt.savefig(os.path.join(report_dir, f'QC_het{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# SEX
def qc_sex(failed_samples: set, file_path: str, report_dir, suffix: str = '') -> None:
    male_threshold = eval(functions_config['parameter']['qc_male_F'])
    female_threshold = eval(functions_config['parameter']['qc_female_F'])
    data = pd.read_csv(os.path.join(file_path, 'qc.sexcheck'), ' ', skipinitialspace=True)
    failed_samples.update(data.loc[(data.F > male_threshold) & (data.F < female_threshold), 'IID'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.F, bins=50)
    plt.text(male_threshold - (plt.axis()[1] - plt.axis()[0]) * 0.17, plt.axis()[3] * 0.5,
             f'F threshold for \nmale: {male_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.text(female_threshold + (plt.axis()[1] - plt.axis()[0]) * 0.02, plt.axis()[3] * 0.816,
             f'F threshold for \nfemale: {female_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(male_threshold, c='r', lw=1.5, ls='--')
    plt.axvline(female_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("F value")
    plt.ylabel("Number of individuals")
    plt.savefig(os.path.join(report_dir, f'QC_sex{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


# RELATEDNESS
def qc_relatedness(failed_samples: set, file_path: str, report_dir, suffix: str = '') -> None:
    pihat_threshold = eval(functions_config['parameter']['qc_pihat'])
    data = pd.read_csv(os.path.join(file_path, 'qc.genome'), ' ', skipinitialspace=True)
    failed_samples.update(data.loc[data.PI_HAT > pihat_threshold, 'IID1'])
    failed_samples.update(data.loc[data.PI_HAT > pihat_threshold, 'IID2'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.PI_HAT, bins=50)
    plt.text(assign_pos(pihat_threshold, plt.axis(), right_offset=0.03, left_offset=0.33), plt.axis()[3] * 0.816,
             f'Pihat threshold:\n {pihat_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(pihat_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Pihat")
    plt.ylabel("Number of individual-individual pairs")
    plt.savefig(os.path.join(report_dir, f'QC_relatedness{"_" + suffix if suffix else ""}.png'), bbox_inches='tight')
    plt.close()


def vep_query(search_ids: list) -> Dict[str, dict]:
    import json
    import requests
    from urllib3.util.retry import Retry
    from requests.adapters import HTTPAdapter
    vep_api = "https://rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = json.dumps({'ids': search_ids, 'variant_class': '1'})
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount('https://', adapter)
    result = {}

    r = session.post(vep_api, headers=headers, data=data)
    if r.ok:
        decoded = r.json()
        result.update({res['id']: {'Variant type': res['most_severe_consequence'].replace('_', ' ').title(),
                                   'Variant class': res['variant_class'].replace('_', ' ')}
                       for res in decoded})
        r.close()
    else:
        r.raise_for_status()
    return result


@progress_value(10, average=True)
def vep_multiprocess(search_ids: list) -> List[dict]:
    queue = [search_ids[i: i + 100] for i in range(0, len(search_ids), 100)]
    process_num = os.cpu_count() if len(queue) > os.cpu_count() else len(queue)
    pool = Pool(processes=process_num)
    res = pool.map(vep_query, queue)
    pool.close()
    pool.join()
    return res


def load_clinvar(clinvar_vcf: str, clinvar_data: dict, human: Human) -> None:
    with gzip_open(clinvar_vcf) as f:
        for line in f:
            if line[0:1] != b'#':
                line = line.decode().strip().split('\t')
                rs_id = search(pattern_clinvar_rs, line[7])
                if rs_id:
                    rs = 'rs' + rs_id[0]
                    if rs in human.gt_data:
                        allele = line[4]
                        allele_type = search(pattern_clinvar, line[7])
                        if allele_type:
                            allele_type = allele_type[0].replace('_', ' ').capitalize()
                            if allele_type != 'Not provided':
                                clinvar_data[rs] = \
                                    {'Allele': allele, 'Allele type': allele_type, 'Clinvar ID': line[2],
                                     'Genotype': trans_gt(human.gt_data[rs], connector='')}


def recognize_color_clinvar(row: pd.Series) -> str:
    risk = True if 'athogenicity' in row['Allele type'] else True if 'Risk' in row['Allele type'] else True \
        if 'Drug response' in row['Allele type'] else None
    return recognize_color(row, risk, 'Allele', 'Genotype', na_value=None)


def recognize_color_pharmgkb(row: pd.Series) -> str:
    risk = True if 'increased' in row['Sentence'] else False if 'decreased' in row['Sentence'] else None \
        if 'not' in row['Sentence'] else True
    return recognize_color(row, risk, 'Alleles', 'Genotype')


def recognize_color(row: pd.Series, risk: bool, allele_col: str, genotype_col: str, na_value=None) -> str:
    if risk is None:
        return na_value
    else:
        if '+' not in row[allele_col]:
            if row[allele_col] in row[genotype_col]:
                return 'Red' if risk else 'Green'
        else:
            genotypes = row[allele_col].split('+')
            for genotype in genotypes:
                genotype = genotype.strip()
                if genotype in row[genotype_col]:
                    return 'Red' if risk else 'Green'
        return 'Green' if risk else 'Red'


def load_pharm(pharm_data: str, human: Human) -> pd.DataFrame:
    pharm_data = pd.read_csv(pharm_data, sep='\t', error_bad_lines=False)
    res = pharm_data.loc[pharm_data.Variant.isin(human.gt_data) & (pharm_data.Significance == 'yes'),
                         ['Variant', 'Alleles', 'Sentence', 'PMID']]
    res['Genotype'] = res.apply(lambda a: trans_gt(human.gt_data[a['Variant']], connector=''), axis=1)
    res['Color'] = res.apply(recognize_color_pharmgkb, axis=1)
    return res.dropna().sort_values(['Color', 'Variant'], ascending=[False, True]).reset_index(drop=True)


@progress_value(10, average=True)
def trans_clinvar_res(clinvar_res: dict) -> pd.DataFrame:
    res = pd.DataFrame(clinvar_res).T
    res.reset_index(inplace=True)
    res.columns = ['Variant'] + list(res.columns)[1:]
    res = res.loc[res.Allele != '.',]
    res['Color'] = res.apply(recognize_color_clinvar, axis=1)
    return res.dropna().sort_values(['Color', 'Variant'], ascending=[False, True]).reset_index(drop=True)


def get_clinvar_data(data_dir: str) -> None:
    # download data
    if not os.path.isfile(os.path.join(data_dir, 'clinvar.vcf.gz')):
        urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz',
                    os.path.join(data_dir, 'clinvar.vcf.gz'))

    if not os.path.isfile(os.path.join(data_dir, 'clinvar_papu.vcf.gz')):
        urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_papu.vcf.gz',
                    os.path.join(data_dir, 'clinvar_papu.vcf.gz'))


@progress_value(10, average=True)
def get_pharmgkb_data(data_dir: str) -> None:
    if not os.path.isfile(os.path.join(data_dir, 'pharmgkb.zip')):
        urlretrieve('https://api.pharmgkb.org/v1/download/file/data/annotations.zip',
                    os.path.join(data_dir, 'pharmgkb.zip'))
    with ZipFile(os.path.join(data_dir, 'pharmgkb.zip')) as fzip:
        if not os.path.isfile(os.path.join(data_dir, 'var_pheno_ann.tsv')):
            fzip.extract('var_pheno_ann.tsv', data_dir)
        if not os.path.isfile(os.path.join(data_dir, 'var_fa_ann.tsv')):
            fzip.extract('var_fa_ann.tsv', data_dir)
        if not os.path.isfile(os.path.join(data_dir, 'var_drug_ann.tsv')):
            fzip.extract('var_drug_ann.tsv', data_dir)


def write_failed_qc_list(failed_snps: set, failed_samples: set, output: str) -> None:
    with open(os.path.join(output, 'population_QC', "failed_qc_samples"), 'w') as f:
        for sample in failed_samples:
            f.write(sample + ' ' + sample + '\n')
    with open(os.path.join(output, 'population_QC', "failed_qc_snps"), 'w') as f:
        for snp in failed_snps:
            f.write(snp + '\n')


@use_time('Get reference result data')
def get_ref_cal(ref_file_list: List[str], type_list: dict, failed_snps: set, failed_samples: set,
                output: str, suffix: str = '') -> Tuple[pd.DataFrame, pd.DataFrame]:
    ref_code_file = os.path.join(output, 'population_QC', f'population_code_res{"_" + suffix if suffix else ""}.csv')
    ref_rs_file = os.path.join(output, 'population_QC', f'population_rs_res{"_" + suffix if suffix else ""}.csv')
    ref_config = os.path.join(output, 'population_QC', 'setting')
    def_config = convert_dict(functions_config['parameter'])
    def_config.update(convert_dict(functions_config['file'], False))
    def_config = {key: def_config[key] for key in def_config if key in ['ref_structure', 'snp', 'ea', 'p',
                                                                        'beta', 'p_threshold', 'clump-p1', 'clump-r2',
                                                                        'clump-kb', 'clump-p2', 'description_suf',
                                                                        'need_suf',
                                                                        'gwas_suf', 'use_qc_ref']}
    use_ref_qc = def_config['use_qc_ref']
    with open(ref_config, 'wb') as f:
        dump(def_config, f)
    if os.path.isfile(ref_code_file) and os.path.isfile(ref_rs_file):
        if verify_data(ref_code_file, ref_file_list + list(type_list) + [ref_config], 'multi'):
            if same_status([ref_code_file, ref_rs_file]):
                return pd.read_csv(ref_rs_file, skiprows=1, index_col=0), \
                       pd.read_csv(ref_code_file, skiprows=1, index_col=0)
    outcomes_rs = pd.DataFrame()
    outcomes_code = pd.DataFrame()
    gt_data = get_ref_gt(ref_file_list, list(type_list), output)
    if use_ref_qc:
        gt_data = gt_data.loc[~gt_data.index.isin(failed_snps), ~gt_data.columns.isin(failed_samples)]
        write_failed_qc_list(failed_snps, failed_samples, output)
    with open(ref_code_file, 'w', newline='', encoding='UTF-8') as fc:
        with open(ref_rs_file, 'w', newline='', encoding='UTF-8') as fr:
            sha = cal_multi_sha(ref_file_list + list(type_list) + [ref_config], vcf_only=False,
                                re=def_config['need_suf'].replace('.', '\\.') + "|" +
                                   def_config['gwas_suf'].replace('.', '\\.') + '|setting|\\.vcf',
                                include=True)
            fc.write(f'###{sha}\n')
            fr.write(f'###{sha}\n')
            try:
                for type_dir in tqdm(list(type_list)):
                    os.chdir(type_dir)
                    for ind_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
                        os.chdir(ind_dir)
                        exist, file_type, file = get_cal_file()
                        if exist:
                            if file_type:
                                res = get_ref_prs(file, ref_file_list, use_ref_qc, output)
                                if outcomes_code.empty:
                                    outcomes_code[ind_dir] = res.apply(lambda x: exp(x))
                                else:
                                    res = res.rename(ind_dir)
                                    outcomes_code = outcomes_code.merge(res.apply(lambda x: exp(x)), left_index=True,
                                                                        right_index=True)
                            else:
                                outcomes_rs, outcomes_code = cal_ref_from_txt(file, type_list[type_dir], ind_dir,
                                                                              gt_data, outcomes_rs, outcomes_code,
                                                                              header=True)
                        os.chdir('..')
                    os.chdir(raw_dir)
            finally:
                os.chdir(raw_dir)
                outcomes_rs.to_csv(fr)
                fr.seek(2)
                fr.write('DONE')
            outcomes_code.to_csv(fc)
            fc.seek(2)
            fc.write('DONE')
    return outcomes_rs, outcomes_code


@progress_value(5)
@use_time()
def recode_and_sex_impute(human: Human, file: str, temp_dir: str):
    """
    recode the sample, determine the sex of the sample
    :param human: Objects Human
    :param file: the sample path
    :param temp_dir: the temp temporary directory path
    :return: a converted sample file in the temporary directory and a sex result
    """
    file_type, file_name = file_recognize(file)
    try:
        run_plink_cmd(f"--{file_type} {file_name} "
                      "--impute-sex y-only "  # Should be reconsider
                      "--recode vcf "
                      f"--out {os.path.join(temp_dir, 'temp')}",
                      plink='plink')
        data = pd.read_csv(os.path.join(temp_dir, 'temp.sexcheck'), ' ', skipinitialspace=True)
    except Exception as e:
        logging.warning('Sex impute failed:\n' + f'{e.args}')
        run_plink_cmd(f"--{file_type} {file_name} "
                      "--recode vcf "
                      f"--out {os.path.join(temp_dir, 'temp')}",
                      plink='plink')
        human.vcf = os.path.join(temp_dir, 'temp.vcf')
    else:
        run_plink_cmd(f"--vcf {os.path.join(temp_dir, 'temp.vcf')} "
                      "--recode vcf "
                      "--rm-dup force-first "
                      f"--out {os.path.join(temp_dir, 'temp')}")
        if file_type == 'vcf':
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.SNPSEX[0] == 1 else 'Female'
        else:
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.PEDSEX[0] == 1 else 'Female'


@use_time('Get population prs result')
def get_ref_prs(prs_data: str, vcf_files: list, qc_ref: bool, output: str) -> pd.DataFrame:
    """
    Get prs result from population for single indicator
    :param prs_data: Information for PRS SNPs
    :param vcf_files: Population genotype data
    :param qc_ref: Whether include QC-failed samples in population data
    :param output: The output directory for whole analysis
    :return: A dataframe which stores prs score result from population for single indicator
    """
    def_config = convert_dict(functions_config['parameter'])
    def_config.update(convert_dict(functions_config['file'], False))
    def_config.update(convert_dict(functions_config['name'], False))
    ref_structure = os.path.realpath(os.path.join(raw_dir, def_config['ref_structure']))
    temp_dir = mkdtemp(suffix='prs')
    save_code = os.path.basename(prs_data).split('.')[0]
    # Clump
    clump_p2 = def_config["clump-p2"] if def_config["clump-p1"] < def_config["clump-p2"] \
        else def_config["clump-p1"]
    run_plink_cmd(f'--vcf {ref_structure} --clump-p1 {str(def_config["clump-p1"])} '
                  f'--clump-p2 {str(clump_p2)} '
                  f'--clump-r2 {str(def_config["clump-r2"])} --clump-kb {str(def_config["clump-kb"])} '
                  f'--clump {prs_data} --clump-snp-field {def_config["snp"]} '
                  f'--clump-field {def_config["p"]} --out {os.path.join(temp_dir, save_code)}', plink='plink')
    # Produce plink need data
    select_columns(os.path.join(temp_dir, save_code + '.clumped'), [2],
                   os.path.join(temp_dir, 'valid.snp'), header=False, skipspace=True, sep=' ')
    select_columns(prs_data, [0, 4], os.path.join(output, 'population_QC', f'{save_code}.SNP.pvalue'))
    with open(os.path.join(output, 'population_QC', 'need_range_list'), 'w') as f:
        f.write(f'{def_config["p_threshold"]} 0 {def_config["p_threshold"]}\n')
    # PRS calculate
    prs_data = get_score_data(prs_data, def_config, temp_dir)
    qc_str = f"--remove {os.path.join(output, 'population_QC', 'failed_qc_samples')} " \
             f"--exclude {os.path.join(output, 'population_QC', 'failed_qc_snps')} " if qc_ref else ""
    for fid, file in enumerate(vcf_files):
        run_plink_cmd(f'--vcf {file} {qc_str}--score {prs_data} 1 2 3 header '
                      f'--q-score-range {os.path.join(output, "population_QC", "need_range_list")} '
                      f'{os.path.join(output, "population_QC", save_code + ".SNP.pvalue")} '
                      f'--extract {os.path.join(temp_dir, "valid.snp")} --out {os.path.join(temp_dir, "result_prs")}',
                      plink='plink')
        if fid == 0:
            res_data = pd.read_csv(os.path.join(temp_dir, f"result_prs.{def_config['p_threshold']}.profile"), ' ',
                                   skipinitialspace=True, index_col=1)['SCORE']
        else:
            res_data = res_data + pd.read_csv(os.path.join(temp_dir,
                                                           f"result_prs.{def_config['p_threshold']}.profile"),
                                              ' ', skipinitialspace=True, index_col=1)['SCORE']
    rm_dir(temp_dir)
    return res_data


def risk_cal(string: str, risk_score: float) -> str:
    string = string.split('(')[0].strip()
    num = 0
    comma = False
    for i in range(len(string)):
        if string[i] == ',':
            comma = True
        elif not string[i].isnumeric() and string[i] != '.':
            break
        num += 1
    risk = float(string[:num].replace(',', '')) * risk_score
    return f'{risk:.2f}{string[num:]}' if not comma else f'{risk:,.2f}{string[num:]}'


def read_line_num(gwas_file: str, header: bool = False) -> int:
    nopen = open_gz(gwas_file)
    line_num = -1 if header else 0
    with nopen(gwas_file, 'rb') as f:
        for _ in f:
            line_num += 1
    return line_num


def chromo_info(human: Human, temp_dir: str, img_dir: str) -> None:
    run_plink_cmd(f'--vcf {human.vcf} --make-just-bim --sort-vars --out {os.path.join(temp_dir, "sample_qc")}')
    data = pd.read_csv(os.path.join(temp_dir, "sample_qc.bim"),
                       '\t', skipinitialspace=True, header=None, dtype={0: str})
    data.columns = ['Chromosome', 'ID', 'Unknown', 'Position', 'Ref', 'ALT']

    pos_count = data.Chromosome.value_counts(sort=False)
    pos_count.sort_index(key=lambda a: [chrom_key(i) for i in a], inplace=True)
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.bar(pos_count.index, pos_count)
    plt.grid(axis='y')
    plt.xlabel("Chromosome")
    plt.ylabel("Variations count")
    plt.tight_layout()
    plt.savefig(os.path.join(img_dir, f'chromosome_count.png'))
    plt.close()

    # if convert_dict(functions_config['module'])['vep']:
    #     vep_res = vep_multiprocess(list(data.ID))
    #     vep_res = {res: queue[res] for queue in vep_res for res in queue}
    #     vep_res = pd.DataFrame(vep_res).T
    #     vep_result = pd.merge(data.loc[:, 'ID'], vep_res, right_index=True, left_on='ID', how='left')
    #     vep_class_result = vep_result['Variant class'].fillna('Unknown').value_counts()
    #     vep_type_result = vep_result['Variant type'].fillna('Unknown').value_counts()
    #     arrow_class = False if vep_result['Variant class'].value_counts(normalize=True)[0] > 0.80 else True
    #     arrow_type = False if vep_result['Variant type'].value_counts(normalize=True)[0] > 0.80 else True
    #     pie_plot(vep_class_result, vep_class_result.index, os.path.join(img_dir, './vep_class_res.png'),
    #              arrow=arrow_class)
    #     pie_plot(vep_type_result, vep_type_result.index, os.path.join(img_dir, './vep_type_res.png'), arrow=arrow_type,
    #              legend_ratio=2)


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


def minor_freq_cal(i_dat: str or float) -> None or float:
    if type(i_dat) == str:
        if ',' in i_dat:
            return None
        else:
            minor = float(i_dat)
            return minor if minor <= 0.5 else 1 - minor
    elif type(i_dat) == float:
        return i_dat


def qc_sample_maf(freq_data: pd.DataFrame, report_dir: str) -> None:
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(freq_data['ALT_FREQS'], bins=50)
    plt.xlabel("Minor allele frequency")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_SAMPLE_MAF.png'), bbox_inches='tight')
    plt.close()


def geno_diff(raw_human: Human, ref_human: Human, intersection: Optional[set] = None) -> int:
    same_num = 0
    if not intersection:
        intersection = set(raw_human.gt_data.keys()) & set(ref_human.gt_data.keys())
    for snps in intersection:
        if equal_gt(raw_human.gt_data[snps], ref_human.gt_data[snps]):
            same_num += 1
    return same_num


def get_maf(human: Human, temp_dir: str, maf_ref: str, img_dir: str) -> str:
    file_prefix, file_name = file_recognize(maf_ref)
    snp_list_rw(temp_dir, 'raw', human.snps)
    run_plink_cmd(f'--{file_prefix} {file_name} --extract {os.path.join(temp_dir, "need_snp_list")} '
                  f'--freq --out {os.path.join(temp_dir, "sample_qc")}', plink='plink2')
    freq_right = pd.read_csv(os.path.join(temp_dir, "sample_qc.afreq"), sep='\t', usecols=[1, 4])
    freq_right['ALT_FREQS'] = freq_right['ALT_FREQS'].apply(minor_freq_cal)
    freq_left = pd.DataFrame({'ID': human.snps})
    freq_data = pd.merge(freq_left, freq_right, how='left')
    freq_data['ALT_FREQS'] = freq_data['ALT_FREQS'].apply(lambda a: a if a <= 0.5 else 1 - a)
    qc_sample_maf(freq_data, img_dir)
    res = freq_data.count().to_list()
    return f'{res[1]} of {res[0]} variants were found in MAF reference data.'


def process_line(line: bytes) -> Tuple[np.ndarray, str]:
    line = line.strip().split(b'\t')
    geno_dat = line[6:]
    return np.array([int(geno.decode()) if geno != b'NA' else None for geno in geno_dat], dtype='f'), line[1].decode()


def add_res(geno_res: np.ndarray, iid: str, data: list, iid_list: list, pbar: tqdm) -> None:
    data.append(geno_res)
    iid_list.append(iid)
    pbar.update(1)


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
    if os.path.normcase(temp_dir) != os.path.normcase(save_dir):
        if prune:
            copy(os.path.join(temp_dir, 'prune_pop.prune.in'), save_dir)
        copy(os.path.join(temp_dir, 'ref_pca.eigenvec.allele'), save_dir)
        copy(os.path.join(temp_dir, 'ref_pca.eigenvec'), save_dir)
        copy(os.path.join(temp_dir, 'ref_pca.eigenval'), save_dir)
        copy(os.path.join(temp_dir, 'ref_pca.acount'), save_dir)


def convert_data(file_prefix: str, file_name: str, temp_dir: str, prune: bool, mode: str) -> None:
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
    # UMAP
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
    ref_data = pd.DataFrame([iid.split('_') for iid in iid_list])
    ref_data.columns = ['FID', 'IID'] if ref_data.shape[1] > 1 else ['IID']
    ref_data.to_csv(os.path.join(save_dir, "sample_info"), index=False)
    umap_fit = umap.UMAP(random_state=42)
    umap_fit.fit(dat)
    with gzip_open(os.path.join(save_dir, 'umap_fit'), 'wb', compresslevel=6) as f:
        f.write(dumps(umap_fit))


def get_metadata(metadata: Optional[str], ref_data: pd.DataFrame,
                 pop_id: Optional[str] = None, pop_col: Optional[str] = None) -> pd.DataFrame:
    if not pop_id:
        pop_id = functions_config['name']['population_id']
    if not pop_col:
        pop_col = functions_config['name']['population_col']
    if metadata:
        with open(metadata) as f:
            population_sep = detect(f.readline(10))
        population_data = pd.read_csv(metadata, sep=population_sep)
        return pd.merge(ref_data, population_data, left_on='IID', right_on=pop_id)
    else:
        ref_data.insert(0, pop_col, 'Reference')
        return ref_data


def pca_sample(sample_vcf: str, work_dir: str, temp_dir: str, prune: bool) -> Tuple[float, float]:
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


def pca_plot(sample_vcf: Optional[str], metadata: Optional[str], work_dir: str, temp_dir: str, img_dir: str,
             prune: bool, pop_id: Optional[str] = None, pop_col: Optional[str] = None) -> None:
    if not pop_col:
        pop_col = functions_config['name']['population_col']
    pca_data = pd.read_csv(os.path.join(work_dir, 'ref_pca.eigenvec'), sep='\t')
    if '#IID' in pca_data.columns:
        pca_data['IID'] = pca_data['#IID'].apply(lambda a: a.split('_')[-1])
    ref_pca_data = get_metadata(metadata, pca_data, pop_id, pop_col)
    sample_pc = pca_sample(sample_vcf, work_dir, temp_dir, prune) if sample_vcf else None
    draw_plot(ref_pca_data, sample_pc, pop_col, img_dir, 'PC')


def draw_plot(ref_data: pd.DataFrame, sample: Optional[Tuple[float, float]],
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


def sample_umap(sample_vcf: str, umap_fit: umap.UMAP, work_dir: str, temp_dir: str,
                prune: bool) -> Tuple[float, float]:
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


def umap_plot(sample_vcf: Optional[str], metadata: Optional[str], work_dir: str, temp_dir: str, img_dir: str,
              prune: bool, pop_id: Optional[str] = None, pop_col: Optional[str] = None) -> None:
    if not pop_col:
        pop_col = functions_config['name']['population_col']
    with gzip_open(os.path.join(work_dir, 'umap_fit')) as f:
        umap_fit = load(f)
    umap_data = pd.read_csv(os.path.join(work_dir, "sample_info"))
    umap_data = pd.concat([umap_data, pd.DataFrame(umap_fit.embedding_, columns=['UMAP1', 'UMAP2'])], axis=1)
    sample_value = sample_umap(sample_vcf, umap_fit, work_dir, temp_dir, prune) if sample_vcf else None
    ref_umap_data = get_metadata(metadata, umap_data, pop_id, pop_col)
    draw_plot(ref_umap_data, sample_value, pop_col, img_dir, 'UMAP')


def ps_plot(sample_vcf: Optional[str], metadata: Optional[str], work_dir: str, temp_dir: str, img_dir: str,
            prune: bool, pop_id: Optional[str] = None, pop_col: Optional[str] = None) -> None:
    pca_plot(sample_vcf, metadata, work_dir, temp_dir, img_dir, prune, pop_id, pop_col)
    umap_plot(sample_vcf, metadata, work_dir, temp_dir, img_dir, prune, pop_id, pop_col)


def ps_verify(ps_file: str, work_dir: str) -> bool:
    ps_save = os.path.join(work_dir, 'ps_data')
    if os.path.isfile(ps_save):
        if verify_data(ps_save, ps_file, 'single'):
            return True
    return False


def ps_analyse(ps_ref: str, sample_vcf: Optional[str], metadata: Optional[str],
               temp_dir: str, save_dir: str,
               img_dir: Optional[str] = None, prune: Optional[bool] = None, **kwargs) -> str:
    if not prune:
        prune = convert_dict(functions_config['parameter'])['ps_prune']
    if not img_dir:
        img_dir = save_dir
    if not ps_verify(ps_ref, save_dir):
        get_ps_need(ps_ref, save_dir, temp_dir, prune)
    if 'pop_id' not in kwargs:
        kwargs['pop_id'] = None
    if 'pop_col' not in kwargs:
        kwargs['pop_col'] = None
    ps_plot(sample_vcf, metadata, save_dir, temp_dir, img_dir,
            prune, pop_id=kwargs['pop_id'], pop_col=kwargs['pop_col'])
    return 'Yes'



def chrom_key(chromosome: str) -> int:
    key_dict = {'X': 23, 'Y': 24, 'MT': 25}
    if chromosome.isnumeric():
        return int(chromosome)
    elif chromosome in key_dict:
        return key_dict[chromosome]
    else:
        return 99


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


@auto_sep
@auto_encoding
def cal_ref_from_txt(file: str, type_: str, ind_dir: str,
                     gt_data: pd.DataFrame, outcomes_rs: pd.DataFrame, outcomes_code: pd.DataFrame,
                     *, encoding: str = 'UTF-8', sep: str = '\t',
                     header: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    columns = convert_dict(functions_config['columns'])
    with open(file, 'r', encoding=encoding) as f:
        first = True
        for line in f:
            if line.strip("\n\r"):
                line = line.strip("\n\r").split(sep)
                if first:
                    assert len(line) > 1, sep.join(line)
                    first = False
                if header:
                    columns_idx = get_columns_index(line, list(columns['columns_' + type_].values()))
                    if type_ == 'qual':
                        assert None not in columns_idx,\
                            f"Cannot find {list(columns['columns_qual'].keys())[columns_idx.index(None)]} " \
                            f"in the header of algorithm file ({file})."
                    elif type_ == 'quan':
                        assert None not in columns_idx[:-2], \
                            f"Cannot find {list(columns['columns_quan'].keys())[columns_idx[:-2].index(None)]} " \
                            f"in the header of algorithm file ({file})."
                        assert any(columns_idx[-2:]),\
                            f"Cannot find statistic (OR or Beta) in the header of algorithm file ({file})."
                        statistic = 'OR' if columns_idx[-2] else 'Beta'
                        columns_idx.pop(-1 if statistic == 'OR' else -2)
                    header = False
                else:
                    rs = line[columns_idx[0]]
                    if rs in gt_data.index:
                        if type_ == 'quan':
                            outcome = gt_data.loc[rs].apply(
                                lambda gt: cal_effect(gt, *select_list(line, columns_idx[1:3]), statistic))
                            code_rs = f'{ind_dir}_{rs}'
                            if code_rs not in outcomes_rs.columns:
                                outcomes_rs[code_rs] = outcome
                                if ind_dir in outcomes_code.columns:
                                    outcomes_code[ind_dir] = outcomes_code[ind_dir] * outcome
                                else:
                                    outcomes_code[ind_dir] = outcome
                        else:
                            if ind_dir not in outcomes_code.columns:
                                outcomes_code[ind_dir] = pd.Series('Normal',
                                                                   index=gt_data.columns)
                            temp = list(map(lambda trait, gt: cal_trait(trait, gt, *select_list(line, columns_idx[1:4])),
                                            outcomes_code[ind_dir], gt_data.loc[rs]))
                            outcomes_code[ind_dir] = temp
    return outcomes_rs, outcomes_code


@auto_encoding
@auto_sep
def add_snp_list(file: str, snp_list: Set[str] or List[str], *,
                 encoding: str = 'UTF-8', sep: str = '\t', header=False) -> Set[str] or List[str]:
    first = True
    with open(file, 'r', encoding=encoding) as f:
        for line in f:
            if first:
                assert len(line.strip('\n\r').split(sep)) > 1, line
                first = False
            if header:
                columns = convert_dict(functions_config['columns'])
                line = line.strip('\n\r').split(sep)
                columns_idx = filter(lambda a: a in columns['columns_quan']['SNP'] + columns['columns_qual']['SNP'],
                                     line)
                columns_idx = line.index(next(columns_idx))
                header = False
            else:
                snp_list.add(line.strip('\n\r').split(sep)[columns_idx])
    return snp_list


@auto_sep
def load_prs_txt(f: BinaryIO, snp_list: List[str] or Set[str], def_config: dict, *,
                 sep: str = '\t') -> List[str] or Set[str]:
    header = True
    snp_column = None
    for line in f:
        if header:
            if not snp_column:
                assert len(line.split(sep.encode())) > 1, line
                snp_column = line.split(sep.encode()).index(def_config['snp'].encode())
            header = False
            continue
        else:
            snp_list.add(line.split(sep.encode())[snp_column].decode())
    return snp_list


@auto_encoding
def load_other_database(database: str, human: Human) -> pd.DataFrame:
    with open(database) as f:
        sep = detect(f.readline(10), whitelist=['\t', ',', ' '], default='\t')
    data = pd.read_csv(database, sep=sep, dtype=object, error_bad_lines=False)
    res = data.loc[data[functions_config['name']['database_snp']].isin(human.gt_data), :]
    res.loc['Genotype'] = res.apply(lambda a: trans_gt(human.gt_data[a[functions_config['name']['database_snp']]],
                                                       connector=''), axis=1)
    return res.dropna().sort_values([functions_config['name']['database_snp']], ascending=[True]).reset_index(drop=True)


def ref_data_qc(temp_dir: str, vcf_list: list, img_dir: str, failed_snps: set, failed_samples: set,
                suffix: str = '') -> None:
    assert len(vcf_list) == 1, 'Reference QC can only execute with single vcf'
    run_plink_cmd(f'--vcf {vcf_list[0]} --hardy --het --check-sex --missing --freq --test-missing --genome '
                  f'--out {os.path.join(temp_dir, "qc")}', plink='plink')
    qc_maf(failed_snps, temp_dir, img_dir, suffix=suffix)
    qc_vmiss(failed_snps, temp_dir, img_dir, suffix=suffix)
    qc_smiss(failed_samples, temp_dir, img_dir, suffix=suffix)
    qc_het(failed_samples, temp_dir, img_dir, suffix=suffix)
    qc_hardy(failed_snps, temp_dir, img_dir, suffix=suffix)
    qc_sex(failed_samples, temp_dir, img_dir, suffix=suffix)
    qc_relatedness(failed_samples, temp_dir, img_dir, suffix=suffix)


# TODO: Add merge function
def merge_ref():
    pass


def get_vcf_files(directory: str) -> List[str]:
    return [os.path.abspath(f"{os.path.join(directory, file)}") for file in os.listdir(directory) if
            'vcf' in file and 'idi' not in file and 'tbi' not in file]


def get_ref_data_type() -> List[List[str]]:
    quan_dir = functions_config['file']['quan_ref']
    qual_dir = functions_config['file']['qual_ref']
    if quan_dir and qual_dir:
        if os.path.realpath(quan_dir) == os.path.realpath(qual_dir):
            return [get_vcf_files(quan_dir)]
        else:
            return list(map(get_vcf_files, [qual_dir, quan_dir]))
    elif qual_dir:
        return [get_vcf_files(qual_dir)]
    elif quan_dir:
        return [get_vcf_files(quan_dir)]
    else:
        raise FileNotFoundError('No reference data was found!')


def mkdir(directories: List[str], pre_dir: str = '') -> None:
    for directory in directories:
        if not os.path.isdir(os.path.join(pre_dir, directory)):
            os.mkdir(os.path.join(pre_dir, directory))


@auto_sep
def get_score_data(file: str, def_config, temp_dir: str, sep: str = '\t') -> str:
    score_file = os.path.join(temp_dir, os.path.basename(file).split('.')[0])
    if not os.path.isfile(score_file):
        open_ = open_gz(file)
        with open_(file, 'rb') as fr:
            header = fr.readline()
            header = header.strip().split(sep.encode())
            header = [col.decode() for col in header]
            assert len(header) > 1, sep.join(header)
            if def_config['beta'] in header:
                need_idx = [header.index(def_config[i]) for i in ['snp', 'ea', 'beta']]
                mode = 'beta'
            elif def_config['or'] in header:
                need_idx = [header.index(def_config[i]) for i in ['snp', 'ea', 'or']]
                mode = 'or'
                from math import log
            else:
                raise Exception(f"Cannot find statistic (Beta or OR) in the file ({file}).")
            with open(score_file, 'wb') as fw:
                for line in fr:
                    line = line.strip().split(sep.encode())
                    fw.write(b'\t'.join(select_list(line, need_idx[:-1])))
                    beta = line[need_idx[-1]]
                    beta = beta if mode == 'beta' else str(log(float(beta.decode()))).encode()
                    fw.write(b'\t' + beta + b'\n')
    return score_file


def add_header(result: dict) -> str or None:
    high_risk_items = []
    for item in result['Quantitative_traits']:
        if item['Distribution']['Low_p'] > 75:
            high_risk_items.append(item['Name'])
    if high_risk_items:
        res = '<b>High risk disease:</b><br>&nbsp;&nbsp;&nbsp;&nbsp;' + '<br>&nbsp;&nbsp;&nbsp;&nbsp;'.join(high_risk_items)
        return res
    return None


def find_samples_file(dir_: str) -> str:
    for file in os.listdir(dir_):
        if sample_file_pat.search(file):
            return os.path.join(dir_, file)
    return ''
