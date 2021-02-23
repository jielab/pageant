from math import exp
from pickle import dump, load
from multiprocessing import Pool
from functools import partial
from gzip import open as gzip_open
from re import compile, sub
from tempfile import mkdtemp
from urllib.request import urlretrieve
from zipfile import ZipFile
from tqdm import tqdm
from src.objects import *

pattern = compile(rb'(?<=[\t])rs[0-9]*(?=[\t;])')
pat_header1 = compile(rb'^##')
pat_header2 = compile(rb'^#')
pattern_clinvar = compile(r'(?<=CLNSIG=)\w+(?=[,;])')
pattern_clinvar_rs = compile(r'(?<=RS=)\d+')
pattern_load_num = compile(r'\d+(?= valid predictors loaded)')


# Sub_functions
def is_excel(file: str) -> bool:
    return True if file.split('.')[-1] == 'xlsx' or file.split('.')[-1] == 'xls' else False


def select_list(ob_list: list, index):
    # assert len(ob_list) < max(index), 'list index out of range'
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def sub_col(ob_list: list, index: OrderedDict, code=True):
    if not code:
        index.pop('code')
    return {key: ob_list[index[key]] if index[key] < len(ob_list) else None for key in index}


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


def read_file_flow(obj: str or list, vcf_only: bool, re: str, include=True):
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


# @use_time('Extract_snp from reference data')
def extract_snp(vcf_file: str, temp_dir: str):
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


def merge_files(file_dir: str, out: str, valid=None):
    with open(out, 'wb') as fw:
        if valid:
            fw.write(b'#' + valid.encode() + b'\n')
        for file in os.listdir(file_dir):
            with open(os.path.join(file_dir, file), 'rb') as fr:
                for line in fr:
                    if line:
                        fw.write(line)
        if valid:
            fw.seek(0)
            fw.write(b'DONE')


def get_gt_iid(vcf_list: list) -> pd.DataFrame:
    for vcf_file in vcf_list:
        nopen = open_gz(vcf_file)
        if nopen:
            with nopen(vcf_file, 'rb') as f:
                for line in f:
                    if pat_header2.search(line):
                        if not pat_header1.search(line):
                            return pd.DataFrame(columns=[iid.decode() for iid in line.strip().split(b'\t')[9:]])


@progress_value(5)
@use_time('Get reference GT data')
def get_ref_gt(vcf_list: list, data_dir: str, filter_vcf=os.path.join('ref_res', 'filter.vcf')) -> pd.DataFrame:
    if os.path.isfile(filter_vcf):
        if not verify_data(filter_vcf, vcf_list, 'multi'):
            get_ref_vcf(vcf_list, get_snp_list(data_dir))
    else:
        get_ref_vcf(vcf_list, get_snp_list(data_dir))
    gt_df = get_gt_iid(vcf_list)
    with open(filter_vcf, 'rb') as f:
        f.readline()
        for line in f:
            temp = line.strip().split(b'\t')
            try:
                gt_df.loc[temp[2].decode()] = [get_gt(gt.decode(), temp[3].decode(), temp[4].decode(), warning=False)
                                      for gt in temp[9:]]
            except IndexError:
                raise
    return gt_df


def snp_list_rw(tempdir: str, method='w', need_snp_list=None):
    with open(os.path.join(tempdir, 'need_snp_list'), method + 'b' if method != 'raw' else 'w') as f:
        if method == 'w':
            dump(need_snp_list, f)
        elif method == 'r':
            return load(f)
        elif method == 'raw':
            for snp in need_snp_list:
                f.write(snp + '\n')


def rm_dir(rmdir: str):
    # danger function, only can delete one level directory
    if os.listdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
    os.rmdir(rmdir)


def copy_files(files: list, from_dir: str, to_dir: str):
    for file in files:
        copy(os.path.join(from_dir, file), to_dir)


@use_time('Get reference subset vcf')
def get_ref_vcf(vcf_files: list, need_snp_list: set or list, out_name='filter.vcf'):
    temp_dir = mkdtemp(suffix='vcf')
    snp_list_rw(tempdir=temp_dir, need_snp_list=need_snp_list)  # temp
    func = partial(extract_snp, temp_dir=temp_dir)
    pro_num = len(vcf_files) if len(vcf_files) < os.cpu_count() else os.cpu_count()
    pool = Pool(processes=pro_num)  # todo: cpu_count need warp
    pool.map(func, vcf_files)
    pool.close()
    pool.join()
    os.remove(os.path.join(temp_dir, 'need_snp_list'))
    merge_files(temp_dir, os.path.join('ref_res', out_name), valid=cal_multi_sha(vcf_files))
    rm_dir(temp_dir)


def verify_data(data_file: str, source_files: str or list, mode: str) -> bool:
    config = load_config()['parameter']
    sha_cal = cal_sha if mode == 'single' else \
        partial(cal_multi_sha, vcf_only=False,
                re=config['need_suf'].replace('.', '\\.') + "|" +
                   config['gwas_suf'].replace('.', '\\.') + '|setting',
                include=True)
    with open(data_file) as f:
        if f.read(4) == 'DONE':
            if f.read(61) == sha_cal(source_files)[3:]:
                return True
    return False


def select_columns(file: str, sep: str, need_col: list, out_file: str, method='number', header=True, skipspace=False):
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
                    first = False


def same_status(files: list):
    row = None
    for file in files:
        with open(file, 'rb') as f:
            if row:
                if row != f.readline():
                    return False
            else:
                row = f.readline()
    return True


def trans_gz(gz_file: str, out_dir: str):
    if gz_file.split('.')[-1] == 'gz':
        with gzip_open(gz_file) as fr:
            with open(os.path.join(out_dir, 'prs_data'), 'w') as fw:
                for line in fr:
                    fw.write(line.decode())
        return os.path.join(out_dir, 'prs_data')
    else:
        return gz_file


def get_ref_index(vcf_file: str):
    id_index = vcf_file + '.idi'
    if os.path.isfile(id_index):
        # with open(id_index, 'r') as f:
        #     if f.read(4) == 'DONE':
        #         if f.read(61) == cal_sha(vcf_file)[3:]:
        #             return id_index
        # if verify_data(id_index, vcf_file, 'single'):
        return id_index
    with open(id_index, 'w') as fw:
        nopen = open_gz(vcf_file)
        fw.write(f'#{cal_sha(vcf_file)}\n')
        with nopen(vcf_file, 'rb') as fr:
            global pattern
            i = 0
            for line in fr:
                i += 1
                if pattern.search(line):
                    fw.write(f'{pattern.search(line).group().decode()}\t{i}\n')
        fw.seek(0)
        fw.write('DONE')
    return id_index


def xlim_adjust(mean: float, x_range: tuple):
    left, right = (mean - x_range[0], x_range[1] - mean)
    return (x_range[0], mean + left) if left > right else (mean - right, x_range[1])


def get_type_list(data_dir: str):
    try:
        return dict(load_txt(os.path.join(data_dir, 'Type.txt')))
    except FileNotFoundError:
        raise Exception("Can not find type list file.")


def get_gt(vcf_gt: str, ref: str, alt: str, snp=None, warning=True):
    # if sum(not i.isnumeric() for i in vcf_gt) != 1:
    #     logging.warning('Not biallelic alleles appear')
    for i in vcf_gt:
        if i == '.':
            logging.warning(f'Find "." in "{snp}", with which it may be difficult to compare the genotypes.')
            return {'*'}
        if i == '*':
            if warning:
                logging.warning(f'Find "*" in "{snp}", with which it may be difficult to compare the genotypes.')
        elif not i.isnumeric():
            gt_n = vcf_gt.split(i)
            return set(alt.split(',')[int(i) - 1] if int(i) else ref for i in gt_n)
    return set(alt.split(',')[int(vcf_gt) - 1] if int(vcf_gt) else ref)


def load_txt(description_text: str):
    with open(description_text, encoding=load_config()['parameter']['encoding']) as f:
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


def get_cal_file(code: str):
    config = load_config()['parameter']
    if os.path.isfile(code + config['need_suf']):
        return True, False, code + config['need_suf']
    elif os.path.isfile(code + config['gwas_suf']):
        return True, True, code + config['gwas_suf']
    elif os.path.isfile(code + config['gwas_suf'] + '.gz'):
        return True, True, code + config['gwas_suf'] + '.gz'
    else:
        return False, None, None


def get_snp_list(data_dir: str, cal_ind=True, prs_ind=False):
    config = load_config()['parameter']
    snp_list = set([])
    snp_column = None
    try:
        os.chdir(data_dir)
        for type_dir in [item for item in os.listdir('.') if os.path.isdir(item)]:
            os.chdir(type_dir)
            for ind_dir in [item for item in os.listdir('.') if os.path.isdir(item)]:
                os.chdir(ind_dir)
                exist, file_type, file = get_cal_file(ind_dir)
                if exist:
                    if cal_ind and not file_type:
                        with open(file, 'r', encoding=config['encoding']) as f:
                            for line in f:
                                snp_list.add(line.strip('\n\r').split(config['text_sep'])[0])
                    elif prs_ind and file_type:
                        header = True
                        if file:
                            nopen = open_gz(file)
                            with nopen(file, 'rb') as f:
                                for line in f:
                                    if header:
                                        if not snp_column:
                                            snp_column = line.split(config['sep'].encode()).index(
                                                config['SNP'].encode())
                                        header = False
                                        continue
                                    else:
                                        snp_list.add(line.split(config['sep'].encode())[snp_column].decode())
                else:
                    logging.warning(f'Code({ind_dir}) has no corresponding algorithm!')
                os.chdir('..')
            os.chdir('..')
    finally:
        os.chdir(raw_dir)
    return snp_list


def get_columns_num(file: str, need_columns: list):
    sep = load_config()['parameter']['sep']
    nopen = open_gz(file)
    with nopen(file, 'rb') as f:
        header = f.readline()
    header = header.strip().split(sep.encode())
    return [header.index(i.encode()) + 1 for i in need_columns]


def get_prs_res(result_file: str):
    with open(result_file, 'r') as f:
        temp = f.readlines()
    os.remove(result_file)
    return exp(float(temp[-1].strip().split(' ')[-1]))


def cal_trait(trait: str or None, gt: set, ref: str, inter: str, no_inter: str or None):
    lan_lib = load_config()['lan_lib']
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
                if trait == 'And_yes':
                    return outcome
                else:
                    return 'And_no'
            else:
                outcome = inter if equal_gt(gt, trans_gt(ref)) else inter \
                    if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
                return 'Normal' if outcome in lan_lib[None] else outcome
        elif trait == 'And_yes':
            outcome = inter if equal_gt(gt, trans_gt(ref)) else inter \
                if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
            return 'Normal' if outcome in lan_lib[None] else outcome
        elif trait == 'And_no':
            outcome = no_inter
            return 'Normal' if outcome in lan_lib[None] else outcome
        else:
            return trait


# QC
# MAF
def qc_maf(failed_snps: set, file_path: str, report_dir: str):
    maf_threshold = load_config()['parameter']['qc_maf']
    data = pd.read_csv(os.path.join(file_path, 'qc.frq'), ' ', skipinitialspace=True)
    failed_snps.update(data.loc[data.MAF < maf_threshold, 'SNP'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.MAF, bins=50)
    plt.text(assign_pos(maf_threshold, plt.axis(), right_offset=0.02, left_offset=0.24), plt.axis()[3] * 0.184,
             f'MAF threshold: {maf_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(maf_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Minor allele frequency")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_MAF.png'), bbox_inches='tight')
    plt.close()


# SNPs_MISS
def qc_vmiss(failed_snps: set, file_path: str, report_dir):
    snp_miss_threshold = load_config()['parameter']['qc_vmiss']
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
    plt.savefig(os.path.join(report_dir, f'QC_var_miss.png'), bbox_inches='tight')
    plt.close()


# FMISS
def qc_smiss(failed_samples: set, file_path: str, report_dir):
    individual_miss_threshold = load_config()['parameter']['qc_smiss']
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
    plt.savefig(os.path.join(report_dir, f'QC_ind_miss.png'), bbox_inches='tight')
    plt.close()


# HARDY
def qc_hardy(failed_snps: set, file_path: str, report_dir):
    hardy_threshold = load_config()['parameter']['qc_hardy']
    data = pd.read_csv(os.path.join(file_path, 'qc.hwe'), ' ', skipinitialspace=True)
    failed_snps.update(data.loc[data.P < 10 ** -hardy_threshold, 'SNP'])
    data.P[data.P == 0] = 1e-295
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(-log10(data.P), bins=50)
    plt.text(assign_pos(hardy_threshold, plt.axis(), right_offset=0.02, left_offset=0.285), plt.axis()[3] * 0.816,
             f'HWE exact test p-value\nthreshold: 1e-{hardy_threshold}', bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(hardy_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("-log$_{10}$(HWE exact test p-value)")
    plt.ylabel("Number of variants")
    plt.savefig(os.path.join(report_dir, f'QC_hardy.png'), bbox_inches='tight')
    plt.close()


# HET
def qc_het(failed_samples: set, file_path: str, report_dir):
    het_threshold = load_config()['parameter']['qc_het']
    data = pd.read_csv(os.path.join(file_path, 'qc.het'), ' ', skipinitialspace=True)
    het_mean, het_sd = data.F.mean(), data.F.std()
    failed_samples.update(data.loc[(data.F - het_mean > het_sd * het_threshold) |
                                   (data.F - het_mean < -het_sd * het_threshold), 'IID'])
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.hist(data.F, bins=50)
    plt.text(assign_pos(het_mean - het_sd * het_threshold, plt.axis(), right_offset=0.02, left_offset=0.285),
             plt.axis()[3] * 0.816,
             f'Heterozygosity rate\nthreshold: mean ±{het_threshold} sd',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.axvline(het_mean + het_sd * het_threshold, c='r', lw=1.5, ls='--')
    plt.axvline(het_mean - het_sd * het_threshold, c='r', lw=1.5, ls='--')
    plt.xlabel("Heterozygosity rate")
    plt.ylabel("Number of individuals")
    plt.xlim(*xlim_adjust(het_mean, plt.xlim()))
    plt.savefig(os.path.join(report_dir, f'QC_het.png'), bbox_inches='tight')
    plt.close()


# SEX
def qc_sex(failed_samples: set, file_path: str, report_dir):
    male_threshold = load_config()['parameter']['qc_male_F']
    female_threshold = load_config()['parameter']['qc_female_F']
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
    plt.savefig(os.path.join(report_dir, f'QC_sex.png'), bbox_inches='tight')
    plt.close()


# RELATEDNESS
def qc_relatedness(failed_samples: set, file_path: str, report_dir,):
    pihat_threshold = load_config()['parameter']['qc_pihat']
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
    plt.savefig(os.path.join(report_dir, f'QC_relatedness.png'), bbox_inches='tight')
    plt.close()


@progress_value(int(10 * average_progress()))
def vep_query(search_ids: list):
    import json
    import requests
    vep_api = "https://rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    result = {}
    for i in tqdm(range(0, len(search_ids), 200)):
        data = json.dumps({'ids': search_ids[i: i + 200]})
        r = requests.post(vep_api, headers=headers, data=data)
        if r.ok:
            decoded = r.json()
            result.update({res['id']: {'Variant type': res['most_severe_consequence']}
                           for res in decoded})
        else:
            r.raise_for_status()


def load_clinvar(clinvar_vcf: str, clinvar_data: dict, human: Human):
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


def recognize_color_clinvar(row: pd.Series):
    risk = True if 'athogenicity' in row['Allele type'] else True if 'Risk' in row['Allele type'] else True \
        if 'Drug response' in row['Allele type'] else None
    return recognize_color(row, risk, 'Allele', 'Genotype', na_value=None)


def recognize_color_pharmgkb(row: pd.Series):
    risk = True if 'increased' in row['Sentence'] else False if 'decreased' in row['Sentence'] else None \
        if 'not' in row['Sentence'] else True
    return recognize_color(row, risk, 'Alleles', 'Genotype')


def recognize_color(row: pd.Series, risk: bool, allele_col: str, genotype_col: str, na_value=None):
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


def load_pharm(pharm_data: str, human: Human):
    pharm_data = pd.read_csv(pharm_data, sep='\t', error_bad_lines=False)
    res = pharm_data.loc[pharm_data.Variant.isin(human.gt_data) & (pharm_data.Significance == 'yes'),
                         ['Variant', 'Alleles', 'Sentence', 'PMID']]
    res['Genotype'] = res.apply(lambda a: trans_gt(human.gt_data[a['Variant']], connector=''), axis=1)
    res['Color'] = res.apply(recognize_color_pharmgkb, axis=1)
    return res.dropna().sort_values(['Color', 'Variant'], ascending=[False, True]).reset_index(drop=True)


@progress_value(int(10 * average_progress()))
def trans_clinvar_res(clinvar_res: dict) -> pd.DataFrame:
    res = pd.DataFrame(clinvar_res).T
    res.reset_index(inplace=True)
    res.columns = ['Variant'] + list(res.columns)[1:]
    res = res.loc[res.Allele != '.',]
    res['Color'] = res.apply(recognize_color_clinvar, axis=1)
    return res.dropna().sort_values(['Color', 'Variant'], ascending=[False, True]).reset_index(drop=True)


def get_clinvar_data(data_dir: str):
    # download data
    if not os.path.isfile(os.path.join(data_dir, 'clinvar.vcf.gz')):
        urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz',
                    os.path.join(data_dir, 'clinvar.vcf.gz'))

    if not os.path.isfile(os.path.join(data_dir, 'clinvar_papu.vcf.gz')):
        urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_papu.vcf.gz',
                    os.path.join(data_dir, 'clinvar_papu.vcf.gz'))


@progress_value(int(10 * average_progress()))
def get_pharmgkb_data(data_dir: str):
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


def write_failed_qc_list(failed_snps: set, failed_samples: set):
    with open(os.path.join('ref_res', "failed_qc_samples"), 'w') as f:
        for sample in failed_samples:
            f.write(sample + ' ' + sample + '\n')
    with open(os.path.join('ref_res', "failed_qc_snps"), 'w') as f:
        for snp in failed_snps:
            f.write(snp + '\n')


@use_time('Get reference result data')
def get_ref_cal(data_dir: str, vcf_files: list, type_list: dict, failed_snps: set, failed_samples: set,
                data_path='ref_res'):
    ref_code_file = os.path.join(data_path, 'population_code_res.csv')
    ref_rs_file = os.path.join(data_path, 'population_rs_res.csv')
    ref_config = os.path.join(data_path, 'setting')
    config = load_config()['parameter']
    config = {key: config[key] for key in config if key not in ['file_type', 'logo_dir', 'show_code', 'qc_maf',
                                                                'qc_vmiss', 'qc_hardy', 'qc_het', 'qc_male_F',
                                                                'qc_female_F', 'qc_pihat']}
    use_ref_qc = config['use_qc_ref']
    with open(ref_config, 'wb') as f:
        dump(config, f)
    if os.path.isfile(ref_code_file) and os.path.isfile(ref_rs_file):
        if verify_data(ref_code_file, vcf_files + [data_dir] + [ref_config], 'multi'):
            if same_status([ref_code_file, ref_rs_file]):
                return pd.read_csv(ref_rs_file, skiprows=1, index_col=0), \
                       pd.read_csv(ref_code_file, skiprows=1, index_col=0)
    outcomes_rs = pd.DataFrame()
    outcomes_code = pd.DataFrame()
    gt_data = get_ref_gt(vcf_files, data_dir)
    if use_ref_qc:
        gt_data = gt_data.loc[~gt_data.index.isin(failed_snps), ~gt_data.columns.isin(failed_samples)]
        write_failed_qc_list(failed_snps, failed_samples)
    with open(ref_code_file, 'w', newline='', encoding='UTF-8') as fc:
        with open(ref_rs_file, 'w', newline='', encoding='UTF-8') as fr:
            sha = cal_multi_sha(vcf_files + [data_dir] + [ref_config], vcf_only=False,
                                re=config['need_suf'].replace('.', '\\.') + "|" +
                                   config['gwas_suf'].replace('.', '\\.') + '|setting',
                                include=True)
            fc.write(f'#{sha}\n')
            fr.write(f'#{sha}\n')
            try:
                os.chdir(data_dir)
                for type_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
                    os.chdir(type_dir)
                    for ind_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
                        os.chdir(ind_dir)
                        exist, file_type, file = get_cal_file(ind_dir)
                        if exist:
                            if file_type:
                                res = get_ref_prs(file, vcf_files, use_ref_qc)
                                if outcomes_code.empty:
                                    outcomes_code[ind_dir] = res
                                else:
                                    res = res.rename(ind_dir)
                                    outcomes_code = outcomes_code.merge(res.apply(lambda x: exp(x)), left_index=True,
                                                                        right_index=True)
                            else:
                                with open(file, 'r', encoding=config['encoding']) as f:
                                    for line in f:
                                        if line.strip("\n\r"):
                                            line = line.strip("\n\r").split(config['text_sep'])
                                            rs = line[0]
                                            if rs in gt_data.index:
                                                if type_list[type_dir] == 'quan':
                                                    outcome = gt_data.loc[rs].apply(
                                                        lambda gt: cal_beta(gt, *line[1:3]))
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
                                                    temp = list(map(lambda trait, gt: cal_trait(trait, gt, *line[1:4]),
                                                                    outcomes_code[ind_dir], gt_data.loc[rs]))
                                                    outcomes_code[ind_dir] = temp
                        os.chdir('..')
                    os.chdir('..')
            finally:
                os.chdir(raw_dir)
                outcomes_rs.to_csv(fr)
                fr.seek(0)
                fr.write('DONE')
            outcomes_code.to_csv(fc)
            fc.seek(0)
            fc.write('DONE')
    return outcomes_rs, outcomes_code


@use_time('Get population prs result')
def get_ref_prs(prs_data: str, vcf_files: list, qc_ref: bool):
    config = load_config()['parameter']
    ref_structure = config['ref_structure']
    temp_dir = mkdtemp(suffix='prs')
    code = os.path.basename(prs_data).split('.')[0]
    # Clump
    clump_p2 = config["clump-p2"] if config["clump-p1"] < config["clump-p2"] \
        else config["clump-p1"]
    run_plink_cmd(f'--vcf {ref_structure} --clump-p1 {str(config["clump-p1"])} '
                  f'--clump-p2 {str(clump_p2)} '
                  f'--clump-r2 {str(config["clump-r2"])} --clump-kb {str(config["clump-kb"])} '
                  f'--clump {prs_data} --clump-snp-field {config["SNP"]} '
                  f'--clump-field {config["P"]} --out {os.path.join(temp_dir, code)}', plink='plink')
    # Produce plink need data
    select_columns(os.path.join(temp_dir, code + '.clumped'), ' ', [2],
                   os.path.join(temp_dir, 'valid.snp'), header=False, skipspace=True)
    select_columns(prs_data, config['sep'], [0, 4], os.path.join(raw_dir, 'ref_res', f'{code}.SNP.pvalue'))
    with open(os.path.join(raw_dir, 'ref_res', 'need_range_list'), 'w') as f:
        f.write(f'{config["p_threshold"]} 0 {config["p_threshold"]}\n')
    # PRS calculate
    columns_num = get_columns_num(prs_data, [config[i] for i in ['SNP', 'EA', 'BETA']])
    prs_data = trans_gz(prs_data, temp_dir)
    qc_str = f"--remove {os.path.join(raw_dir, 'ref_res', 'failed_qc_samples')} " \
             f"--exclude {os.path.join(raw_dir, 'ref_res', 'failed_qc_snps')} " if qc_ref else ""
    for fid, file in enumerate(vcf_files):
        run_plink_cmd(f'--vcf {file} {qc_str}--score {prs_data} '
                      f'{" ".join([str(i) for i in columns_num])} header '
                      f'--q-score-range {os.path.join(raw_dir, "ref_res", "need_range_list")} '
                      f'{os.path.join(raw_dir, "ref_res", code + ".SNP.pvalue")} '
                      f'--extract {os.path.join(temp_dir, "valid.snp")} --out {os.path.join(temp_dir, "result_prs")}',
                      plink='plink')
        if fid == 0:
            res_data = pd.read_csv(os.path.join(temp_dir, f"result_prs.{config['p_threshold']}.profile"), ' ',
                                   skipinitialspace=True, index_col=1)['SCORE']
        else:
            res_data = res_data + pd.read_csv(os.path.join(temp_dir,
                                                           f"result_prs.{config['p_threshold']}.profile"),
                                              ' ', skipinitialspace=True, index_col=1)['SCORE']
    rm_dir(temp_dir)
    return res_data


def risk_cal(string: str, risk_score: float):
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


def read_line_num(gwas_file: str, header=False) -> int:
    nopen = open_gz(gwas_file)
    line_num = -1 if header else 0
    with nopen(gwas_file, 'rb') as f:
        for line in f:
            line_num += 1
    return line_num
