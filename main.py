import os
import getopt
from sys import exit as sys_exit
import pandas as pd
import numpy as np
from xlrd import open_workbook
import logging
from collections import OrderedDict
from gzip import open as gzip_open
from re import compile, sub
from time import time, strftime
from tempfile import mkdtemp
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from jinja2 import Environment, FileSystemLoader
from shutil import copy
import subprocess
from math import exp
from pickle import dump, load
from multiprocessing import Pool
from functools import partial


description = "Usage: python perhaps.py -i INPUT_FILE [--ftype FILE_TYPE] -o OUTPUT_DIR [--trait TRAIT_FILE]\n" \
              "[--qual QUALITATIVE_FILE] [--quan QUANTITATIVE] [--ref-dir REFERENCE_DIR] [--file-dir REFERENCE_DIR]\n" \
              "[--sep GWAS_SEPARATOR] [-p GWAS_PTHRESHOLD] [--ref-dir REFERENCE_DIR] [--GWAS_config KEY=VALUE ...]"

columns = {'columns_ind': OrderedDict({'code': 0, 'name': 1, 'sex': 3, 'Disease_description': 5}),
           'columns_quan': OrderedDict({'code': 0, 'snp': 1, 'reference': 2, 'beta': 3}),
           'columns_qual': OrderedDict({'code': 0, 'snp': 1, 'reference': 2, 'interpretation': 3, 'reverse-inter': 4})}
lan_lib = {'sex': {'male': ['male', 'Male', '男', '男性', '1', 1], 'female': ['female', 'Female', '女', '女性', '2', 2]},
           None: ['', '无', 'None', 'No', None], 'if': ['条件', 'if']}
GWAS_setting = {'SNP': 'SNP', 'EA': 'EA', 'P': 'P', 'BETA': 'BETA', 'sep': '\t', 'p_threshold': 1e-5,
                'clump-p1': 1, 'clump-r2': 0.1, 'clump-kb': 250, 'clump-p2': 0.01}
config = {'file_type': 'vcf', 'ref_structure': f'{os.path.join("reference", "ld_ref", "hapmap3.vcf.gz")}',
          'logo_dir': f'{os.path.join("database", "logo")}'}

pattern = compile(rb'(?<=[\t])rs[0-9]*(?=[\t;])')
pat_header1 = compile(rb'^##')
pat_header2 = compile(rb'^#')
ref_code = pd.DataFrame()
ref_rs = pd.DataFrame()
logo_code = [file.split('.')[0] for file in os.listdir(config['logo_dir'])]

if not os.path.isdir('log'):
    os.mkdir('log')
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
logger.addHandler(ch)


def log_start():
    log_name = f'log/{strftime("%Y%m%d%H%M%S")}.log'
    fh = logging.FileHandler(log_name)
    formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info('Logger start.')
    return log_name


def use_time(process_name=None):
    """
    Record the running time of the function
    :param process_name: the name of this process
    :return: A informational log which records the time
    """

    def decorator(func):
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time()
            logger.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                        f' start')
            fun_res = func(*args, **kwargs)
            logger.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                        f' used time: {time() - start:.2f}s')
            return fun_res

        return wrapper

    return decorator


@use_time()
def run_plink_cmd(cmd: str, plink='plink2', delete_log=True) -> None:
    """
    Creat a subprocess which not display in GUI.
    :param delete_log:
    :param plink:
    :param cmd: cmd string
    :return: A finished subprocess with no result
    """
    if '--out' in cmd:
        out_file = cmd.split('--out ')[-1].split(' ')[0]
    else:
        out_file = os.path.join('.', plink)
    cmd = f'{os.path.join(".", "bin", plink)} ' + cmd
    a = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, sterr = a.communicate()
    a.wait()
    if delete_log:
        os.remove(out_file + '.log')
    if b'Error' in sterr:
        raise Exception("PlinkError:\n" + sterr.decode())


def initial_res_dir(output: str):
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, 'html_files')
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, 'img')
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, '..', 'dist_plot')
    if not os.path.isdir(output):
        os.mkdir(output)


def initial_ref_data(qual_list: list, quan_list: list, vcf_list: list) -> None:
    """
    Get reference population data from existing data or generating newly
    :param quan_list: Quantitative files
    :param vcf_list: Vcf files in reference population directory
    :return: Two dataframe: one recorded the population snps result, the other recorded the code result
    """
    global ref_rs, ref_code
    if not os.path.isdir('ref_res'):
        os.mkdir('ref_res')
    assert vcf_list, 'There is no vcf file in population reference directory.'
    ref_rs, ref_code = get_ref_cal(qual_list, quan_list, vcf_list)


class Human(object):
    def __init__(self, name: str, file: str, temp_dir: str):
        self.name = name
        self.gt_data = {}
        self.ind_data = {}
        self.vcf, self.sex = recode_and_sex_impute(file, temp_dir)

    def set_gt(self, snp_id: str, gt: set):
        self.gt_data[snp_id] = gt

    def set_ind(self, code: str, name: str, sex: str, itype: str, **other):
        if trans_sex(sex) is None or trans_sex(sex) == self.sex:
            self.ind_data[code] = Ind(code, name, sex, itype, **other)

    def cal_ind(self, ftype: str, code: str, snp: str, *cal_col):
        if code in self.ind_data:
            self.ind_data[code].set_ftype(ftype)
            if snp in self.gt_data:
                self.ind_data[code].cal(snp, self.gt_data[snp], *cal_col)
            else:
                self.ind_data[code].mis_snp(snp, cal_col[1])
        else:
            logger.warning(f'Code ({code}) is missing in index excel!')

    def get_prs_data(self, code: str, res: float, detail: str):
        if code in self.ind_data:
            self.ind_data[code].add_prs_res(res, detail)
        else:
            logger.warning(f'Code ({code}) is missing in index excel!')

    @use_time('Add population distribution')
    def add_distribution(self, output: str):
        global ref_code
        assert not ref_code.empty, "Can not get reference data"
        output = os.path.join(output, 'html_files', 'dist_plot')
        for ind in self.ind_data.values():
            ind.add_quan_dist(output)

    @use_time('Export report result')
    def export_res(self, output):
        result = {}
        for ind in self.ind_data.values():
            if not ind.status:
                if not ind.detail:
                    logger.warning(f'Code ({ind.code}) has no corresponding algorithm!')
                    ind.detail.append(f'Code ({ind.code}) has no corresponding algorithm!')
                ind.outcome = 'Undected'
            if ind.itype not in result:
                result[ind.itype] = []
            result[ind.itype].append(ind.report(output))
        return result


class Ind(object):
    def __init__(self, code: str, name: str, sex: str, itype: str, **other):
        self.code = code
        self.name = name
        self.sex = sex
        self.itype = itype
        self.ftype = None
        self.outcome = None
        self.other = other
        self.detail = []
        self.status = False
        self.decide = False
        self.distribution = None
        self.picture = True if self.code in logo_code else False

    def set_ftype(self, ftype: str):
        if not self.ftype:
            if ftype == 'quan':
                self.outcome = 1
                self.distribution = {'Average': None, 'Percentage': None, 'Plot': None, 'Status': False}
            else:
                self.distribution = {'Plot': None, 'Status': False}
            self.ftype = ftype

    def add_prs_res(self, res: float, detail: str):
        self.ftype = 'quan'
        self.outcome = res
        self.distribution = {'Average': None, 'Percentage': None, 'Plot': None, 'Status': False}
        self.detail = [detail]
        self.status = True

    def cal(self, snp: str, gt: set, ref: str, inter: str or float or int, no_inter=None):
        if self.ftype == 'qual':
            if not self.decide:
                if inter not in lan_lib['if']:
                    if self.decide is not None:
                        self.outcome = inter if gt == trans_gt(ref) else \
                            inter if gt == trans_gt(gene_filp(ref)) else no_inter
                        self.outcome = 'normal' if self.outcome in lan_lib[None] else self.outcome
                        if {snp: trans_gt(gt)} not in self.detail:
                            self.detail.append({snp: trans_gt(gt)})
                        if self.outcome != 'normal':
                            self.decide = True
                        self.status = True
                else:
                    self.decide = False if gt == trans_gt(ref) else None
                    if self.decide is False:
                        if {snp: trans_gt(gt)} not in self.detail:
                            self.detail.append({snp: trans_gt(gt)})
        elif self.ftype == 'quan':
            outcome = cal_beta(gt, ref, inter)
            self.outcome *= outcome
            for item in self.detail:
                if type(item) == dict:
                    if snp in item.keys():
                        logger.warning(f"There are duplicate snps ({snp}) when calculate the quantitative trait!")
                        self.outcome /= item['outcome']
                        item['beta'] = inter
                        item['effect allele'] = ref
                        item['outcome'] = outcome
            if {snp: trans_gt(gt), 'effect allele': ref, 'beta': inter, 'outcome': outcome} not in self.detail:
                self.detail.append({snp: trans_gt(gt), 'effect allele': ref, 'beta': inter, 'outcome': outcome})
            self.status = True

    def mis_snp(self, snp: str, inter: str or float or int):
        warning = f'Snp ({snp}) can not find in vcf file!'
        if self.ftype == 'quan':
            global ref_rs
            if f'{self.code}_{snp}' in ref_rs.columns:
                beta = ref_rs[f"{self.code}_{snp}"].mean()
                warning += f' Using reference population average: {beta}'
                self.outcome *= beta
        elif self.ftype == 'qual':
            if inter in lan_lib['if']:
                self.decide = None
                if warning in self.detail:
                    self.detail.remove(warning)
                warning += f' Code {self.code} algorithm can not process further.'
        logger.warning(warning)
        if warning not in self.detail:
            self.detail.append(warning)
        self.status = False

    def add_quan_dist(self, output: str):
        if self.code in ref_code.columns:
            if self.ftype == 'quan':
                self.distribution['Average'] = ref_code[self.code].mean()
                self.distribution['Percentage'] = sum(ref_code[self.code] < self.outcome) / ref_code.shape[0] * 100
                self.distribution['Plot'] = quan_dist_plot(self.code, self.outcome, ref_code[self.code],
                                                           self.distribution['Percentage'], output)
                self.distribution['Status'] = True
            else:
                self.distribution['Plot'] = qual_dist_plot(self.code, self.outcome, self.name, output)
                self.distribution['Status'] = True

    def report(self, output):
        if self.picture:
            self.picture = os.path.join('html_files', 'img', os.listdir(config['logo_dir'])[logo_code.index(self.code)])
            'no_pic.jpg'
            copy(os.path.join(config['logo_dir'], os.path.basename(self.picture)), os.path.join(output, self.picture))
        else:
            self.picture = os.path.join('html_files', 'img', 'no_pic.jpg')
        if self.distribution['Status']:
            del self.distribution['Status']
            return {'type': self.ftype, 'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail,
                    'Distribution': self.distribution, **self.other, 'Status': self.status, 'Picture': self.picture}
        return {'type': self.ftype, 'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail, **self.other,
                'Status': self.status,
                'Picture': self.picture}


def assign_pos(value: float, plt_size: tuple):
    ran = plt_size[1] - plt_size[0]
    p = (value - plt_size[0]) / ran
    return value + 0.04 * ran if p < 0.5 else value - 0.29 * ran


def risk_level(per: float):
    assert 0 <= per <= 100, 'Error percentage!'
    cuts = [0, 1, 5, 10, 20, 40, 60, 80, 90, 95, 99, 100]
    for level, cut in enumerate(cuts):
        if per < cut:
            return level
    return 11


def is_excel(file: str) -> bool:
    return True if file.split('.')[-1] == 'xlsx' or file.split('.')[-1] == 'xls' else False


def qual_dist_plot(code: str, res: str, name: str, report_dir: str):
    population = ref_code[code].value_counts()
    labels = population.index
    labels = [label.split(":")[0] for label in labels]

    plt.figure(figsize=(8, 6), dpi=400)
    ax1, ax2 = plt.subplot(1, 3, (1, 2)), plt.subplot(1, 3, 3)
    patches, texts = ax1.pie(population, radius=1.2, pctdistance=1.13)
    ax1.axis('equal')
    ax2.axis('off')
    ax2.legend(patches, labels, loc='center')
    plt.suptitle(f'Distribution of {name} in reference population', size=12, weight="bold", y=0.95)

    props = population.apply(lambda a: f'{a / sum(population) * 100:.2f}%')
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="->"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(patches):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        if 'ang_r' not in locals().keys():
            y_dif = 0
        elif (ang - 180) * (ang_r - 180) >= 0:
            if ang - ang_r < 10:
                y_dif += 0.3 / (ang - ang_r)
            else:
                y_dif -= 1 / (ang - ang_r)
                y_dif = y_dif if y_dif >= 0 else 0
        else:
            y_dif = 0
        ang_r = ang
        x = np.cos(np.deg2rad(ang))
        y = np.sin(np.deg2rad(ang))
        y_text = y + y_dif
        y_text = (y_text - 0.8) * 0.33 + 0.8 if y_text > 0.8 else y_text
        horizontalalignment = "right" if ang <= 180 else "left"
        connectionstyle = f"angle, angleA={0 + np.rad2deg(np.arctan(1.5 * y_dif))}, " \
                          f"angleB={ang + np.rad2deg(np.arcsin(1.5 * y_dif))}" if abs(y) > 0.5 else 'arc3'
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax1.annotate(props[i], xy=(1.18 * x, 1.18 * y), xytext=(1.5 * np.sign(x), round(1.5 * y_text, 2)),
                     horizontalalignment=horizontalalignment, **kw)

    plt.tight_layout()
    plt.savefig(os.path.join(report_dir, f'{code}.png'))
    plt.close()
    return os.path.join('html_files', 'dist_plot', f'{code}.png')


def quan_dist_plot(code: str, value: float, ref_data: pd.Series, per: float, report_dir: str):
    from numpy import log10
    plt.hist(log10(ref_data), bins=100, weights=[1. / len(ref_data)] * len(ref_data))
    plt.axvline(log10(value), c='r', lw=2, ls='--')
    plt.text(assign_pos(log10(value), plt.axis()), plt.axis()[3] * 0.816,
             f'Percentage: {per:.0f}% \nLevel: {risk_level(per)}',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda a, b: f'{a * 100:.0f}%'))
    plt.xlabel("Logarithm of Risk score")
    plt.ylabel("Percentage")
    plt.savefig(os.path.join(report_dir, f'{code}.png'))
    plt.close()
    return os.path.join('html_files', 'dist_plot', f'{code}.png')


def cal_beta(gt: set, ref: str, beta: int or float):
    return 1 if ref not in gt else beta ** 2 if len(gt) == 1 else beta


def cal_trait(trait: str or None, gt: set, ref: str, inter: str, no_inter: str or None):
    if trait is None:
        if inter in lan_lib['if']:
            return 'Normal' if gt == trans_gt(ref) else 'Normal' if gt == trans_gt(gene_filp(ref)) else None
        else:
            return None
    else:
        if trait == 'Normal':
            if inter in lan_lib['if']:
                return 'Normal' if gt == trans_gt(ref) else 'Normal' if gt == trans_gt(gene_filp(ref)) else None
            else:
                outcome = inter if gt == trans_gt(ref) else inter if gt == trans_gt(gene_filp(ref)) else no_inter
                return 'Normal' if outcome in lan_lib[None] else outcome
        else:
            return trait


def trans_sex(sex: str or bool):
    if type(sex) == str:
        return True if sex in lan_lib['sex']['male'] else False if sex in lan_lib['sex']['female'] else None
    elif type(sex) == bool:
        return 'male' if sex else 'female'


def trans_gt(gt: set or str):
    if type(gt) == set:
        return '/'.join(list(gt)[0] * 2 if len(gt) == 1 else list(gt))  # Todo: not all biallelic
    elif type(gt) == str:
        for i in gt:
            if not i.isalpha():
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def gene_filp(raw: str):
    base_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    return ''.join([base_dict[char] if char in base_dict.keys() else char for char in raw])


def select_list(ob_list: list, index):
    # assert len(ob_list) < max(index), 'list index out of range'
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def sub_col(ob_list: list, index: OrderedDict, code=True):
    if not code:
        index.pop('code')
    return {key: ob_list[index[key]] if index[key] < len(ob_list) else None for key in index}


@use_time()
def recode_and_sex_impute(file: str, temp_dir: str, file_type=config['file_type']):
    if file_type == 'vcf' or file_type == 'vcf.gz':
        read_txt = f"--vcf {file} "
    elif file_type == '23andme':
        read_txt = f"--23file {file} "
    else:
        raise TypeError("Unsupported input file")
    try:
        run_plink_cmd(read_txt +
                      "--impute-sex y-only "  # Should be reconsider
                      "--recode vcf "
                      f"--out {os.path.join(temp_dir, 'temp')}",
                      plink='plink')
        data = pd.read_csv(os.path.join(temp_dir, 'temp.sexcheck'), ' ', skipinitialspace=True)
    except Exception as e:
        logger.warning('Sex impute failed:\n' + f'{e.args}')
        run_plink_cmd(f"--vcf {file} "
                      "--recode vcf "
                      f"--out {os.path.join(temp_dir, 'temp')}",
                      plink='plink')
        return os.path.join(temp_dir, 'temp.vcf'), None
    else:
        # os.system('rm sex_impute.*')
        return os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.SNPSEX[0] == 1 else 'Female'


def cal_sha(file: str):
    from hashlib import sha256
    with open(file, 'rb') as f:
        return sha256(f.read()).hexdigest()


@use_time("Calculate files' sha value")
def cal_multi_sha(file_list: str or list, vcf_only=True):
    from hashlib import sha256
    value = sha256(b'')
    if type(file_list) == str:
        file_list = os.listdir('dir')
    for file in file_list:
        if file.split('.')[-1] == 'gz' or file.split('.')[-1] == 'vcf' or not vcf_only:
            with open(file, 'rb') as f:
                value.update(f.read())
    return value.hexdigest()


def open_gz(file: str):
    if file.split('.')[-1] == 'gz':
        return gzip_open
    else:
        return open


@use_time('Get reference index')
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


@use_time('Extract_snp from reference data')
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
    if nopen:
        logger.info(f"Step 1: Get {vcf_file} index")
        idi_file = get_ref_index(vcf_file)
        logger.info(f"Step 2: Get {vcf_file} need snps' rows' number")
        with open(idi_file) as f:
            for line in f:
                temp = line.strip().split('\t')
                if temp[0] in need_snp_list:
                    need_lines.append(int(temp[1]))
                    need_snp_list.remove(temp[0])
        logger.info(f"Step 3: Extract {vcf_file} need snps")
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
        logger.info(f"{vcf_file} finish extract.")


@use_time('Merge vcf files')
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


@use_time('Get reference GT data')
def get_ref_gt(vcf_list: list, excel_files: list,
               filter_vcf=os.path.join('ref_res', 'filter.vcf')) -> pd.DataFrame:
    if os.path.isfile(filter_vcf):
        if not verify_data(filter_vcf, vcf_list, 'multi'):
            get_ref_vcf(vcf_list, get_snp_list(*excel_files))
    else:
        get_ref_vcf(vcf_list, get_snp_list(*excel_files))
    gt_df = get_gt_iid(vcf_list)
    with open(filter_vcf, 'r') as f:
        f.readline()
        for line in f:
            temp = line.strip().split('\t')
            try:
                gt_df.loc[temp[2]] = [get_gt(gt, temp[3], temp[4]) for gt in temp[9:]]
            except IndexError:
                print(temp)
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


@use_time()
def verify_data(data_file: str, source_files: str or list, mode: str) -> bool:
    sha_cal = cal_sha if mode == 'single' else partial(cal_multi_sha, vcf_only=False)
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


@use_time('Get population prs result')
def get_ref_prs(prs_data: str, vcf_files: list, ref_structure=config['ref_structure']):
    temp_dir = mkdtemp(suffix='prs')
    code = os.path.basename(prs_data).split('.')[0]
    # Clump
    clump_p2 = GWAS_setting["clump-p2"] if GWAS_setting["clump-p1"] < GWAS_setting["clump-p2"] \
        else GWAS_setting["clump-p1"]
    run_plink_cmd(f'--vcf {ref_structure} --clump-p1 {str(GWAS_setting["clump-p1"])} '
                  f'--clump-p2 {str(clump_p2)} '
                  f'--clump-r2 {str(GWAS_setting["clump-r2"])} --clump-kb {str(GWAS_setting["clump-kb"])} '
                  f'--clump {prs_data} --clump-snp-field {GWAS_setting["SNP"]} '
                  f'--clump-field {GWAS_setting["P"]} --out {os.path.join(temp_dir, code)}', plink='plink')
    # Produce plink need data
    select_columns(os.path.join(temp_dir, code + '.clumped'), ' ', [2],
                   os.path.join(temp_dir, 'valid.snp'), header=False, skipspace=True)
    select_columns(prs_data, GWAS_setting['sep'], [0, 4], os.path.join('ref_res', f'{code}.SNP.pvalue'))
    with open(os.path.join('ref_res', 'need_range_list'), 'w') as f:
        f.write(f'{GWAS_setting["p_threshold"]} 0 {GWAS_setting["p_threshold"]}\n')
    # PRS calculate
    columns_num = get_columns_num(prs_data, [GWAS_setting[i] for i in ['SNP', 'EA', 'BETA']])
    prs_data = trans_gz(prs_data, temp_dir)
    for fid, file in enumerate(vcf_files):
        run_plink_cmd(f'--vcf {file} --score {prs_data} '
                      f'{" ".join([str(i) for i in columns_num])} header '
                      f'--q-score-range {os.path.join("ref_res", "need_range_list")} '
                      f'{os.path.join("ref_res", code + ".SNP.pvalue")} '
                      f'--extract {os.path.join(temp_dir, "valid.snp")} --out {os.path.join(temp_dir, "result_prs")}',
                      plink='plink')
        if fid == 0:
            res_data = pd.read_csv(os.path.join(temp_dir, f"result_prs.{GWAS_setting['p_threshold']}.profile"), ' ',
                                   skipinitialspace=True, index_col=1)['SCORE']
        else:
            res_data = res_data + pd.read_csv(os.path.join(temp_dir,
                                                           f"result_prs.{GWAS_setting['p_threshold']}.profile"),
                                              ' ', skipinitialspace=True, index_col=1)['SCORE']
    rm_dir(temp_dir)
    return res_data


@use_time('Get reference quantitative result data')
def get_ref_cal(qual_file: list, quan_file: list, vcf_files: list, data_path='ref_res'):
    ref_code_file = os.path.join(data_path, 'population_code_res.csv')
    ref_rs_file = os.path.join(data_path, 'population_rs_res.csv')
    ref_config = os.path.join(data_path, 'setting')
    with open(ref_config, 'wb') as f:
        dump(GWAS_setting, f)
    if os.path.isfile(ref_code_file) and os.path.isfile(ref_rs_file):
        if verify_data(ref_code_file, vcf_files + qual_file + quan_file + [ref_config], 'multi'):
            if same_status([ref_code_file, ref_rs_file]):
                return pd.read_csv(ref_rs_file, skiprows=1, index_col=0), \
                       pd.read_csv(ref_code_file, skiprows=1, index_col=0)
    outcomes_rs = pd.DataFrame()
    outcomes_code = pd.DataFrame()
    gt_data = get_ref_gt(vcf_files, qual_file + quan_file)
    with open(ref_code_file, 'w', newline='', encoding='UTF-8') as fc:
        with open(ref_rs_file, 'w', newline='', encoding='UTF-8') as fr:
            sha = cal_multi_sha(vcf_files + qual_file + quan_file + [ref_config], vcf_only=False)
            fc.write(f'#{sha}\n')
            fr.write(f'#{sha}\n')
            for file in qual_file + quan_file:
                if is_excel(file):
                    workbook = open_workbook(file)
                    for sheet in workbook.sheets():
                        for row in range(1, sheet.nrows):
                            code = sheet.row_values(row)[0] if sheet.row_types(row)[0] == 1 \
                                else f'{sheet.row_values(row)[0]:.0f}'
                            rs = sheet.row_values(row)[1]
                            if rs in gt_data.index:
                                if file in quan_file:
                                    outcome = gt_data.loc[rs].apply(
                                        lambda gt: cal_beta(trans_gt(gt), *sheet.row_values(row)[2:4]))
                                    code_rs = f'{code}_{rs}'
                                    if code_rs not in outcomes_rs.columns:
                                        outcomes_rs[code_rs] = outcome
                                        if code in outcomes_code.columns:
                                            outcomes_code[code] = outcomes_code[code] * outcome
                                        else:
                                            outcomes_code[code] = outcome
                                else:
                                    if code not in outcomes_code.columns:
                                        outcomes_code[code] = pd.Series('Normal', index=gt_data.columns)
                                    temp = list(map(lambda trait, gt: cal_trait(trait, gt, *sheet.row_values(row)[2:5]),
                                        outcomes_code[code], gt_data.loc[rs]))
                                    outcomes_code[code] = temp
            outcomes_rs.to_csv(fr)
            fr.seek(0)
            fr.write('DONE')
        for file in quan_file:
            if not is_excel(file):
                code = os.path.basename(file).split('.')[0]
                res = get_ref_prs(file, vcf_files)
                res = res.rename(code)
                outcomes_code = outcomes_code.merge(res.apply(lambda x: exp(x)), left_index=True, right_index=True)
        outcomes_code.to_csv(fc)
        fc.seek(0)
        fc.write('DONE')
    return outcomes_rs, outcomes_code


def arg(args):
    # todo: Update
    try:
        opts, temp = getopt.getopt(args, "hi:d:r:o:", ['help', 'input=', 'data-excel=', 'ref-dir=', 'output='])
    except getopt.GetoptError:
        print(description)
        sys_exit(2)
    else:
        ref_dir = r'.\1KG'
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                print(description)
                sys_exit()
            elif opt in ('-i', '--input'):
                finput = arg
            elif opt in ('-d', '--data-excel'):
                data_excel = arg
            elif opt in ('-r', '--ref-dir'):
                ref_dir = arg
            elif opt in ('-o', '--output'):
                output = arg
        try:
            return [finput, data_excel, ref_dir, output]
        except NameError:
            raise Exception('Not complete arguments!')


def get_gt(vcf_gt: str, ref: str, alt: str):
    # if sum(not i.isnumeric() for i in vcf_gt) != 1:
    #     logger.warning('Not biallelic alleles appear')
    for i in vcf_gt:
        if i == '.':
            return {'*'}
        if not i.isnumeric():
            gt_n = vcf_gt.split(i)
            return set(alt.split(',')[int(i) - 1] if int(i) else ref for i in gt_n)
    return set(alt.split(',')[int(vcf_gt) - 1] if int(vcf_gt) else ref)


@use_time('Load sample vcf')
def load_vcf(human: Human, data_snps: set):
    with open(human.vcf, 'r') as f:
        for line in f:
            line = line.strip()
            if line[:2] == '##':
                continue
            elif line[0] == '#':
                header = line[1:].split('\t')
                assert len(header) == 10, 'VCF file contain multiple samples!'
                assert list(map(header.index, ['ID', 'REF', 'ALT', 'FORMAT'])) == [2, 3, 4, 8], 'Not standard VCF'
            else:
                snp = line.split('\t')
                snp_id = snp[2]
                assert snp[8] == 'GT', 'Not standard VCF'
                if snp_id in data_snps:
                    human.set_gt(snp_id, get_gt(snp[9], snp[3], snp[4]))


@use_time('Load indicator data')
def load_ind(human: Human, excel_ind_file: str):
    workbook = open_workbook(excel_ind_file)
    for sheet in workbook.sheets():
        itype = sheet.name
        for row in range(1, sheet.nrows):
            if sheet.row_types(row)[0] != 1:
                logger.warning(f'Not standard code format, indicator '
                               f'({sheet.row_values(row)[columns["columns_ind"]["name"]]}) may work improperly. ')
                human.set_ind(f'{sheet.row_values(row)[columns["columns_ind"]["code"]]:.0f}',
                              **sub_col(sheet.row_values(row), columns['columns_ind'], code=False), itype=itype)
            else:
                human.set_ind(**sub_col(sheet.row_values(row), columns['columns_ind']), itype=itype)


def get_snp_list(*algorithm_files: str, excel_files=True, prs_files=False):
    snp_list = set([])
    snp_column = None
    for file in algorithm_files:
        try:
            if is_excel(file):
                if excel_files:
                    workbook = open_workbook(file)
                    for sheet in workbook.sheets():
                        for row in range(1, sheet.nrows):
                            snp_list.add(sheet.row_values(row)[1])
            elif prs_files:
                header = True
                nopen = open_gz(file)
                with nopen(file, 'rb') as f:
                    for line in f:
                        if header:
                            if not snp_column:
                                snp_column = line.split(GWAS_setting['sep'].encode()).index(
                                    GWAS_setting['SNP'].encode())
                            header = False
                            continue
                        else:
                            snp_list.add(line.split(GWAS_setting['sep'].encode())[snp_column].decode())
        except Exception:
            logger.warning(f'{file} has some proble when getting snp list from it.')
            raise
    return snp_list


def get_columns_num(file: str, need_columns: list, sep=GWAS_setting['sep']):
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


@use_time('Loading algorithm data and calculate the result')
def load_cal(human: Human, qual_file: list, quan_file: list, temp_dir: str):
    # todo(Sheng): Many codes are similar to get_ref_res, maybe it can be simplified
    for file in qual_file + quan_file:
        if is_excel(file):
            workbook = open_workbook(file)
            ftype = 'quan' if file in quan_file else 'qual'
            for sheet in workbook.sheets():
                for row in range(1, sheet.nrows):
                    if sheet.row_types(row)[0] != 1:
                        human.cal_ind(ftype, f'{sheet.row_values(row)[0]:.0f}',
                                      *select_list(sheet.row_values(row), columns['columns_' + ftype].values())[1:])
                    else:
                        human.cal_ind(ftype, *select_list(sheet.row_values(row), columns['columns_' + ftype].values()))
        else:
            code = os.path.basename(file).split('.')[0]
            columns_num = get_columns_num(file, [GWAS_setting[i] for i in ['SNP', 'EA', 'BETA']])
            file = trans_gz(file, temp_dir)
            run_plink_cmd(f"--vcf {human.vcf} --score {file} {' '.join([str(i) for i in columns_num])} header "
                          f"--q-score-range {os.path.join('ref_res', 'need_range_list')}"
                          f" {os.path.join('ref_res', code + '.SNP.pvalue')} "
                          f"--out {os.path.join(temp_dir, 'result_prs')}",
                          plink='plink', delete_log=False)
            res = get_prs_res(os.path.join(temp_dir, "result_prs." + str(GWAS_setting["p_threshold"]) + ".profile"))
            with open(os.path.join(temp_dir, 'result_prs.log'), 'r') as f:
                result_text = f.read()
            detail = result_text.split('\n\n')[2]
            human.get_prs_data(code, res, detail)


def get_ftype(result: dict):
    # type_list = OrderedDict()
    return {key: result[key][0]['type'] for key in result}


@use_time('Whole process')
def main(name: str, input_file: str, ind_file: str, qual_files: list, quan_files: list, ref: str, output: str):
    res_str = ''
    log_name = log_start()
    temp_dir = mkdtemp(suffix='pageant')
    human = Human(name, input_file, temp_dir)
    initial_res_dir(output)
    try:
        vcf_files = [f"{os.path.join(ref, file)}" for file in os.listdir(ref) if
                     'vcf' in file and 'idi' not in file and 'tbi' not in file]
        initial_ref_data(qual_files, quan_files, vcf_files)
        load_vcf(human, get_snp_list(*qual_files, *quan_files))
        load_ind(human, ind_file)
        load_cal(human, qual_files, quan_files, temp_dir)
        human.add_distribution(output)
        res_dict = human.export_res(output=output)
    except Exception as e:
        res_str += f'Error: {str(e)}, analysis falied.'
        res_dict = {}
        logger.error(f'Error: {str(e)}, analysis falied.', exc_info=True)
        return res_str
    else:
        res_str += f'Analysis runs successfully!'
        logger.info(res_str)
    finally:
        env = Environment(loader=FileSystemLoader(os.getcwd()))
        template = env.get_template('./bin/template.html')
        with open(log_name) as fl:
            log_text = fl.read()
        log_text = log_text.replace('\n', '<br>')
        t = template.render(human=human, res=res_dict, type=get_ftype(res_dict),
                            time=strftime('%Y-%m-%d %H:%M'), config=locals(),
                            log=log_text)
        copy_files(['go_top.jpg', 'no_pic.jpg'], 'bin', os.path.join(output, 'html_files', 'img'))
        copy('./bin/Setting.css', os.path.join(output, 'html_files'))
        with open(os.path.join(output, 'Report.html'), 'w', encoding="UTF-8") as fhtml:
            fhtml.write(t)
        rm_dir(temp_dir)
        return res_dict
        # return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
    # vcf_file = r'.\test.vcf'
    # ref = r'.\reference'
    vcf_file = r'.\test\sample.vcf.gz'
    ref = r'.\reference'
    ind_file = r'.\database\001_Traits.xlsx'
    qual_file = [r'.\database\002_Qualitative.xlsx']
    quan_file = [r'.\database\003_Quantitative.xlsx', r'.\database\R0002.tsv.gz', r'.\database\R0004.tsv.gz']
    vcf_files = [f"{os.path.join(ref, file)}" for file in os.listdir(ref) if
                 'vcf' in file and 'idi' not in file and 'tbi' not in file]
    # vcf_file = r'E:\Data\Database\1KGP\sample\sample_1.vcf'
    # initial_ref_data(quan_file, vcf_files)
    # a = Human('User', vcf_file)
    # load_vcf(a, get_snp_list(*qual_file, *quan_file))
    # load_ind(a, excel_file)
    # load_cal(a, qual_file, quan_file)
    # a.add_quan_dis()
    # res = a.export_res()
    # env = Environment(loader=FileSystemLoader(os.getcwd()))
    # template = env.get_template('template.html')
    # t = template.render(human=a, res=res, type=get_ftype(res)[0], undected=get_ftype(res)[1], time="2020-10-04 19:00")
    # with open('temp.html', 'w', encoding="UTF-8") as f:
    #     f.write(t)
    test = main('Test', vcf_file, ind_file, qual_file, quan_file, ref, output=r'.\res')
