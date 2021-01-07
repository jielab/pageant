import os
import getopt
from sys import exit as sys_exit
import pandas as pd
from numpy import cos, sin, arcsin, deg2rad, rad2deg, sign, log10
from xlrd import open_workbook
import logging
from collections import OrderedDict
from gzip import open as gzip_open
from re import compile, sub, search
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
from tqdm import tqdm


description = "Usage: python perhaps.py -i INPUT_FILE [--ftype FILE_TYPE] -o OUTPUT_DIR [--trait TRAIT_FILE]\n" \
              "[--qual QUALITATIVE_FILE] [--quan QUANTITATIVE] [--ref-dir REFERENCE_DIR] [--file-dir REFERENCE_DIR]\n" \
              "[--sep GWAS_SEPARATOR] [-p GWAS_PTHRESHOLD] [--ref-dir REFERENCE_DIR] [--GWAS_config KEY=VALUE ...]"

columns = {'columns_ind': OrderedDict({'code': 0, 'name': 1, 'sex': 3, 'Disease_description': 5}),
           'columns_quan': OrderedDict({'snp': 0, 'reference': 1, 'beta': 2}),
           'columns_qual': OrderedDict({'snp': 0, 'reference': 1, 'interpretation': 2, 'reverse-inter': 3}),
           'names_ind': {'name': 'Name', 'Disease_description': 'Disease description', 'sex': 'Sex'}}
lan_lib = {'sex': {'male': ['male', 'Male', '男', '男性', '1', 1], 'female': ['female', 'Female', '女', '女性', '2', 2]},
           None: ['', '无', 'None', 'No', None], 'if': ['条件', 'if', 'IF'], 'and': ['并列', 'and', 'AND', '&', 'And']}
config = {'file_type': 'vcf', 'ref_structure':
          f'{os.path.abspath(os.path.join("reference", "ld_ref", "hapmap3.vcf.gz"))}',
          'logo_dir': f'{os.path.join("database", "logo")}', 'text_sep': '\t', 'encoding': 'UTF-8',
          'SNP': 'SNP', 'EA': 'EA', 'P': 'P', 'BETA': 'BETA', 'sep': '\t', 'p_threshold': 1e-5,
          'clump-p1': 1, 'clump-r2': 0.1, 'clump-kb': 250, 'clump-p2': 0.01
          }
type_list = {}

pattern = compile(rb'(?<=[\t])rs[0-9]*(?=[\t;])')
pat_header1 = compile(rb'^##')
pat_header2 = compile(rb'^#')
ref_code = pd.DataFrame()
ref_rs = pd.DataFrame()

progress_bar = 0

columns['names_ind'] = {value: key for key, value in columns['names_ind'].items()}
raw_dir = os.getcwd()

logger = logging.getLogger()
logger.setLevel(logging.INFO)
if __name__ == '__main__':
    ch = logging.StreamHandler()
    logger.addHandler(ch)


def log_start():
    if not os.path.isdir('log'):
        os.mkdir('log')
    log_name = f'log/{strftime("%Y%m%d%H%M%S")}.log'
    fh = logging.FileHandler(log_name)
    formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logging.info('Logging start.')
    return log_name


def gui_log_setup(hander: logging.Handler):
    logger.addHandler(hander)


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
            logging.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                         f' start')
            fun_res = func(*args, **kwargs)
            logging.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                         f' used time: {time() - start:.2f}s')
            return fun_res

        return wrapper

    return decorator


def progress_value(process_value: int):

    def decorator(func):
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            fun_res = func(*args, **kwargs)
            global progress_bar
            progress_bar += process_value
            logging.info(f"Progress of the analysis: {progress_bar}%")
            return fun_res

        return wrapper

    return decorator


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
    cmd = f'{os.path.join(raw_dir, "bin", plink)} ' + cmd
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


@progress_value(30)
@use_time('Initial reference data')
def initial_ref_data(data_dir: str, vcf_list: list) -> None:
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
    ref_rs, ref_code = get_ref_cal(data_dir, vcf_list)


class Human(object):
    def __init__(self, name: str, file: str, temp_dir: str):
        self.name = name
        self.gt_data = {}
        self.ind_data = {}
        self.vcf, self.sex = recode_and_sex_impute(file, temp_dir, config['file_type'])

    def set_gt(self, snp_id: str, gt: set):
        self.gt_data[snp_id] = gt

    def set_ind(self, code: str, name: str, sex: str, itype: str, ind_dir: str, **other):
        if trans_sex(sex) is None or trans_sex(sex) == self.sex:
            self.ind_data[code] = Ind(code, name, sex, itype, ind_dir, **other)

    def cal_ind(self, ftype: str, code: str, snp: str, *cal_col):
        if code in self.ind_data:
            self.ind_data[code].set_ftype(ftype)
            if snp in self.gt_data:
                self.ind_data[code].cal(snp, self.gt_data[snp], *cal_col)
            else:
                self.ind_data[code].mis_snp(snp, cal_col[1])

    def get_prs_data(self, code: str, res: float, detail: str):
        if code in self.ind_data:
            self.ind_data[code].add_prs_res(res, detail)

    @progress_value(10)
    @use_time('Add population distribution')
    def add_distribution(self, output: str):
        global ref_code
        assert not ref_code.empty, "Can not get reference data"
        output = os.path.join(output, 'html_files', 'dist_plot')
        for ind in self.ind_data.values():
            ind.add_dist(output)

    @progress_value(10)
    @use_time('Export report result')
    def export_res(self, output):
        result = {}
        for ind in self.ind_data.values():
            if ind.ftype:
                if not ind.status:
                    if not ind.detail:
                        logging.warning(f'Code ({ind.code}) has no corresponding algorithm!')
                        ind.detail.append(f'Code ({ind.code}) has no corresponding algorithm!')
                    ind.outcome = 'Undected'
                if ind.itype not in result:
                    result[ind.itype] = []
                result[ind.itype].append(ind.report(output))
        return result


class Ind(object):
    def __init__(self, code: str, name: str, sex: str, itype: str, ind_dir, **other):
        self.code = code
        self.name = name
        self.sex = sex
        self.itype = itype
        self.ftype = None
        self.outcome = None
        self.detail = []
        self.status = False
        self.decide = False
        self.distribution = None
        self.and_status = True
        # if 'ftype' in other:
        #     self.set_ftype(other['ftype'])
        #     other.pop('ftype')
        self.other = other
        try:
            self.picture = next(os.path.abspath(os.path.join(ind_dir, file)) for file in os.listdir(ind_dir)
                                if search(rf'^{self.code}\.(jpeg|jpg|png|bmp)$', file))
        except StopIteration:
            self.picture = None

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
                if inter in lan_lib['if']:
                    # once in if condition, it can not exit this status
                    # now it only supports one if condition once time
                    self.decide = False if equal_gt(gt, trans_gt(ref)) else \
                            False if equal_gt(gt, trans_gt(gene_filp(ref))) else None
                    if self.decide is False:
                        if {snp: trans_gt(gt)} not in self.detail:
                            self.detail.append({snp: trans_gt(gt)})
                else:
                    if self.decide is not None:
                        if inter in lan_lib['and']:
                            outcome = True if equal_gt(gt, trans_gt(ref)) else\
                                True if equal_gt(gt, trans_gt(gene_filp(ref))) else False
                            self.and_status = self.and_status & outcome
                        else:
                            if self.and_status:
                                self.outcome = inter if equal_gt(gt, trans_gt(ref)) else \
                                    inter if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
                            else:
                                self.outcome = no_inter
                                self.and_status = True
                            self.outcome = 'normal' if self.outcome in lan_lib[None] else self.outcome
                            if self.outcome != 'normal':
                                self.decide = True
                            self.status = True
                        if {snp: trans_gt(gt)} not in self.detail:
                            self.detail.append({snp: trans_gt(gt)})
        elif self.ftype == 'quan':
            outcome = cal_beta(gt, ref, inter)
            self.outcome *= outcome
            for item in self.detail:
                if type(item) == dict:
                    if snp in item.keys():
                        logging.warning(f"There are duplicate snps ({snp}) when calculate the quantitative trait!")
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
            if inter in lan_lib['if'] + lan_lib['and']:
                self.decide = None
                if warning in self.detail:
                    self.detail.remove(warning)
                warning += f' Code {self.code} algorithm can not process further.'
                self.status = False
        logging.warning(warning)
        if warning not in self.detail:
            self.detail.append(warning)

    def add_dist(self, output: str):
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
            copy(self.picture, os.path.join(output, 'html_files', 'img', os.path.basename(self.picture)))
        else:
            self.picture = os.path.join('html_files', 'img', 'no_pic.jpg')
        if self.status:
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
        x = cos(deg2rad(ang))
        y = sin(deg2rad(ang))
        y_text = y + y_dif
        y_text = (y_text - 0.8) * 0.33 + 0.8 if y_text > 0.8 else y_text
        horizontalalignment = "right" if ang <= 180 else "left"
        connectionstyle = f"angle, angleA={0 + rad2deg(arcsin(1.5 * y_dif))}, " \
                          f"angleB={ang + rad2deg(arcsin(1.5 * y_dif))}" if abs(y) > 0.5 else 'arc3'
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax1.annotate(props[i], xy=(1.18 * x, 1.18 * y), xytext=(1.5 * sign(x), round(1.5 * y_text, 2)),
                     horizontalalignment=horizontalalignment, **kw)

    plt.tight_layout()
    plt.savefig(os.path.join(report_dir, f'{code}.png'))
    plt.close()
    return os.path.join('html_files', 'dist_plot', f'{code}.png')


def quan_dist_plot(code: str, value: float, ref_data: pd.Series, per: float, report_dir: str):
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
    if type(beta) == str:
        beta = float(beta)
    return 1 if ref not in gt else beta ** 2 if len(gt) == 1 else beta


def cal_trait(trait: str or None, gt: set, ref: str, inter: str, no_inter: str or None):
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
            if not i.isalpha() and i != '*':
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def equal_gt(gt1: set, gt2: set):
    base = ('A', 'G', 'C', 'T')
    if gt1 == gt2:
        return True
    elif 'D' in gt1.union(gt2) or 'I' in gt1.union(gt2):
        if 'D' in gt2 or 'I' in gt2:
            gt1, gt2 = gt2, gt1
        gt2 = {i if set(i).difference(base) != set() else 'D' if len(i) == 1 else 'I' for i in gt2}
        return gt1 == gt2
    else:
        return False


def gene_filp(raw: str):
    # Todo: add on / off buttons
    base_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    if raw[0] + raw[-1] == '||':
        return raw[1:-1]
    else:
        return ''.join([base_dict[char] if char in base_dict.keys() else char for char in raw])


def select_list(ob_list: list, index):
    # assert len(ob_list) < max(index), 'list index out of range'
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def sub_col(ob_list: list, index: OrderedDict, code=True):
    if not code:
        index.pop('code')
    return {key: ob_list[index[key]] if index[key] < len(ob_list) else None for key in index}


@progress_value(10)
@use_time()
def recode_and_sex_impute(file: str, temp_dir: str, file_type: str):
    if file_type == 'vcf':
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
        logging.warning('Sex impute failed:\n' + f'{e.args}')
        run_plink_cmd(read_txt +
                      "--recode vcf "
                      f"--out {os.path.join(temp_dir, 'temp')}",
                      plink='plink')
        return os.path.join(temp_dir, 'temp.vcf'), None
    else:
        run_plink_cmd(f"--vcf {os.path.join(temp_dir, 'temp.vcf')} "
                      "--recode vcf "
                      "--rm-dup force-first "
                      f"--out {os.path.join(temp_dir, 'temp')}")
        # os.system('rm sex_impute.*')
        if file_type == 'vcf':
            return os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.SNPSEX[0] == 1 else 'Female'
        elif file_type == '23andme':
            return os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.PEDSEX[0] == 1 else 'Female'


def cal_sha(file: str):
    from hashlib import sha256
    with open(file, 'rb') as f:
        return sha256(f.read()).hexdigest()


def read_file_flow(obj: str or list, vcf_only: bool, re: str, include=True):
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


@use_time("Calculate files' sha")
def cal_multi_sha(file_list: str or list, vcf_only=True, re=r'', include=True):
    from hashlib import sha256
    value = sha256(b'')
    for file_flow in read_file_flow(file_list, vcf_only, re, include):
        value.update(file_flow)
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
def get_ref_gt(vcf_list: list, data_dir: str, filter_vcf=os.path.join('ref_res', 'filter.vcf')) -> pd.DataFrame:
    if os.path.isfile(filter_vcf):
        if not verify_data(filter_vcf, vcf_list, 'multi'):
            get_ref_vcf(vcf_list, get_snp_list(data_dir))
    else:
        get_ref_vcf(vcf_list, get_snp_list(data_dir))
    gt_df = get_gt_iid(vcf_list)
    with open(filter_vcf, 'r') as f:
        f.readline()
        for line in f:
            temp = line.strip().split('\t')
            try:
                gt_df.loc[temp[2]] = [get_gt(gt, temp[3], temp[4], warning=False) for gt in temp[9:]]
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


@use_time()
def verify_data(data_file: str, source_files: str or list, mode: str) -> bool:
    sha_cal = cal_sha if mode == 'single' else partial(cal_multi_sha, vcf_only=False, re=r'\.|_', include=False)
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
    for fid, file in enumerate(vcf_files):
        run_plink_cmd(f'--vcf {file} --score {prs_data} '
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


@use_time('Get reference quantitative result data')
def get_ref_cal(data_dir: str, vcf_files: list, data_path='ref_res'):
    ref_code_file = os.path.join(data_path, 'population_code_res.csv')
    ref_rs_file = os.path.join(data_path, 'population_rs_res.csv')
    ref_config = os.path.join(data_path, 'setting')
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
    with open(ref_code_file, 'w', newline='', encoding='UTF-8') as fc:
        with open(ref_rs_file, 'w', newline='', encoding='UTF-8') as fr:
            sha = cal_multi_sha(vcf_files + [data_dir] + [ref_config], vcf_only=False, re=r'\.|_', include=False)
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
                                res = get_ref_prs(file, vcf_files)
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


def get_type_list(data_dir: str):
    try:
        global type_list
        type_list = dict(load_txt(os.path.join(data_dir, 'Type.txt')))
    except FileNotFoundError:
        raise Exception("Can not find type list file.")


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


@progress_value(10)
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
                    human.set_gt(snp_id, get_gt(snp[9], snp[3], snp[4], snp_id))


@use_time('Load indicator data')
def load_ind_old(human: Human, excel_ind_file: str):
    workbook = open_workbook(excel_ind_file)
    for sheet in workbook.sheets():
        itype = sheet.name
        for row in range(1, sheet.nrows):
            if sheet.row_types(row)[0] != 1:
                logging.warning(f'Not standard code format, indicator '
                               f'({sheet.row_values(row)[columns["columns_ind"]["name"]]}) may work improperly. ')
                human.set_ind(f'{sheet.row_values(row)[columns["columns_ind"]["code"]]:.0f}',
                              **sub_col(sheet.row_values(row), columns['columns_ind'], code=False), itype=itype)
            else:
                human.set_ind(**sub_col(sheet.row_values(row), columns['columns_ind']), itype=itype)


def load_txt(description_text: str):
    with open(description_text, encoding=config['encoding']) as f:
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


@use_time('Load indicator data')
def load_data(human: Human, data_dir: str, temp_dir: str):
    os.chdir(data_dir)
    try:
        for type_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
            os.chdir(type_dir)
            for ind_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
                os.chdir(ind_dir)
                sex = None
                try:
                    info = {columns['names_ind'][key] if key in columns['names_ind'] else key: value
                            for key, value in load_txt(ind_dir + '_description')}
                    if sex not in info:
                        info['sex'] = ''
                    human.set_ind(code=ind_dir, **info, ind_dir=os.getcwd(), itype=type_dir)
                except NameError as e:
                    raise Exception(e.args[0] + f'in code {ind_dir}')
                except FileNotFoundError:
                    logging.warning(f'There is no description file in code {ind_dir}')
                exist, file_type, file_name = get_cal_file(ind_dir)
                if exist:
                    if file_type:
                        columns_num = get_columns_num(file_name, [config[i] for i in ['SNP', 'EA', 'BETA']])
                        file = trans_gz(file_name, temp_dir)
                        run_plink_cmd(
                            f"--vcf {human.vcf} --score {file} {' '.join([str(i) for i in columns_num])} header "
                            f"--q-score-range {os.path.join(raw_dir, 'ref_res', 'need_range_list')}"
                            f" {os.path.join(raw_dir, 'ref_res', ind_dir + '.SNP.pvalue')} "
                            f"--out {os.path.join(temp_dir, 'result_prs')}",
                            plink='plink', delete_log=False)
                        res = get_prs_res(
                            os.path.join(temp_dir, "result_prs." + str(config["p_threshold"]) + ".profile"))
                        with open(os.path.join(temp_dir, 'result_prs.log'), 'r') as f:
                            result_text = f.read()
                        detail = result_text.split('\n\n')[2]
                        human.get_prs_data(ind_dir, res, detail)
                    else:
                        with open(file_name, encoding=config['encoding']) as f:
                            for line in f:
                                line = line.strip('\n\r')
                                if line:
                                    line = line.split(config['text_sep'])
                                    human.cal_ind(type_list[type_dir], ind_dir,
                                                  *select_list(line,
                                                               columns['columns_' + type_list[type_dir]].values()))
                os.chdir('..')
            os.chdir('..')
    finally:
        os.chdir(raw_dir)


def get_snp_list_old(*algorithm_files: str, excel_files=True, prs_files=False, ):
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
                                snp_column = line.split(config['sep'].encode()).index(
                                    config['SNP'].encode())
                            header = False
                            continue
                        else:
                            snp_list.add(line.split(config['sep'].encode())[snp_column].decode())
        except Exception:
            logging.warning(f'{file} has some problem when getting snp list from it.')
            raise
    return snp_list


def get_cal_file(code: str):
    if os.path.isfile(code):
        return True, False, code
    elif os.path.isfile(code + '.tsv'):
        return True, True, code + '.tsv'
    elif os.path.isfile(code + '.tsv.gz'):
        return True, True, code + '.tsv.gz'
    else:
        return False, None, None


def get_snp_list(data_dir: str, cal_ind=True, prs_ind=False):
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


def get_columns_num(file: str, need_columns: list, sep=config['sep']):
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


@progress_value(30)
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
            columns_num = get_columns_num(file, [config[i] for i in ['SNP', 'EA', 'BETA']])
            file = trans_gz(file, temp_dir)
            run_plink_cmd(f"--vcf {human.vcf} --score {file} {' '.join([str(i) for i in columns_num])} header "
                          f"--q-score-range {os.path.join('ref_res', 'need_range_list')}"
                          f" {os.path.join('ref_res', code + '.SNP.pvalue')} "
                          f"--out {os.path.join(temp_dir, 'result_prs')}",
                          plink='plink', delete_log=False)
            res = get_prs_res(os.path.join(temp_dir, "result_prs." + str(config["p_threshold"]) + ".profile"))
            with open(os.path.join(temp_dir, 'result_prs.log'), 'r') as f:
                result_text = f.read()
            detail = result_text.split('\n\n')[2]
            human.get_prs_data(code, res, detail)


@use_time('Whole process')
def main(name: str, input_file: str, data_dir: str, ref: str, output: str, **other_parameters):
    global progress_bar
    progress_bar = 0
    if other_parameters:
        global config
        for key in other_parameters:
            config[key] = other_parameters[key]
    res_str = ''
    log_name = log_start()
    temp_dir = mkdtemp(suffix='pageant')
    get_type_list(data_dir)
    human = Human(name, input_file, temp_dir)
    initial_res_dir(output)
    try:
        vcf_files = [os.path.abspath(f"{os.path.join(ref, file)}") for file in os.listdir(ref) if
                     'vcf' in file and 'idi' not in file and 'tbi' not in file]
        initial_ref_data(data_dir, vcf_files)
        load_vcf(human, get_snp_list(data_dir))
        load_data(human, data_dir, temp_dir)
        human.add_distribution(output)
        res_dict = human.export_res(output=output)
    except Exception as e:
        res_str += f'Error: {str(e)}, analysis falied.'
        res_dict = {}
        logging.error(f'Error: {str(e)}, analysis falied.', exc_info=True)
        return res_str
    else:
        res_str += f'Analysis runs successfully!'
        logging.info(res_str)
    finally:
        env = Environment(loader=FileSystemLoader(os.getcwd()))
        template = env.get_template('./bin/template.html')
        with open(log_name) as fl:
            log_text = fl.read()
        log_text = log_text.replace('\n', '<br>')
        t = template.render(human=human, res=res_dict, type=type_list,
                            time=strftime('%Y-%m-%d %H:%M'), config=locals(),
                            log=log_text)
        copy_files(['go_top.jpg', 'no_pic.jpg'], 'bin', os.path.join(output, 'html_files', 'img'))
        copy('./bin/Setting.css', os.path.join(output, 'html_files'))
        with open(os.path.join(output, 'Report.html'), 'w', encoding="UTF-8") as fhtml:
            fhtml.write(t)
        rm_dir(temp_dir)
        # return human
        return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
    # vcf_file = r'.\test.vcf'
    # ref = r'.\reference'
    # vcf_file = r'./test/sample.vcf.gz'
    # file_type = 'vcf'
    # vcf_file = r'.\test\test.vcf.gz'
    # file_type = 'vcf'
    vcf_file = '../../111-1642-2174.download.snp.txt'
    file_type = '23andme'
    ref = r'.\reference'
    data_dir = './database'
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
    test = main('Test', vcf_file, data_dir, ref, output=r'.\res', file_type=file_type)
