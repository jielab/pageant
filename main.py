import os
import getopt
import sys
import pandas as pd
import xlrd
import logging
from collections import OrderedDict
import gzip
import re
from time import time, strftime
from multiprocessing import Pool
import tempfile
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from jinja2 import Environment, FileSystemLoader
from shutil import copy


description = "python perhaps.py -i <sample IID> -d <sample file location> -s <the chr and positions of SNPs>"

columns = {'columns_ind': OrderedDict({'code': 0, 'name': 2, 'sex': 3}),
           'columns_quan': OrderedDict({'code': 0, 'snp': 1, 'reference': 2, 'beta': 3}),
           'columns_qual': OrderedDict({'code': 0, 'snp': 1, 'reference': 2, 'interpretation': 3, 'reverse-inter': 4})}
lan_lib = {'sex': {'male': ['male', '男', '男性'], 'female': ['female', '女', '女性']},
           None: ['', '无', 'None', 'No', None], 'if': ['条件', 'if']}

pattern = re.compile(rb'(?<=[\t])rs[0-9]*(?=[\t;])')
pat_header1 = re.compile(rb'^##')
pat_header2 = re.compile(rb'^#')
ref_code = pd.DataFrame()
ref_rs = pd.DataFrame()


if not os.path.isdir('log'):
    os.mkdir('log')
logger = logging.getLogger()
logger.setLevel(logging.INFO)
log_name = f'log/{strftime("%Y%m%d%H%M%S")}.log'
fh = logging.FileHandler(log_name)
formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)
logger.info('Log start.')


def use_time(process_name=None):
    def decorator(func):
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time()
            fun_res = func(*args, **kwargs)
            logger.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                        f' used time: {time() - start:.2f}s')
            return fun_res

        return wrapper
    return decorator


def initial_ref_data(quan_list: list, vcf_list: list):
    global ref_rs, ref_code
    if not os.path.isdir('./Ref_res'):
        os.mkdir('./Ref_res')
    ref_rs, ref_code = get_ref_cal(quan_list, vcf_list)


class Human(object):
    def __init__(self, name: str):
        self.name = name
        self.gt_data = {}
        self.ind_data = {}
        self.sex = None

    def set_gt(self, snp_id: str, gt: set):
        self.gt_data[snp_id] = gt

    def set_ind(self, code: str, name: str, sex: str, itype: str, **other):
        if trans_sex(sex) is None or trans_sex(sex) == self.sex:
            self.ind_data[code] = Ind(code, name, sex, itype, **other)

    def set_sex(self, sex: str or bool):
        self.sex = sex if type(sex) == bool else True if sex == 'Male' else False

    def cal_ind(self, ftype: str, code: str, snp: str, *cal_col):
        if code in self.ind_data:
            self.ind_data[code].set_ftype(ftype)
            if snp in self.gt_data:
                self.ind_data[code].cal(snp, self.gt_data[snp], *cal_col)
            else:
                self.ind_data[code].mis_snp(snp)
        else:
            logger.warning(f'Code ({code}) is missing in index excel!')

    @use_time('Add population distribution')
    def add_quan_dis(self, output):
        global ref_code
        assert not ref_code.empty, "Can not get reference data"
        for ind in self.ind_data.values():
            if ind.ftype == 'quan':
                ind.add_quan_dist(output)

    @use_time('Export report result')
    def export_res(self):
        result = {}
        for ind in self.ind_data.values():
            if not ind.status:
                if ind.outcome is None:
                    logger.warning(f'Code ({ind.code}) has no corresponding algorithm!')
                if 'Undected' not in result:
                    result['Undected'] = {}
                if ind.itype not in result['Undected']:
                    result['Undected'][ind.itype] = []
                result['Undected'][ind.itype].append(ind.name)
            else:
                if ind.itype not in result:
                    result[ind.itype] = []
                result[ind.itype].append(ind.report())
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
        self.quan_dist = None

    def set_ftype(self, ftype: str):
        if not self.ftype:
            if ftype == 'quan':
                self.outcome = 1
                self.quan_dist = {'Average': None, 'Percentage': None, 'Plot': None, 'Status': False}
            self.ftype = ftype

    def cal(self, snp: str, gt: set, ref: str, inter: str or float or int, no_inter=None):
        if self.ftype == 'qual':
            if not self.decide:
                if inter not in lan_lib['if']:
                    self.outcome = inter if gt == trans_gt(ref) else 'normal' if no_inter in lan_lib[None] else no_inter
                    if {snp: trans_gt(gt)} not in self.detail:
                        self.detail.append({snp: trans_gt(gt)})
                    if no_inter not in lan_lib[None]:
                        self.decide = None if not self.decide else True
                else:
                    self.decide = False if gt == trans_gt(ref) else None
        elif self.ftype == 'quan':
            self.outcome *= cal_beta(gt, ref, inter)
            if {snp: trans_gt(gt), 'beta': inter} not in self.detail:
                self.detail.append({snp: trans_gt(gt), 'beta': inter, 'outcome': cal_beta(gt, ref, inter)})
            else:
                logger.warning(f"There are duplicate snps ({snp}) when calculate the quantitative trait!")
        self.status = True

    def mis_snp(self, snp: str):
        warning = f'Snp ({snp}) can not find in vcf file!'
        if self.ftype == 'quan':
            global ref_rs
            if f'{self.code}_{snp}' in ref_rs.columns:
                beta = ref_rs[f"{self.code}_{snp}"].mean()
                warning += f' Using reference population average: {beta}'
                self.outcome *= beta
        logger.warning(warning)
        if warning not in self.detail:
            self.detail.append(warning)

    def add_quan_dist(self, output):
        assert self.ftype == 'quan', 'Not quantitative trait!'
        global ref_code
        if self.code in ref_code.columns:
            self.quan_dist['Average'] = ref_code[self.code].mean()
            self.quan_dist['Percentage'] = sum(ref_code[self.code] < self.outcome) / ref_code.shape[0] * 100
            self.quan_dist['Plot'] = quan_dist_plot(self.code, self.outcome, ref_code[self.code],
                                                    self.quan_dist['Percentage'], output)
            self.quan_dist['Status'] = True

    def report(self):
        if self.ftype == 'quan':
            if self.quan_dist['Status']:
                del self.quan_dist['Status']
                return {'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail,
                        'Distribution': self.quan_dist, **self.other}
        return {'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail, **self.other}


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


def quan_dist_plot(code: str, value: float, ref_data: pd.Series, per: float, report_dir):
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    report_dir = os.path.join(report_dir, 'img')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
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
    return os.path.join(report_dir, f'{code}.png')


def cal_beta(gt: set, ref: str, beta: int or float):
    return 1 if ref not in gt else beta ** 2 if len(gt) == 1 else beta


def trans_sex(sex: str or bool):
    if type(sex) == str:
        return True if sex in lan_lib['sex']['male'] else False if sex in lan_lib['sex']['female'] else None
    elif type(sex) == bool:
        return 'male' if sex else 'female'


def trans_gt(gt: set or str):
    if type(gt) == set:
        return '/'.join(list(gt)[0] * 2 if len(gt) == 1 else list(gt))
    elif type(gt) == str:
        for i in gt:
            if not i.isalpha():
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def select_list(ob_list: list, index):
    # assert len(ob_list) < max(index), 'list index out of range'
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def cal_sha(file: str):
    from hashlib import sha256
    with open(file, 'rb') as f:
        return sha256(f.read()).hexdigest()


@use_time("Calculate files' sha value")
def cal_multi_sha(file_list: str or list):
    from hashlib import sha256
    value = sha256(b'')
    if type(file_list) == str:
        file_list = os.listdir('dir')
    for file in file_list:
        if file.split('.')[-1] == 'gz' or file.split('.')[-1] == 'vcf':
            with open(file, 'rb') as f:
                value.update(f.read())
    return value.hexdigest()


def open_vcf(vcf_file: str):
    if vcf_file.split('.')[-1] == 'gz':
        return gzip.open
    elif vcf_file.split('.')[-1] == 'vcf':
        return open
    else:
        logger.warning('Unsupported file.')
        return


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
        nopen = open_vcf(vcf_file)
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
def extract_snp(vcf_file: str, temp_dir):
    nopen = open_vcf(vcf_file)
    need_lines = []
    need_snp_list = snp_list_rw('r')
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
def merge_files(file_dir: str, out: str):
    with open(out, 'ab') as fw:
        for file in os.listdir(file_dir):
            with open(os.path.join(file_dir, file), 'rb') as fr:
                for line in fr:
                    if line:
                        fw.write(line)


def get_gt_iid(vcf_list: list) -> pd.DataFrame:
    for vcf_file in vcf_list:
        nopen = open_vcf(vcf_file)
        if nopen:
            with nopen(vcf_file, 'rb') as f:
                global pat_header1
                global pat_header2
                for line in f:
                    if pat_header2.search(line):
                        if not pat_header1.search(line):
                            return pd.DataFrame(columns=[iid.decode() for iid in line.strip().split(b'\t')[9:]])


@use_time('Get reference GT data')
def get_ref_gt(vcf_list: list, excel_files: list,
               filter_vcf=os.path.join('Ref_res', 'filter.vcf')) -> pd.DataFrame:
    if os.path.isfile(filter_vcf):
        if verify_data(filter_vcf, vcf_list, 'multi'):
            get_ref_vcf(vcf_list, excel_files)
    else:
        get_ref_vcf(vcf_list, excel_files)
    gt_df = get_gt_iid(vcf_list)
    with open(filter_vcf, 'r') as f:
        for line in f:
            temp = line.strip().split('\t')
            gt_df.loc[temp[2]] = [get_gt(gt, temp[3], temp[4]) for gt in temp[9:]]
    return gt_df


def snp_list_rw(method: str, need_snp_list=None):
    import pickle
    with open('need_snp_list', method + 'b') as f:
        if method == 'w':
            pickle.dump(need_snp_list, f)
        elif method == 'r':
            return pickle.load(f)


def rm_dir(rmdir: str):
    # danger function, only can delete one level directory
    if os.listdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
    os.rmdir(rmdir)


@use_time('Get reference subset vcf')
def get_ref_vcf(vcf_files: list, excel_files: list):
    from functools import partial
    need_snp_list = get_snp_list(*excel_files)
    temp_dir = tempfile.mkdtemp(suffix='vcf')
    snp_list_rw('w', need_snp_list)  # temp
    func = partial(extract_snp, temp_dir=temp_dir)
    pool = Pool(processes=os.cpu_count() - 2)  # cpu_count need warp
    pool.map(func, vcf_files)
    pool.close()
    pool.join()
    merge_files(temp_dir, os.path.join('Ref_res', 'filter.vcf'))
    rm_dir(temp_dir)


@use_time()
def verify_data(data_file: str, source_files: str or list, mode: str):
    sha_cal = cal_sha if mode == 'single' else cal_multi_sha
    with open(data_file) as f:
        if f.read(4) == 'DONE':
            if f.read(61) == sha_cal(source_files)[3:]:
                return True
    return False


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


@use_time('Get reference quantitative result data')
def get_ref_cal(quan_excel: list, vcf_files: list, data_path='Ref_res'):
    ref_code_file = os.path.join(data_path, 'population_code_res.csv')
    ref_rs_file = os.path.join(data_path, 'population_rs_res.csv')
    if os.path.isfile(ref_code_file) and os.path.isfile(ref_rs_file):
        if verify_data(ref_code_file, vcf_files + quan_excel, 'multi'):
            if same_status([ref_code_file, ref_rs_file]):
                return pd.read_csv(ref_rs_file, skiprows=1, index_col=0), \
                       pd.read_csv(ref_code_file, skiprows=1, index_col=0)
    outcomes_rs = pd.DataFrame()
    outcomes_code = pd.DataFrame()
    gt_data = get_ref_gt(vcf_files, quan_excel)
    with open(ref_code_file, 'w', newline='') as fc:
        with open(ref_rs_file, 'w', newline='') as fr:
            sha = cal_multi_sha(vcf_files + quan_excel)
            fc.write(f'#{sha}\n')
            fr.write(f'#{sha}\n')
            for excel_file in quan_excel:
                workbook = xlrd.open_workbook(excel_file)
                for sheet in workbook.sheets():
                    for row in range(1, sheet.nrows):
                        code = sheet.row_values(row)[0] if sheet.row_types(row)[0] == 1 \
                            else f'{sheet.row_values(row)[0]:.0f}'
                        rs = sheet.row_values(row)[1]
                        if rs in gt_data.index:
                            outcome = gt_data.loc[rs].apply(
                                lambda gt: cal_beta(trans_gt(gt), *sheet.row_values(row)[2:4]))
                            code_rs = f'{code}_{rs}'
                            if code_rs not in outcomes_rs.columns:
                                outcomes_rs[code_rs] = outcome
                                if code in outcomes_code.columns:
                                    outcomes_code[code] = outcomes_code[code] * outcome
                                else:
                                    outcomes_code[code] = outcome
            outcomes_rs.to_csv(fr)
            fr.seek(0)
            fr.write('DONE')
        outcomes_code.to_csv(fc)
        fc.seek(0)
        fc.write('DONE')
    return outcomes_rs, outcomes_code


def arg(args):
    try:
        opts, temp = getopt.getopt(args, "hi:d:r:o:", ['help', 'input=', 'data-excel=', 'ref-dir=', 'output='])
    except getopt.GetoptError:
        print(description)
        sys.exit(2)
    else:
        ref_dir = r'.\1KG'
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                print(description)
                sys.exit()
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
    if sum(not i.isnumeric() for i in vcf_gt) != 1:
        logger.warning('Not biallelic alleles appear')
    for i in vcf_gt:
        if not i.isnumeric():
            gt_n = vcf_gt.split(i)
            return set(alt.split(',')[int(i) - 1] if int(i) else ref for i in gt_n)
    return set(alt.split(',')[int(vcf_gt) - 1] if int(vcf_gt) else ref)


@use_time('Load sample vcf')
def load_vcf(human: Human, vcf: str, data_snps: set):
    with open(vcf, 'r') as f:
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
    workbook = xlrd.open_workbook(excel_ind_file)
    for sheet in workbook.sheets():
        itype = sheet.name
        for row in range(1, sheet.nrows):
            if sheet.row_types(row)[0] != 1:
                logger.warning(f'Not standard code format, indicator '
                                f'({sheet.row_values(row)[columns["columns_ind"]["name"]]}) may work improperly. ')
                human.set_ind(f'{sheet.row_values(row)[0]:.0f}',
                              *select_list(sheet.row_values(row), columns['columns_ind'].values())[1:], itype)
            else:
                human.set_ind(*select_list(sheet.row_values(row), columns['columns_ind'].values()), itype)


def get_snp_list(*excel_files: str):
    snp_list = set([])
    for file in excel_files:
        workbook = xlrd.open_workbook(file)
        for sheet in workbook.sheets():
            for row in range(1, sheet.nrows):
                snp_list.add(sheet.row_values(row)[1])
    return snp_list


@use_time('Loading calculation data')
def load_cal(human: Human, qual_file: list, quan_file: list):
    for file in qual_file + quan_file:
        try:
            workbook = xlrd.open_workbook(file)
        except xlrd.XLRDError:
            raise Exception("Unsupported file!")
        else:
            ftype = 'quan' if file in quan_file else 'qual'
            for sheet in workbook.sheets():
                for row in range(1, sheet.nrows):
                    if sheet.row_types(row)[0] != 1:
                        human.cal_ind(ftype, f'{sheet.row_values(row)[0]:.0f}',
                                      *select_list(sheet.row_values(row), columns['columns_' + ftype].values())[1:])
                    else:
                        human.cal_ind(ftype, *select_list(sheet.row_values(row), columns['columns_' + ftype].values()))


def get_ftype(result: dict):
    undect = False
    type_list = OrderedDict()
    for key in result:
        if key == 'Undected':
            undect = True
        else:
            type_list[key] = 'qual' if len(result[key][0]) == 3 else 'quan'
    return type_list, undect


@use_time('Whole process')
def main(name: str, vcf_file: str, ind_file: str, qual_files: list, quan_files: list, ref: str, output: str):
    res_str = ''
    human = Human(name)
    try:
        vcf_files = [f"{os.path.join(ref, file)}" for file in os.listdir(ref) if
                     'genotypes.vcf.gz' in file and 'idi' not in file and 'tbi' not in file]
        initial_ref_data(quan_files, vcf_files)
        load_vcf(human, vcf_file, get_snp_list(*quan_files, *qual_files))
        load_ind(human, ind_file)
        load_cal(human, qual_files, quan_files)
        human.add_quan_dis(output)
        res_dict = human.export_res()
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
        template = env.get_template('./BIN/template.html')
        with open(log_name) as fl:
            log_text = fl.read()
        log_text = log_text.replace('\n', '<br>')
        t = template.render(human=human, res=res_dict, type=get_ftype(res_dict)[0],
                            undected=get_ftype(res_dict)[1], time=strftime('%Y-%m-%d %H:%M'), config=locals(),
                            log=log_text)
        copy('./BIN/go_top.jpg', f'{output}/img')
        copy('./BIN/Setting.css', output)
        with open(os.path.join(output, 'Report.html'), 'w', encoding="UTF-8") as fhtml:
            fhtml.write(t)
        if os.path.isfile('./need_snp_list'):
            os.remove('./need_snp_list')
        return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
    vcf_file = r'.\test.vcf'
    ref = r'E:/Data/Database/1KGP/raw'
    excel_file = r'.\Database\001指标.xlsx'
    excel_file2 = [r'.\Database\002算法.xlsx']
    excel_file3 = [r'.\Database\003算法.xlsx']
    vcf_files = [f"{ref}\\{file}" for file in os.listdir(ref) if
                 'genotypes.vcf.gz' in file and 'idi' not in file and 'tbi' not in file]
    # vcf_file = r'E:\Data\Database\1KGP\sample\sample_1.vcf'
    # initial_ref_data(excel_file3, vcf_files)
    # a = Human('User')
    # load_vcf(a, vcf_file, get_snp_list(*excel_file2, *excel_file3))
    # load_ind(a, excel_file)
    # load_cal(a, excel_file2, excel_file3)
    # a.add_quan_dis()
    # res = a.export_res()
    # env = Environment(loader=FileSystemLoader(os.getcwd()))
    # template = env.get_template('template.html')
    # t = template.render(human=a, res=res, type=get_ftype(res)[0], undected=get_ftype(res)[1], time="2020-10-04 19:00")
    # with open('temp.html', 'w', encoding="UTF-8") as f:
    #     f.write(t)
    main('Test', vcf_file, excel_file, excel_file2, excel_file3, ref, output='.')
