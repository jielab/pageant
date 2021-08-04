import subprocess
from re import search
from shutil import copy
import matplotlib.pyplot as plt
import pandas as pd
from numpy import cos, sin, arcsin, deg2rad, rad2deg, sign, log10, linspace
from scipy.stats import gaussian_kde
from src.decorator import *


plink_dir = ''
objects_config = configparser.ConfigParser()


def load_config_obj(main_config: configparser.ConfigParser):
    global objects_config
    objects_config = main_config


class PLINKERROR(Exception):
    def __init__(self, *args):
        super(PLINKERROR, self).__init__(*args)


class RUNPLINKERROR(Exception):
    def __init__(self, *args):
        super(RUNPLINKERROR, self).__init__(*args)


class Human(object):
    def __init__(self, name: str):
        self.name = name
        self.gt_data = {}
        self.ind_data = {}
        self.vcf, self.sex, self.missing, self.snps = None, None, None, None

    def set_gt(self, snp_id: str, gt: set):
        self.gt_data[snp_id] = gt

    def set_ind(self, code: str, name: str, sex: str, itype: str, ind_dir: str, **other):
        if trans_sex(sex) is None or trans_sex(sex) == trans_sex(self.sex):
            self.ind_data[code] = Ind(code, name, sex, itype, ind_dir, **other)

    def cal_ind(self, ftype: str, code: str, snp: str, *cal_col, ref_rs: pd.DataFrame):
        if code in self.ind_data:
            self.ind_data[code].set_ftype(ftype)
            if snp in self.gt_data:
                self.ind_data[code].cal(snp, self.gt_data[snp], *cal_col)
            else:
                self.ind_data[code].mis_snp(snp, cal_col[1], ref_rs)

    @auto_encoding
    @auto_sep
    def cal_ind_from_txt(self, file: str, type: str, ind_dir: str, ref_rs: pd.DataFrame, columns: dict, *,
                         encoding: str = 'UTF-8', sep: str = '\t') -> None:
        with open(file, encoding=encoding) as f:
            for line in f:
                line = line.strip('\n\r')
                if line:
                    line = line.split(sep)
                    assert len(line) > 1, sep.join(line)
                    self.cal_ind(type, ind_dir,
                                 *select_list(line, columns['columns_' + type].values()), ref_rs=ref_rs)

    def get_prs_data(self, code: str, res: float, detail: str):
        if code in self.ind_data:
            self.ind_data[code].add_prs_res(res, detail)

    def sample_qc(self, temp_dir: str):
        run_plink_cmd(f'--vcf {self.vcf} --missing --out {os.path.join(temp_dir, "sample_qc")}')
        qc_data = pd.read_csv(os.path.join(temp_dir, 'sample_qc.smiss'), '\t', skipinitialspace=True)
        self.snps = pd.read_csv(os.path.join(temp_dir, 'sample_qc.vmiss'), '\t', skipinitialspace=True)['ID'].to_list()
        self.missing = qc_data.iloc[-1, 1:].to_list()
        self.missing[0] = f'{self.missing[0]:,d}'
        self.missing[1] = f'{self.missing[1]:,d}'
        self.missing[2] = f'{self.missing[2] * 100:.2f}%'

    @progress_value(10)
    @use_time('Export report result')
    def export_res(self, output) -> Dict[str, list]:
        result = {}
        for ind in self.ind_data.values():
            if ind.ftype:
                if not ind.status:
                    if not ind.detail:
                        logging.warning(f'Code ({ind.code}) has no corresponding algorithm!')
                        ind.detail.append(f'Code ({ind.code}) has no corresponding algorithm!')
                    ind.outcome = 'Undetected'
                if ind.itype not in result:
                    result[ind.itype] = []
                result[ind.itype].append(ind.report(output))
        return result


class Ind(object):
    def __init__(self, code: str, name: str, sex: str, itype: str, ind_dir, **other):
        self.code = code
        self.name = name if eval(objects_config['parameter']['show_code']) else f'{name} ({code})'
        self.sex = sex
        self.itype = os.path.basename(itype)
        self.ftype = None
        self.outcome = None
        self.detail = []
        self.status = False
        self.decide = False
        self.distribution = None
        self.and_status = True
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
        self.detail = detail
        self.status = True

    def cal(self, snp: str, gt: set, ref: str, inter: str or float or int, no_inter=None):
        lan_lib = convert_dict(objects_config['lan_lib'])
        if self.ftype == 'qual':
            if not self.decide:
                if inter in lan_lib['if']:
                    # once in if condition, it cannot exit this status
                    # now it only supports one if condition once time
                    self.decide = False if equal_gt(gt, trans_gt(ref)) else \
                        False if equal_gt(gt, trans_gt(gene_filp(ref))) else None
                    if self.decide is False:
                        if {snp: trans_gt(gt)} not in self.detail:
                            self.detail.append({snp: trans_gt(gt)})
                else:
                    if self.decide is not None:
                        if inter in lan_lib['and']:
                            outcome = True if equal_gt(gt, trans_gt(ref)) else \
                                True if equal_gt(gt, trans_gt(gene_filp(ref))) else False
                            self.and_status = self.and_status & outcome
                        else:
                            if self.and_status:
                                self.outcome = inter if equal_gt(gt, trans_gt(ref)) else \
                                    inter if equal_gt(gt, trans_gt(gene_filp(ref))) else no_inter
                            else:
                                self.outcome = no_inter
                                self.and_status = True
                            self.outcome = 'normal' if self.outcome in lan_lib['none'] else self.outcome
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
                        self.outcome /= outcome
                        item['beta'] = inter
                        item['effect allele'] = ref
            if {snp: trans_gt(gt), 'effect allele': ref, 'beta': inter} not in self.detail:
                self.detail.append({snp: trans_gt(gt), 'effect allele': ref, 'beta': inter})
            self.status = True

    def mis_snp(self, snp: str, inter: str or float or int, ref_rs: pd.DataFrame):
        lan_lib = convert_dict(objects_config['lan_lib'])
        warning = f'Snp ({snp}) cannot be found in your genotype file!'
        if self.ftype == 'quan':
            if f'{self.code}_{snp}' in ref_rs.columns:
                beta = ref_rs[f"{self.code}_{snp}"].mean()
                warning += f' Using reference population average: {beta}'
                self.outcome *= beta
        elif self.ftype == 'qual':
            if inter in lan_lib['if'] + lan_lib['and']:
                self.decide = None
                if warning in self.detail:
                    self.detail.remove(warning)
                warning += f' Code {self.code} algorithm cannot process further.'
                self.status = False
        logging.warning(warning)
        if warning not in self.detail:
            self.detail.append(warning)

    def add_dist(self, output: str, ref_code: pd.DataFrame):
        if self.code in ref_code.columns:
            if self.ftype == 'quan':
                self.distribution['All'] = len(ref_code[self.code].dropna())
                self.distribution['Average'] = ref_code[self.code].mean()
                self.distribution['Low'] = sum(ref_code[self.code] < self.outcome)
                self.distribution['Up'] = self.distribution['All'] - self.distribution['Low']
                self.distribution['Low_p'] = self.distribution['Low'] / self.distribution['All'] * 100
                self.distribution['Up_p'] = self.distribution['Up'] / self.distribution['All'] * 100
                self.distribution['Plot'], self.distribution['SGRS'] = \
                    quan_dist_plot(self.code, self.outcome, ref_code[self.code], self.distribution['Low_p'], output)
                self.distribution['Status'] = True
            else:
                self.distribution['Plot'] = qual_dist_plot(self.code, self.name, output, ref_code)
                self.distribution['Status'] = True

    def report(self, output):
        if self.picture:
            copy(self.picture,
                 os.path.join(output, 'genetic_report', 'html_files', 'img', os.path.basename(self.picture)))
            self.picture = os.path.join('html_files', 'img', os.path.basename(self.picture))
        else:
            self.picture = os.path.join('html_files', 'img', 'no_pic.jpg')
        if self.status:
            if self.distribution['Status']:
                del self.distribution['Status']
                return {'type': self.ftype, 'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail,
                        'Distribution': self.distribution, 'Other': self.other, 'Status': self.status,
                        'Picture': self.picture}
        return {'type': self.ftype, 'Name': self.name, 'Outcome': self.outcome, 'Detail': self.detail,
                'Other': self.other, 'Status': self.status, 'Picture': self.picture}


# objects need functions
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


def get_plink_dir() -> None:
    global plink_dir
    plink_dir = os.path.abspath(objects_config['file']['plink_dir'])


def detect_plink(plink: str = 'plink' or 'plink2'):
    a = subprocess.Popen(f'which {plink}', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, sterr = a.communicate()
    a.wait()
    return True if stdout else False


def select_list(ob_list: list, index) -> list:
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def equal_gt(gt1: set, gt2: set) -> bool:
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


def trans_gt(gt: set or str, connector: str = '/') -> Set[str] or str:
    if type(gt) == set:
        return connector.join(list(gt)[0] * 2 if len(gt) == 1 else list(gt))  # Todo: not all biallelic
    elif type(gt) == str:
        for i in gt:
            if not i.isalpha() and i != '*':
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def gene_filp(raw: str) -> str:
    """
    Filp the genotype
    :param raw: raw genotype
    :return: filp genotype
    """
    # Todo: add on / off buttons
    base_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    if raw[0] + raw[-1] == '||':
        return raw[1:-1]
    else:
        return ''.join([base_dict[char] if char in base_dict.keys() else char for char in raw])


def cal_beta(gt: set, ref: str, beta: int or float) -> float:
    if type(beta) == str:
        beta = float(beta)
    return 1 if ref not in gt else beta ** 2 if len(gt) == 1 else beta


def trans_sex(sex: str or bool) -> bool or str:
    lan_lib = convert_dict(objects_config['lan_lib'])
    if type(sex) == str:
        return True if sex in lan_lib['male'] else False if sex in lan_lib['female'] else None
    elif type(sex) == bool:
        return 'male' if sex else 'female'


def assign_pos(value: float, plt_size: tuple, left_offset: float = 0.29, right_offset: float = 0.04):
    ran = plt_size[1] - plt_size[0]
    p = (value - plt_size[0]) / ran
    return value + right_offset * ran if p < 0.5 else value - left_offset * ran


def risk_level(per: float):
    assert 0 <= per <= 100, 'Error percentage!'
    cuts = [0, 1, 5, 10, 20, 40, 60, 80, 90, 95, 99, 100]
    for level, cut in enumerate(cuts):
        if per < cut:
            return level
    return 11


def pie_plot(count: list or pd.Series, labels: list or pd.Index, save_path: str, title: str or None = None,
             arrow: bool = True, legend_ratio: int = 3, **other_kwargs):
    plt.figure(figsize=(8, 6), dpi=400)
    ax1, ax2 = plt.subplot(1, legend_ratio, 1 if legend_ratio == 2 else tuple(range(1, legend_ratio))), \
               plt.subplot(1, legend_ratio, legend_ratio)
    patches, texts = ax1.pie(count, radius=1.2, pctdistance=1.13, **other_kwargs)
    ax1.axis('equal')
    ax2.axis('off')
    ax2.legend(patches, labels, loc='lower left', bbox_to_anchor=(0, 0.65))
    if title:
        plt.suptitle(title, size=12, weight="bold", y=0.95)

    if arrow:
        props = count.apply(lambda a: f'{a / sum(count) * 100:.2f}%')
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
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()


def qual_dist_plot(code: str, name: str, report_dir: str, ref_code: pd.DataFrame) -> str:
    population = ref_code[code].value_counts()
    labels = population.index
    labels = [label.split(":")[0] for label in labels]
    pie_plot(population, labels, os.path.join(report_dir, f'{code}.png'),
             title=f'Distribution of {name} in reference population')
    return os.path.join('html_files', 'dist_plot', f'{code}.png')


def quan_dist_plot(code: str, value: float, ref_data: pd.Series, per: float, report_dir: str) -> Tuple[str, float]:
    prs_data = log10(ref_data) if all(ref_data > 0) else ref_data
    ref_mean = prs_data.mean()
    ref_sd = prs_data.var() ** 0.5
    prs_data = prs_data.apply(lambda a: (a - ref_mean) / ref_sd)
    density = gaussian_kde(prs_data)
    density.covariance_factor = lambda: .25
    density._compute_covariance()
    left, right = min(prs_data), max(prs_data)
    xs = linspace(left, right, 300)
    own_value = ((log10(value) if all(ref_data > 0) else value) - ref_mean) / ref_sd

    plt.figure(dpi=400)
    if own_value > left and own_value < right:
        xs1 = linspace(left, own_value, 150)
        xs2 = linspace(own_value, right, 150)
        plt.fill_between(xs1, density(xs1), color='tab:blue', alpha=0.85)
        plt.fill_between(xs2, density(xs2), color='tab:blue', alpha=0.25)
    else:
        plt.fill_between(xs, density(xs), alpha=1)
    plt.plot(xs, density(xs), c='tab:blue', label='Reference population distribution')
    plt.axvline(own_value, c='r', lw=2, ls='--', label='Your standardized genetic risk score')
    plt.text(assign_pos(own_value, plt.axis(), left_offset=0.35, right_offset=0.03), plt.axis()[3] * 0.184,
             f'Percentage: {per:.2f}% '
             # f'\nLevel: {risk_level(per)}'
             , bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.xlabel("Standardized genetic rise score")
    plt.yticks([])
    plt.xlim(symmertry_range(plt.xlim()))
    plt.ylim(bottom=0, top=plt.ylim()[1] * 1.15)
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(report_dir, f'{code}.png'), bbox_inches='tight')
    plt.close()
    return os.path.join('html_files', 'dist_plot', f'{code}.png'), own_value


def symmertry_range(limit: Tuple[float, float], center: float = 0) -> Tuple[float, float]:
    if center - limit[0] > limit[1] - center:
        return limit[0], 2 * center - limit[0]
    else:
        return 2 * center - limit[1], limit[1]
