import matplotlib.pyplot as plt
from numpy import cos, sin, arcsin, deg2rad, rad2deg, sign, log10
import pandas as pd
from matplotlib.ticker import FuncFormatter
import subprocess
from re import search
from shutil import copy
from src.decorator import *
from src.config import *


class Human(object):
    def __init__(self, name: str):
        self.name = name
        self.gt_data = {}
        self.ind_data = {}
        self.vcf, self.sex, self.missing = None, None, None

    def set_gt(self, snp_id: str, gt: set):
        self.gt_data[snp_id] = gt

    def set_ind(self, code: str, name: str, sex: str, itype: str, ind_dir: str, **other):
        if trans_sex(sex) is None or trans_sex(sex) == self.sex:
            self.ind_data[code] = Ind(code, name, sex, itype, ind_dir, **other)

    def cal_ind(self, ftype: str, code: str, snp: str, *cal_col, ref_rs: pd.DataFrame):
        if code in self.ind_data:
            self.ind_data[code].set_ftype(ftype)
            if snp in self.gt_data:
                self.ind_data[code].cal(snp, self.gt_data[snp], *cal_col)
            else:
                self.ind_data[code].mis_snp(snp, cal_col[1], ref_rs)

    def get_prs_data(self, code: str, res: float, detail: str):
        if code in self.ind_data:
            self.ind_data[code].add_prs_res(res, detail)

    def sample_qc(self, temp_dir: str):
        run_plink_cmd(f'--vcf {self.vcf} --missing --out {os.path.join(temp_dir, "sample_qc")}')
        qc_data = pd.read_csv(os.path.join(temp_dir, 'sample_qc.smiss'), '\t', skipinitialspace=True)
        self.missing = qc_data.iloc[-1, 1:].to_list()
        self.missing[0] = f'{self.missing[0]:,d}'
        self.missing[1] = f'{self.missing[1]:,d}'
        self.missing[2] = f'{self.missing[2] * 100:.2f}%'

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
        self.name = name if load_config()['parameter']['show_code'] == 0 else f'{name} ({code})'
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
        self.detail = detail
        self.status = True

    def cal(self, snp: str, gt: set, ref: str, inter: str or float or int, no_inter=None):
        lan_lib = load_config()['lan_lib']
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

    def mis_snp(self, snp: str, inter: str or float or int, ref_rs: pd.DataFrame):
        lan_lib = load_config()['lan_lib']
        warning = f'Snp ({snp}) can not find in vcf file!'
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
                warning += f' Code {self.code} algorithm can not process further.'
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
                self.distribution['Plot'] = quan_dist_plot(self.code, self.outcome, ref_code[self.code],
                                                           self.distribution['Low_p'], output)
                self.distribution['Status'] = True
            else:
                self.distribution['Plot'] = qual_dist_plot(self.code, self.name, output, ref_code)
                self.distribution['Status'] = True

    def report(self, output):
        if self.picture:
            copy(self.picture, os.path.join(output, 'html_files', 'img', os.path.basename(self.picture)))
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
    cmd = f'{os.path.join(raw_dir, "bin", plink)} ' + cmd
    a = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, sterr = a.communicate()
    a.wait()
    if delete_log:
        os.remove(out_file + '.log')
    if b'Error' in sterr:
        raise Exception("PlinkError:\n" + sterr.decode())


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


def trans_gt(gt: set or str, connector='/'):
    if type(gt) == set:
        return connector.join(list(gt)[0] * 2 if len(gt) == 1 else list(gt))  # Todo: not all biallelic
    elif type(gt) == str:
        for i in gt:
            if not i.isalpha() and i != '*':
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def gene_filp(raw: str):
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


def cal_beta(gt: set, ref: str, beta: int or float):
    if type(beta) == str:
        beta = float(beta)
    return 1 if ref not in gt else beta ** 2 if len(gt) == 1 else beta


def trans_sex(sex: str or bool):
    lan_lib = load_config()['lan_lib']
    if type(sex) == str:
        return True if sex in lan_lib['sex']['male'] else False if sex in lan_lib['sex']['female'] else None
    elif type(sex) == bool:
        return 'male' if sex else 'female'


def assign_pos(value: float, plt_size: tuple, left_offset=0.29, right_offset=0.04):
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


def pie_plot(count: list or pd.Series, labels: list or pd.Index, save_path: str, title=None, arrow=True):
    plt.figure(figsize=(8, 6), dpi=400)
    ax1, ax2 = plt.subplot(1, 3, (1, 2)), plt.subplot(1, 3, 3)
    patches, texts = ax1.pie(count, radius=1.2, pctdistance=1.13)
    ax1.axis('equal')
    ax2.axis('off')
    ax2.legend(patches, labels, loc='center')

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


def qual_dist_plot(code: str, name: str, report_dir: str, ref_code: pd.DataFrame):
    population = ref_code[code].value_counts()
    labels = population.index
    labels = [label.split(":")[0] for label in labels]
    pie_plot(population, labels, os.path.join(report_dir, f'{code}.png'),
             title=f'Distribution of {name} in reference population')
    return os.path.join('html_files', 'dist_plot', f'{code}.png')


def quan_dist_plot(code: str, value: float, ref_data: pd.Series, per: float, report_dir: str):
    plt.figure(dpi=400)
    plt.hist(log10(ref_data), bins=100, weights=[1. / len(ref_data)] * len(ref_data))
    plt.axvline(log10(value), c='r', lw=2, ls='--')
    plt.text(assign_pos(log10(value), plt.axis(), left_offset=0.33, right_offset=0.03), plt.axis()[3] * 0.816,
             f'Percentage: {per:.2f}% \nLevel: {risk_level(per)}',
             bbox={'facecolor': 'white', 'alpha': 0.8})
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda a, b: f'{a * 100:.1f}%'))
    plt.xlabel("Logarithm of Risk score")
    plt.ylabel("Percentage")
    plt.savefig(os.path.join(report_dir, f'{code}.png'), bbox_inches='tight')
    plt.close()
    return os.path.join('html_files', 'dist_plot', f'{code}.png')
