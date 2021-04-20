from jinja2 import Environment, FileSystemLoader
from time import strftime
from src.functions import *


ref_code = pd.DataFrame()
ref_rs = pd.DataFrame()
failed_snps = set()
failed_samples = set()
img_dir = ''


# Main function
def initial_res_dir(output: str) -> None:
    """
    Initial result directory
    :param output: output directory
    :return: A result directory which has created the needed folder
    """
    global img_dir
    img_dir = os.path.join(output, 'genetic_report', 'html_files', 'img')
    if not os.path.isdir(output):
        os.mkdir(output)
    if not os.path.isdir(os.path.join(output, 'log_files')):
        os.mkdir(os.path.join(output, 'log_files'))
    if not os.path.isdir(os.path.join(output, 'population_QC')):
        os.mkdir(os.path.join(output, 'population_QC'))
    if not os.path.isdir(os.path.join(output, 'genetic_report')):
        os.mkdir(os.path.join(output, 'genetic_report'))
    output = os.path.join(output, 'genetic_report', 'html_files')
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, 'dist_plot')
    if not os.path.isdir(output):
        os.mkdir(output)
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)


@progress_value(15)
@use_time('Initial reference data')
def initial_ref_data(data_dir: str, vcf_list: list, type_list: dict, output: str) -> None:
    """
    Get reference population data from existing data or generating newly
    :param data_dir: Data directory
    :param vcf_list: Vcf files in reference population directory
    :param type_list: A text file described the types of indicators which were in the database
    :return: Two dataframe: one recorded the population snps result, the other recorded the code result
    """
    global ref_rs, ref_code
    # if not os.path.isdir('ref_res'):
    #     os.mkdir('ref_res')
    assert vcf_list, 'There is no vcf file in population reference directory.'
    ref_rs, ref_code = get_ref_cal(data_dir, vcf_list, type_list, failed_snps, failed_samples, output)
    get_ref_freq(data_dir, vcf_list, output)


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
        # os.system('rm sex_impute.*')
        if file_type == 'vcf':
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.SNPSEX[0] == 1 else 'Female'
        else:
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.PEDSEX[0] == 1 else 'Female'


@progress_value(10)
@use_time('Load sample vcf')
def load_vcf(human: Human, data_snps=None):
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
                if data_snps:
                    if snp_id in data_snps:
                        human.set_gt(snp_id, get_gt(snp[9], snp[3], snp[4], snp_id))
                else:
                    human.set_gt(snp_id, get_gt(snp[9], snp[3], snp[4], snp_id))


@progress_value(10)
@use_time('Load indicator data')
def load_data(human: Human, data_dir: str, temp_dir: str, type_list: dict, output: str) -> None:
    columns, config = load_config()['columns'], load_config()['parameter']
    columns_name = {value: key for key, value in columns['names_ind'].items()}
    os.chdir(data_dir)
    try:
        for type_dir in tqdm(type_list):
            os.chdir(type_dir)
            for ind_dir in tqdm([item for item in sorted(os.listdir('.')) if os.path.isdir(item)]):
                os.chdir(ind_dir)
                try:
                    info = {columns_name[key] if key in columns_name else key: value
                            for key, value in load_txt(ind_dir + config['description_suf'])}
                    if 'sex' not in info:
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
                        prs_fun = partial(run_plink_cmd,
                                          cmd=f"--vcf {human.vcf} "
                                              f"--read-freq {os.path.join(output, 'population_QC', 'prs.ref.afreq')} "
                                              f"--score {file} {' '.join([str(i) for i in columns_num])} header"
                                              f" center "
                                              f"--q-score-range "
                                              f"{os.path.join(output, 'population_QC', 'need_range_list')}"
                                              f" {os.path.join(output, 'population_QC', ind_dir + '.SNP.pvalue')} "
                                              f"--out {os.path.join(temp_dir, 'result_prs')}",
                                          delete_log=False)
                        prs_fun(plink='plink2')
                        res = get_prs_res(
                            os.path.join(temp_dir, "result_prs." + str(config["p_threshold"]) + ".sscore"))
                        with open(os.path.join(temp_dir, 'result_prs.log'), 'r') as f:
                            result_text = f.read()
                        # detail = result_text.split('\n\n')[2].strip()
                        gwas_num = read_line_num(file, header=True)
                        valid_num = read_line_num(os.path.join(output, 'population_QC', ind_dir + '.SNP.pvalue'), True)
                        load_num = search(pattern_load_num, result_text)[0]
                        detail = ['Using "clump" and "score" function in plink to calculate PRS.<br>&nbsp;',
                                  f'There are {gwas_num} variants for this trait in raw GWAS data.<br>&nbsp;',
                                  f'After clumping, there are {valid_num} valid variants for PRS.<br>&nbsp;',
                                  f'In this sample, there are {load_num} loaded to calculate PRS.<br>&nbsp;']
                        # todo: What detail should be telled to the users?
                        human.get_prs_data(ind_dir, res, detail)
                    else:
                        human.cal_ind_from_txt(file_name, type_dir, ind_dir, ref_rs, type_list, columns)
                os.chdir('..')
            os.chdir('..')
    finally:
        os.chdir(raw_dir)


@progress_value(10 * average_progress())
@use_time('Sample QC')
def sample_qc(human: Human, vcf: str, temp_dir: str, maf_ref=None, pca_ref=None, concord_ref=None) -> list:
    res = [None, None, None]
    run_plink_cmd(f'--vcf {vcf} --make-just-bim --sort-vars --out {os.path.join(temp_dir, "sample_qc")}')
    chr_and_vep(os.path.join(temp_dir, "sample_qc.bim"), img_dir)
    if maf_ref:
        res[0] = get_maf(human, temp_dir, maf_ref, img_dir)
    if pca_ref:
        pca_data(human, temp_dir, pca_ref)
        res[1] = pca_plot(temp_dir, img_dir)
    if concord_ref:
        res[2] = check_concordence(human, temp_dir, concord_ref)
    return res


@progress_value(10 * average_progress())
@use_time('Reference data QC')
def ref_qc(temp_dir: str, vcf_list: list) -> None:
    assert len(vcf_list) == 1, 'Reference QC can only execute with single vcf'
    run_plink_cmd(f'--vcf {vcf_list[0]} --hardy --het --check-sex --missing --freq --test-missing --genome '
                  f'--out {os.path.join(temp_dir, "qc")}', plink='plink')
    global failed_snps
    global failed_samples
    qc_maf(failed_snps, temp_dir, img_dir)
    qc_vmiss(failed_snps, temp_dir, img_dir)
    qc_smiss(failed_samples, temp_dir, img_dir)
    qc_het(failed_samples, temp_dir, img_dir)
    qc_hardy(failed_snps, temp_dir, img_dir)
    qc_sex(failed_samples, temp_dir, img_dir)
    qc_relatedness(failed_samples, temp_dir, img_dir)


@use_time('Query database')
def query_database(human: Human, data_dir: str, pharmgkb=load_config()['module']['pharmgkb'],
                   clinvar=load_config()['module']['clinvar']) -> tuple:
    query_database_path = os.path.join(data_dir, 'query_database')
    if not os.path.isdir(query_database_path):
        os.mkdir(query_database_path)

    # clinvar
    if clinvar:
        clinvar_res = {}
        # get clinvar data
        get_clinvar_data(query_database_path)
        # load
        load_clinvar(os.path.join(query_database_path, 'clinvar.vcf.gz'), clinvar_res, human)
        load_clinvar(os.path.join(query_database_path, 'clinvar_papu.vcf.gz'), clinvar_res, human)
        clinvar_res = trans_clinvar_res(clinvar_res)
    else:
        clinvar_res = None

    # pharmgkb
    if pharmgkb:
        # get pharmgkb data
        get_pharmgkb_data(query_database_path)

        pharm_pheno = load_pharm(os.path.join(query_database_path, 'var_pheno_ann.tsv'), human)
        pharm_drug = load_pharm(os.path.join(query_database_path, 'var_drug_ann.tsv'), human)
        pharm_func = load_pharm(os.path.join(query_database_path, 'var_fa_ann.tsv'), human)
        pharm_res = (pharm_pheno, pharm_drug, pharm_func)
    else:
        pharm_res = None

    return clinvar_res, pharm_res


@progress_value(5 * average_progress())
@use_time('Produce QR code')
def produce_qr_code(human: Human, need_list: str) -> list:
    import qrcode
    if os.path.isfile(need_list):
        qr_code_snps = []
        with open(need_list) as f:
            for line in f:
                qr_code_snps.append(line.strip())
    else:
        if '\n' in need_list:
            qr_code_snps = [i.strip() for i in need_list.split('\n')]
        elif ' ' in need_list:
            qr_code_snps = [i.strip() for i in need_list.split(' ') if i]
        else:
            qr_code_snps = [need_list]
    qr_text = ''
    for snps in qr_code_snps:
        if snps in human.gt_data:
            qr_text += trans_gt(human.gt_data[snps], connector='')
        else:
            qr_text += 'NN'
    img = qrcode.make(qr_text)
    img.save(os.path.join(img_dir, 'qr_code.png'))
    return qr_code_snps


@progress_value(5 * average_progress())
@use_time('Add population distribution')
def add_distribution(human: Human, output: str) -> None:
    assert not ref_code.empty, "Can not get reference data"
    output = os.path.join(output, 'genetic_report', 'html_files', 'dist_plot')
    for ind in tqdm(human.ind_data.values()):
        ind.add_dist(output, ref_code)


@progress_value(5)
def export_html(human: Human, type_list: dict, config: dict, log_txt: str, output: str, **html_data):
    env = Environment(loader=FileSystemLoader(os.getcwd()))
    template = env.get_template('./bin/Template.html')
    with open(log_txt) as fl:
        log_text = fl.read()
    log_text = log_text.replace('\n', '<br>')
    t = template.render(risk_cal=risk_cal, human=human, type=type_list, config=config,
                        time=strftime('%Y-%m-%d %H:%M'), log=log_text,
                        **html_data)
    copy_files(['go_top.jpg', 'no_pic.jpg'], 'bin', os.path.join(output, 'genetic_report', 'html_files', 'img'))
    copy('./bin/Setting.css', os.path.join(output, 'genetic_report', 'html_files'))
    with open(os.path.join(output, 'genetic_report', 'Report.html'), 'w', encoding="UTF-8") as fhtml:
        fhtml.write(t)


def check_concordence(human: Human, temp_dir: str, concord_ref: str) -> str:
    from matplotlib_venn import venn2, venn2_unweighted
    from copy import copy as obj_copy
    concord_human = Human('Ref_human')
    temp_concord = os.path.join(temp_dir, 'concordance')
    os.mkdir(temp_concord)
    recode_and_sex_impute.__wrapped__.__wrapped__(concord_human, concord_ref, temp_concord)
    load_vcf.__wrapped__.__wrapped__(concord_human)

    raw_human = obj_copy(human)
    load_vcf.__wrapped__.__wrapped__(raw_human)

    raw_snps = set(raw_human.gt_data.keys())
    ref_snps = set(concord_human.gt_data.keys())
    intersection_snps = raw_snps & ref_snps
    intersection_len = len(intersection_snps)
    raw_diff = len(raw_snps - ref_snps)
    concord_diff = len(ref_snps - raw_snps)
    same_num = geno_diff(raw_human, concord_human, intersection_snps)
    raw_diff_g = raw_diff + intersection_len - same_num
    concord_diff_g = concord_diff + intersection_len - same_num
    ratio = raw_diff_g / concord_diff_g if concord_diff != 0 else 1
    ratio = ratio if ratio >= 1 else 1 / ratio
    venn_plot = venn2 if ratio < 6 else venn2_unweighted
    venn_plot([raw_diff_g, concord_diff_g, same_num],
              set_labels=('Current genotype file', 'Previous genotype file'),
              set_colors=('r', 'g'), alpha=0.5)
    plt.title("Venn diagram of genotype files' variants")
    plt.savefig(os.path.join(img_dir, 'QC_concordence.png'))
    plt.close()

    rm_dir(temp_concord)
    return f'In the current genotype file, {intersection_len}  of {intersection_len + raw_diff} ' \
           f'({intersection_len / (intersection_len + raw_diff) * 100:.2f}%) ' \
           f'variants were found in previous genotype file.\n' \
           f'In these common variants, {same_num} of {intersection_len} ({same_num / intersection_len * 100:.2f}%) ' \
           f'variants were identical.'
