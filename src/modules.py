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
    img_dir = os.path.join(output, 'html_files', 'img')
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, 'html_files')
    if not os.path.isdir(output):
        os.mkdir(output)
    output = os.path.join(output, 'dist_plot')
    if not os.path.isdir(output):
        os.mkdir(output)
    if not os.path.isdir(img_dir):
        os.mkdir(img_dir)


@progress_value(10)
@use_time('Initial reference data')
def initial_ref_data(data_dir: str, vcf_list: list, type_list: dict) -> None:
    """
    Get reference population data from existing data or generating newly
    :param data_dir: Data directory
    :param vcf_list: Vcf files in reference population directory
    :return: Two dataframe: one recorded the population snps result, the other recorded the code result
    """
    global ref_rs, ref_code
    if not os.path.isdir('ref_res'):
        os.mkdir('ref_res')
    assert vcf_list, 'There is no vcf file in population reference directory.'
    ref_rs, ref_code = get_ref_cal(data_dir, vcf_list, type_list, failed_snps, failed_samples)


@progress_value(5)
@use_time()
def recode_and_sex_impute(human: Human, file: str, temp_dir: str):
    """
    recode the sample, determine the sex of the sample
    :param human: Objects Human
    :param file: the sample path
    :param temp_dir: the temp temporary directory path
    :param file_type: the type of sample file
    :return: a converted sample file in the temporary directory and a sex result
    """
    file_type = load_config()['parameter']['file_type']
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
        human.vcf = os.path.join(temp_dir, 'temp.vcf')
    else:
        run_plink_cmd(f"--vcf {os.path.join(temp_dir, 'temp.vcf')} "
                      "--recode vcf "
                      "--rm-dup force-first "
                      f"--out {os.path.join(temp_dir, 'temp')}")
        # os.system('rm sex_impute.*')
        if file_type == 'vcf':
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.SNPSEX[0] == 1 else 'Female'
        elif file_type == '23andme':
            human.vcf, human.sex = os.path.join(temp_dir, 'temp.vcf'), 'Male' if data.PEDSEX[0] == 1 else 'Female'


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


@progress_value(10)
@use_time('Load indicator data')
def load_data(human: Human, data_dir: str, temp_dir: str, type_list: dict) -> None:
    columns, config = load_config()['columns'], load_config()['parameter']
    columns_name = {value: key for key, value in columns['names_ind'].items()}
    os.chdir(data_dir)
    try:
        for type_dir in tqdm(type_list):
            os.chdir(type_dir)
            for ind_dir in tqdm([item for item in os.listdir('.') if os.path.isdir(item)]):
                os.chdir(ind_dir)
                sex = None
                try:
                    info = {columns_name[key] if key in columns_name else key: value
                            for key, value in load_txt(ind_dir + config['description_suf'])}
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
                        # detail = result_text.split('\n\n')[2].strip()
                        gwas_num = read_line_num(file, header=True)
                        valid_num = read_line_num(os.path.join(raw_dir, 'ref_res', ind_dir + '.SNP.pvalue'), True)
                        load_num = search(pattern_load_num, result_text)[0]
                        detail = ['Using "clump" and "score" function in plink to calculate PRS.<br>&nbsp;',
                                  f'There are {gwas_num} variants for this trait in raw GWAS data.<br>&nbsp;',
                                  f'After clumping, there are {valid_num} valid variants for PRS.<br>&nbsp;',
                                  f'In this sample, there are {load_num} loaded to calculate PRS.<br>&nbsp;']
                        # todo: What detail should be telled to the users?
                        human.get_prs_data(ind_dir, res, detail)
                    else:
                        with open(file_name, encoding=config['encoding']) as f:
                            for line in f:
                                line = line.strip('\n\r')
                                if line:
                                    line = line.split(config['text_sep'])
                                    human.cal_ind(type_list[type_dir], ind_dir,
                                                  *select_list(line,
                                                               columns['columns_' + type_list[type_dir]].values()),
                                                  ref_rs=ref_rs)
                os.chdir('..')
            os.chdir('..')
    finally:
        os.chdir(raw_dir)


@progress_value(int(10 * average_progress()))
@use_time('Sample QC')
def sample_qc(vcf: str, temp_dir: str) -> None:
    run_plink_cmd(f'--vcf {vcf} --make-just-bim --sort-vars --out {os.path.join(temp_dir, "sample_qc")}')
    data = pd.read_csv(os.path.join(temp_dir, "sample_qc.bim"), '\t', skipinitialspace=True, header=None, dtype={0: str})
    data.columns = ['Chromosome', 'ID', 'Unknown', 'Position', 'Ref', 'ALT']

    pos_count = data.Chromosome.value_counts(sort=True)
    plt.figure(figsize=[9, 3.6], dpi=300)
    plt.bar(pos_count.index, pos_count)
    plt.grid(axis='y')
    plt.xlabel("Chromosome")
    plt.ylabel("Variations count")
    plt.savefig(os.path.join(img_dir, f'chromosome_count.png'))
    plt.close()
    if load_config()['module']['vep']:
        vep_res = vep_query(list(data.ID))
        vep_res = pd.DataFrame(vep_res).T
        vep_result = pd.merge(data.loc[:, 'ID'], vep_res, right_index=True, left_on='ID', how='left')
        vep_result = vep_result['Variant type'].fillna('Unknown').value_counts()
        arrow = False if vep_result['Variant type'].value_counts(normalize=True) > 0.80 else True
        pie_plot(vep_result, vep_result.index, os.paht.join(img_dir, './vep_res.png'), arrow=arrow)


@progress_value(int(10 * average_progress()))
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


@progress_value(int(5 * average_progress()))
@use_time('Produce QR code')
def produce_qr_code(human: Human, need_list: str) -> list:
    import qrcode
    qr_code_snps = []
    with open(need_list) as f:
        for line in f:
            qr_code_snps.append(line.strip())
    qr_text = ''
    for snps in qr_code_snps:
        if snps in human.gt_data:
            qr_text += trans_gt(human.gt_data[snps], connector='')
        else:
            qr_text += 'NN'
    img = qrcode.make(qr_text)
    img.save(os.path.join(img_dir, 'qr_code.png'))
    return qr_code_snps


@progress_value(int(5 * average_progress()))
@use_time('Add population distribution')
def add_distribution(human: Human, output: str) -> None:
    assert not ref_code.empty, "Can not get reference data"
    output = os.path.join(output, 'html_files', 'dist_plot')
    for ind in human.ind_data.values():
        ind.add_dist(output, ref_code)


@progress_value(5)
def export_html(human: Human, type_list: dict, config: dict, log_txt: str, output: str, **html_data):
    env = Environment(loader=FileSystemLoader(os.getcwd()))
    template = env.get_template('./bin/template.html')
    with open(log_txt) as fl:
        log_text = fl.read()
    log_text = log_text.replace('\n', '<br>')
    t = template.render(risk_cal=risk_cal, human=human, type=type_list, config=config,
                        time=strftime('%Y-%m-%d %H:%M'), log=log_text,
                        **html_data)
    copy_files(['go_top.jpg', 'no_pic.jpg'], 'bin', os.path.join(output, 'html_files', 'img'))
    copy('./bin/Setting.css', os.path.join(output, 'html_files'))
    with open(os.path.join(output, 'Report.html'), 'w', encoding="UTF-8") as fhtml:
        fhtml.write(t)
