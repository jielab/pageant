import getopt
from sys import exit as sys_exit
from src.log import *
from src.modules import *


description = "Usage: python perhaps.py -i INPUT_FILE [--ftype FILE_TYPE] -o OUTPUT_DIR [--trait TRAIT_FILE]\n" \
              "[--qual QUALITATIVE_FILE] [--quan QUANTITATIVE] [--ref-dir REFERENCE_DIR] [--file-dir REFERENCE_DIR]\n" \
              "[-p GWAS_PTHRESHOLD] [--ref-dir REFERENCE_DIR] [--GWAS_config KEY=VALUE ...]"

# stream log start
if __name__ == '__main__':
    stream_log_start()


def arg(args):
    # todo: Update cammand line version
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


@use_time('Whole process')
def main(name: str, input_file: str, data_dir: str, ref: str, output: str, qr_snps_txt: str or None = None,
         maf_ref: str or None = None, pca_ref: str or None = None, concord_ref: str or None = None,
         **kwargs) -> str:
    """
    Main process
    :param name: Name
    :param input_file: Input genotype file
    :param data_dir: Database directory
    :param ref: Reference directory
    :param output: Output direcory
    :param qr_snps_txt: SNPs list used for QC generating
    :param maf_ref: Reference genotype data for MAF QC
    :param pca_ref: Reference genotype data for PCA
    :param concord_ref: Reference genotype data for Concordance check
    :return:
    """
    restart_progress_bar()
    if kwargs:
        set_config(kwargs)
    res_str = ''
    basic_res_dict = {}
    extra_res = []
    success = True
    output = os.path.abspath(output)
    initial_res_dir(output)
    log_name = log_start(output)
    temp_dir = mkdtemp(suffix='pageant')
    type_list = get_type_list(data_dir)
    human = Human(name)
    module = load_config()['module']
    try:
        recode_and_sex_impute(human, input_file, temp_dir)
        if module['sample_qc']:
            human.sample_qc(temp_dir)
            extra_res.append(sample_qc(human, human.vcf, temp_dir, maf_ref, pca_ref, concord_ref))
        else:
            extra_res.append([None, None, None])
        vcf_files = [os.path.abspath(f"{os.path.join(ref, file)}") for file in os.listdir(ref) if
                     'vcf' in file and 'idi' not in file and 'tbi' not in file]
        if module['ref_qc']:
            ref_qc(temp_dir, vcf_files)
        initial_ref_data(data_dir, vcf_files, type_list, output)
        load_vcf(human, get_snp_list(data_dir))
        load_data(human, data_dir, temp_dir, type_list, output)
        if module['ref_dist']:
            add_distribution(human, output)
        if module['query_database']:
            extra_res.append(query_database(human, data_dir))
        if module['qr_code']:
            extra_res.append(produce_qr_code(human, qr_snps_txt))
        basic_res_dict = human.export_res(output=output)
    except Exception as e:
        success = False
        res_str += f'Error: {str(e)}, analysis falied. '
        logging.error(f'Error: {str(e)}, analysis falied.', exc_info=True)
        return res_str
    else:
        res_str += f'Analysis runs successfully! '
        logging.info(res_str)
    finally:
        export_html(human, type_list, locals(), log_name, output, basic_res=basic_res_dict, module=module,
                    extra_res=extra_res, success=success)
        rm_dir(temp_dir)
        # return human, basic_res_dict, extra_res
        return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
    # file_type = 'vcf'
    vcf_file = './personal_genome/sample.vcf.gz'
    # file_type = 'vcf'
    # vcf_file = '../../111-1642-2174.download.snp.txt'
    # vcf_file = 'D:/Users/Sheng/Documents/jie_huang.txt'
    # vcf_file = '../../111-1642-2174_202001221629.txt'
    file_type = 'vcf'
    ref = r'./population_genome'
    data_dir = './algorithm_database'
    vcf_files = [f"{os.path.join(ref, file)}" for file in os.listdir(ref) if
                 'vcf' in file and 'idi' not in file and 'tbi' not in file]
    snp_list = './personal_genome/fingerprint_snps.txt'
    maf_ref = pca_ref = './personal_genome/hapmap3.vcf.gz'
    concord_ref = './personal_genome/concordance.vcf.gz'
    test = main('Test', vcf_file, data_dir, ref, output=r'./output', qr_snps_txt=snp_list,
                # maf_ref=maf_ref,
                pca_ref=pca_ref,
                # concord_ref=concord_ref,
                file_type=file_type,
                qr_code=False, query_database=False, ref_qc=True
                )
