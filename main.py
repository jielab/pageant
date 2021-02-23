import getopt
from sys import exit as sys_exit
from src.log import *
from src.modules import *


description = "Usage: python perhaps.py -i INPUT_FILE [--ftype FILE_TYPE] -o OUTPUT_DIR [--trait TRAIT_FILE]\n" \
              "[--qual QUALITATIVE_FILE] [--quan QUANTITATIVE] [--ref-dir REFERENCE_DIR] [--file-dir REFERENCE_DIR]\n" \
              "[--sep GWAS_SEPARATOR] [-p GWAS_PTHRESHOLD] [--ref-dir REFERENCE_DIR] [--GWAS_config KEY=VALUE ...]"

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
def main(name: str, input_file: str, data_dir: str, ref: str, output: str, qr_snps_txt=None, **other_parameters):
    restart_progress_bar()
    if other_parameters:
        set_config(other_parameters)
    res_str = ''
    basic_res_dict = {}
    extra_res = []
    success = True
    log_name = log_start()
    temp_dir = mkdtemp(suffix='pageant')
    type_list = get_type_list(data_dir)
    human = Human(name)
    initial_res_dir(output)
    module = load_config()['module']
    try:
        recode_and_sex_impute(human, input_file, temp_dir)
        if module['sample_qc']:
            human.sample_qc(temp_dir)
            sample_qc(human.vcf, temp_dir)
        vcf_files = [os.path.abspath(f"{os.path.join(ref, file)}") for file in os.listdir(ref) if
                     'vcf' in file and 'idi' not in file and 'tbi' not in file]
        if module['ref_qc']:
            ref_qc(temp_dir, vcf_files)
        initial_ref_data(data_dir, vcf_files, type_list)
        load_vcf(human, get_snp_list(data_dir))
        load_data(human, data_dir, temp_dir, type_list)
        if module['ref_dist']:
            add_distribution(human, output)
        if module['query_database']:
            extra_res.append(query_database(human, data_dir))
        if module['qr_code']:
            extra_res.append(produce_qr_code(human, qr_snps_txt))
        basic_res_dict = human.export_res(output=output)
    except Exception as e:
        success = False
        res_str += f'Error: {str(e)}, analysis falied.'
        logging.error(f'Error: {str(e)}, analysis falied.', exc_info=True)
        return res_str
    else:
        res_str += f'Analysis runs successfully!'
        logging.info(res_str)
    finally:
        export_html(human, type_list, locals(), log_name, output, basic_res=basic_res_dict, module=module,
                    extra_res=extra_res, success=success)
        rm_dir(temp_dir)
        # return human, basic_res_dict, extra_res
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
    snp_list = './test/fingerprint_snps.txt'
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
    test = main('Test', vcf_file, data_dir, ref, output=r'.\res', qr_snps_txt=snp_list, file_type=file_type)
