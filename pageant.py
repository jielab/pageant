import getopt
from sys import exit as sys_exit
from src.log import *
from src.modules import *


version = '2021-08-10'
description = "Usage: python pageant.py -n --name NAME -i --input INPUT_FILE -o --output OUTPUT_DIR\n" \
              "\t Options [-c --config CONFIG_FILE] [-s --set-config KEY=VALUE ...]"
# todo: add description for function

# stream log start
if __name__ == '__main__':
    stream_log_start()


def get_kwargs(kwargs_str: str) -> dict:
    pairs = kwargs_str.split(' ')
    return dict([pair.split('=') for pair in pairs])


def arg(args: list) -> Tuple[List, Dict]:
    """
    Parses command line options and parameter list.
    :param args: Argument list to be parsed, without the leading reference to the running program.
    :return: Parameters needed for running pageant.
    """
    try:
        opts, temp = getopt.getopt(args, "hn:i:o:c:s:v:",
                                   ['help', 'name=', 'input=', 'output=', 'config=', 'set-config=', 'version'])
    except getopt.GetoptError:
        print(description)
        sys_exit(2)
    else:
        for opt, arg in opts:
            kwargs = {}
            fconfig = './bin/config.ini'
            if opt in ('-h', '--help'):
                print(description)
                sys_exit()
            elif opt in ('-v', '--version'):
                print(version)
                sys_exit()
            elif opt in ('-n', '--name'):
                name = arg
            elif opt in ('-i', '--input'):
                finput = arg
            elif opt in ('-o', '--output'):
                output = arg
            elif opt in ('-c', '--config'):
                fconfig = arg
            elif opt in ('-s', '--set-config'):
                kwargs = get_kwargs(arg)
        try:
            return [name, finput, output, fconfig], kwargs
        except NameError:
            raise Exception('Not complete arguments!')


@use_time('Whole process')
def main(name: str, input_file: str, output: str, config_file: str = './bin/config.ini', **kwargs) -> str:
    """
    Main process
    :param name: Name
    :param input_file: Input genotype file
    :param output: Output direcory
    :param config_file: configuration file
    :param kwargs: Other key parameters in the config file
    :return: A genetic report
    """
    restart_progress_bar()
    res_str = ''
    basic_res_dict = {}
    extra_res = []
    success = True
    temp_dir, type_list, ref_file_list, model, main_config = initial(output, config_file, kwargs)
    output = os.path.abspath(output)
    log_name = log_start(output)
    module = convert_dict(main_config['module'])

    human = Human(name)
    try:
        recode_and_sex_impute(human, input_file, temp_dir)
        if module['sample_qc']:
            human.sample_qc(temp_dir)
            extra_res.append(sample_qc(human, human.vcf, temp_dir))
        else:
            extra_res.append([None, None, None])

        if module['ref_qc']:
            ref_qc(temp_dir, ref_file_list)
        initial_ref_data(ref_file_list, type_list, output, model)
        load_vcf(human, get_snp_list(list(type_list)))
        load_data(human, temp_dir, type_list, output)
        if module['ref_dist']:
            add_distribution(human, output, len(ref_file_list) == 1)
        if module['query_database']:
            extra_res.append(query_database(human))
        if module['qr_code']:
            extra_res.append(produce_qr_code(human, output))
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
        export_html(human, {os.path.basename(itype): type_list[itype] for itype in type_list}, locals(),
                    log_name, output, basic_res=basic_res_dict, module=module, extra_res=extra_res,
                    one=len(ref_file_list) == 1, success=success)
        rm_dir(temp_dir)
        return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
    paras = arg(sys.argv[1:])
    main(*paras[0], **paras[1])
