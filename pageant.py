import getopt
import argparse
import warnings
from sys import exit as sys_exit
from src.log import *
from src.modules import *


version = '2021-09-10'
description = "Usage: python pageant.py -n --name NAME -i --input INPUT_FILE -o --output OUTPUT_DIR\n" \
              "\t Options [-c --config CONFIG_FILE] [-s --set-config KEY=VALUE ...]"
warnings.filterwarnings('ignore')
# todo: add description for function

# stream log start
if __name__ == '__main__':
    stream_log_start()


def get_kwargs(kwargs_list: Optional[List[str]]) -> dict:
    if kwargs_list:
        return dict([kwarg.split('=') for kwarg in kwargs_list])
    else:
        return {}


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
            extra_res.append(sample_qc(human, output, temp_dir))
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
    # paras = arg(sys.argv[1:])
    parser = argparse.ArgumentParser(prog='PAGEANT', description=f'PAGEANT ({version})')
    parser.add_argument('-v', '--version', action='version', version=version)
    sub_parser = parser.add_subparsers(dest='fun', help='The functions in PAGEANT')
    pageant = sub_parser.add_parser('main', help='RUN the whole analysis of PAGEANT')
    pageant.add_argument('-n', '--name', help="User's name", required=True)
    pageant.add_argument('-i', '--input', help="Input genotype file", required=True)
    pageant.add_argument('-o', '--output', help="Output directory", required=True)
    pageant.add_argument('-c', '--config', help="PATH of config file", default='./bin/config.ini')
    pageant.add_argument('-s', '--set-config', dest='kwargs', nargs='+', help="KEY=VALUE Set the specific config")

    umap = sub_parser.add_parser('umap', help='RUN the UMAP and PCA plot analysis')
    umap.add_argument('-s', '--sample', default=None,
                      help='The sample genotype data (VCF (*.vcf *.vcf.gz); 23andme (*.txt); '
                      'PLINK 1 binary (*.bed)); PLINK 2 binary (*.pgen)')
    umap.add_argument('-p', '--population', required=True, dest='population',
                      help='The population genotype data (vcf (*.vcf *.vcf.gz); '
                      'PLINK 1 binary (*.bed); PLINK 2 binary (*.pgen))')
    umap.add_argument('-m', '--metadata', default=None, help='The metadata for population file')
    umap.add_argument('-o', '--output', default='.', help='The directory where the images save')
    umap.add_argument('-t', '--temp', default=None, help='The temporary directory')
    umap.add_argument('--plink-dir', default=None, dest='plink',
                      help='The directory for storing PLINK executable files')
    umap.add_argument('-w', '--workpath', default=None, help='The working path of analysis')
    umap.add_argument('--pop-col', default='population', dest='pop',
                      help="The name of the column that describe population infomation")
    umap.add_argument('--id-col', default='IID', dest='iid',
                      help="The name of the column that describe individual id")
    umap.add_argument('--prune', dest='prune', action='store_true')
    umap.add_argument('--no-prune', dest='prune', action='store_false')
    umap.set_defaults(prune=True)

    add_rsid = sub_parser.add_parser('add_rsid', help='Add rsID in the GWAS file')
    add_rsid.add_argument('-i', '--input', help="Input GWAS file", required=True)
    add_rsid.add_argument('-w', '--workpath', default='./add_rsid', help="The working path of analysis")
    add_rsid.add_argument('-g', '--genome', default='19', help="The build version of genome")
    add_rsid.add_argument('-d', '--dbsnp', default='150', help="The version of dbSNP")
    add_rsid.add_argument('-o', '--output', default='sites-rsids.tsv', help="The output name of file")

    qr_code = sub_parser.add_parser('qr_code', help='Generate the SNP QR CODE')
    qr_code.add_argument('-s', '--snp-list', dest='snp', help='SNP list')
    qr_code.add_argument('-k', '--key', help='The file which contains public key')
    qr_code.add_argument('-o', '--output', default='./SNP_QR_CODE.png', help='The QR CODE file path')
    qr_code.add_argument('-l', '--logo', default='./bin/SNP_logo.png', help='The logo of the QR CODE')

    args = parser.parse_args()
    if args.fun == 'main':
        main(args.name, args.input, args.output, config_file=args.config, **get_kwargs(args.kwargs))
    elif args.fun == 'umap':
        if not args.temp:
            args.temp = mkdtemp(suffix='pageant')
        if not args.workpath:
            args.workpath = args.temp
        if not args.plink:
            args.plink = os.path.abspath('./bin')
        get_plink_dir(args.plink)
        ps_analyse(args.population, args.sample, args.metadata, args.temp, args.workpath, args.output, prune=args.prune,
                   pop_id=args.iid, pop_col=args.pop)
    elif args.fun == 'qr_code':
        import src.qr_code as crypto
        crypto.request(args.key, args.snp, args.output, args.logo)
    elif args.fun == 'add_rsid':
        import src.add_rsid as rsid
        rsid.output_name = args.output
        rsid.data_dir = args.workpath
        rsid.dbsnp_version = args.dbsnp
        rsid.genes_version = args.genome
        rsid.run(args.input)
    print('Finish!')
