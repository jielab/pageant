import argparse
import warnings
from src.log import *
from src.modules import *


version = '2023-06-12'
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
    header = []
    success = True
    temp_dir, type_list, ref_file_list, model, main_config = initial(output, config_file, kwargs)
    output = os.path.abspath(output)
    log_name = log_start(output)
    module = convert_dict(main_config['module'])

    human = Human(name)
    try:
        recode_and_sex_impute(human, input_file, temp_dir)
        if module['sample_qc']:
            extra_res.append(sample_qc(human, output, temp_dir))
        else:
            extra_res.append([None, None, None])
        if module['ref_qc']:
            ref_qc(temp_dir, ref_file_list)
        initial_ref_data(ref_file_list, type_list, output, model)
        load_vcf(human, get_snp_list(list(type_list)))
        load_database(human, temp_dir, type_list, output)
        if module['ref_dist']:
            add_distribution(human, output, len(ref_file_list) == 1)
        if module['query_database']:
            extra_res.append(query_database(human))
        if module['qr_code']:
            with open(main_config['file']['qr_snps']) as f:
                snp_list = set(f.read().strip().split('\n'))
            load_vcf(human, snp_list)
            extra_res.append(produce_qr_code(human, output))
        basic_res_dict = human.export_res(output=output)
        header = add_header(basic_res_dict)
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
                    one=len(ref_file_list) == 1, success=success, header=header)
        rm_dir(temp_dir)
        return res_str + ' Report has saved in output directory.'


if __name__ == '__main__':
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
    umap.add_argument('--prune', dest='prune', action='store_true',
                      help='Using "prune" function before analysis')
    umap.add_argument('--no-prune', dest='prune', action='store_false',
                      help='Not using "prune" function before analysis')
    umap.set_defaults(prune=True)

    add_rsid = sub_parser.add_parser('add_rsid', help='Add rsID in the GWAS file')
    add_rsid.add_argument('-i', '--input', help="Input GWAS file", required=True)
    add_rsid.add_argument('-d', '--database', help="The dbSNP file", required=True)
    add_rsid.add_argument('-s', '--sep', default='\t', help="The delimiter of input file")
    add_rsid.add_argument('-c', '--core', default=None, help="The number of processors used in analysis")
    add_rsid.add_argument('--na', default='NA', help="The string of NA value")
    add_rsid.add_argument('--rsid', default='RS_ID', help="The name of annotated column")
    add_rsid.add_argument('--alt', default='ALT', help="The column's name of reference allele")
    add_rsid.add_argument('--ref', default='REF', help="The column's name of alternative allele")
    add_rsid.add_argument('--pos', default='POS', help="The column's name of allele position")
    add_rsid.add_argument('--chr', default='CHR', help="The column's name of chromosome")
    add_rsid.add_argument('-o', '--output', default='./annotated.tsv.gz', help="The output name of file", required=True)

    qr_code = sub_parser.add_parser('qr_code', help='Generate the SNP QR CODE')
    qr_code.add_argument('-s', '--snp-list', dest='snp', help='SNP list')
    qr_code.add_argument('-k', '--key', help='The file which contains public key')
    qr_code.add_argument('-o', '--output', default='./SNP_QR_CODE.png', help='The QR CODE file path')
    qr_code.add_argument('-l', '--logo', default='./bin/SNP_logo.png', help='The logo of the QR CODE')

    liftover = sub_parser.add_parser('liftover', help='Convert genome coordinates between different assemblies')
    liftover.add_argument('-i', '--input', help="Input GWAS file", required=True)
    liftover.add_argument('-c', '--chain-file', dest='chain_file', help="The chain file", required=True)
    liftover.add_argument('-s', '--sep', default='\t', help="The delimiter of input file")
    liftover.add_argument('-p', '--processor', default=None, help="The number of processors used in analysis")
    liftover.add_argument('--pos', default='POS', help="The column's name of allele position")
    liftover.add_argument('--chr', default='CHR', help="The column's name of chromosome")
    liftover.add_argument('--new-column', dest='new_column', action='store_true',
                          help='Store new genome coordinates in new columns')
    liftover.add_argument('--no-new-column', dest='new_column', action='store_false',
                          help='Replace old genome coordinates into new coordinates')
    liftover.set_defaults(new_column=True)
    liftover.add_argument('--chr-prefix', default='', help="The prefix of chromosome string")
    liftover.add_argument('-o', '--output', default='./liftover.tsv.gz', help="The output name of file", required=True)

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
        print('Finish!')
    elif args.fun == 'qr_code':
        import src.qr_code as crypto
        crypto.request(args.key, args.snp, args.output, args.logo)
        print('Finish!')
    elif args.fun == 'add_rsid':
        import src.add_rsid as add_rsid
        add_rsid.addrsid_run(args.input, args.database, args.output, args.core, args.sep, args.chr,
                             args.pos, args.ref, args.alt, args.rsid, args.na)
    elif args.fun == 'liftover':
        import src.liftover as liftover
        liftover.liftover_run(args.input, args.chain_file, args.output, args.processor, args.sep, args.chr, args.pos,
                              args.new_column, args.chr_prefix)
