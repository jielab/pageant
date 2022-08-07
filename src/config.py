import configparser
import os
import sys
from typing import List, Dict, Set, Tuple, Generator, Iterable, Optional, Callable
from typing.io import IO, TextIO, BinaryIO
from collections import OrderedDict
from platform import system as platform_check
from time import sleep


platform = platform_check()

if platform == 'Darwin':
    if hasattr(sys, 'frozen'):
        raw_dir = os.path.dirname(sys.executable)
    else:
        raw_dir = os.path.dirname(sys.argv[0])
        if not raw_dir:
            raw_dir = os.getcwd()
else:
    raw_dir = os.getcwd()

default_config = {
    'file': {'quan_data': './algorithm/Quan',
             'qual_data': './algorithm/Qual',
             'quan_ref': './population_genome/g1k_ref',
             'qual_ref': './population_genome/g1k_ref',
             'maf_ref': './population_genome/g1k_ps.vcf.gz',
             'ps_ref': './population_genome/g1k_ps.vcf.gz',
             'qr_key': './output/qr_code/doctor/key',
             'qr_snps': './personal_genome/fingerprint_snps.txt',
             'qr_dir': './output/qr_code',
             'qr_user': './output/qr_code/user_qr_code.png',
             'qr_give': './output/qr_code/snp_qr_code.png',
             'population_file': './population_genome/g1k_ps_samples.txt',
             'concord_ref': './personal_genome/concordance.vcf.gz',
             'ref_structure': './population_genome/g1k_ref/g1k_ref.vcf.gz',
             'ref_ethnicity': './population_genome/g1k_ref_samples.txt',
             'query_db': './algorithm/Query_db/Phewas_Catalog.txt',
             'description_suf': "desc.txt", 'need_suf': 'snps.ref', 'gwas_suf': ".tsv",
             'plink_dir': './bin'
             },
    'name': {'SNP': 'SNP', 'EA': 'EA', 'P': 'P', 'BETA': 'BETA', 'OR': 'OR',
             'population_col': 'population', 'population_id': 'IID',
             'database_SNP': 'snp', 'quan_pop': 'EUR'
             },
    'parameter': {'p_threshold': 5e-8, 'clump-p1': 5e-8, 'clump-r2': 0.1, 'clump-kb': 1000, 'clump-p2': 0.01,
                  'show_code': True, 'qc_maf': 0.01, 'qc_vmiss': 0.02, 'qc_smiss': 0.02, 'qc_hardy': 50,
                  'qc_het': 3, 'qc_male_F': 0.4, 'qc_female_F': 0.6, 'qc_pihat': 0.2, 'use_qc_ref': False,
                  'ps_prune': True
                  },
    'module': {'sample_qc': True, 'ref_dist': True, 'ref_qc': True,
               'query_database': True, 'pharmgkb': True, 'clinvar': True, 'qr_code': True
               },
    'columns': {'columns_quan': OrderedDict({'SNP': ['SNP', 'RS_ID'], 'EA': ['EA', 'EFFECT_ALLELE'],
                                             'OR': ['OR', 'ODDS_RATIO'],
                                             'BETA': ['BETA', 'EFFECT_SIZE', 'EFFECT']}),
                'columns_qual': OrderedDict(
                    {'SNP': ['SNP', 'RS_ID'], 'GENOTYPE': ['GENOTYPE', 'GENO', 'GT'],
                     'MATCHED': ['MATCHED'], 'UNMATCHED': ['UNMATCHED']}),
                'names_ind': {'name': 'Name', 'Description': 'Description', 'sex': 'Sex'}
                },
    'lan_lib': {'male': ['male', 'Male', '男', '男性', '1', 1],
                'female': ['female', 'Female', '女', '女性', '2', 2],
                'none': ['', '无', 'None', 'No', None],
                'if': ['条件', 'if', 'IF'],
                'and': ['并列', 'and', 'AND', '&', 'And']
                },
}
config = configparser.ConfigParser()


def set_default_config(config_txt=os.path.join(raw_dir, 'bin', 'config.ini')):
    save_config(default_config, config_txt)


def load_config(config_txt=os.path.join(raw_dir, 'bin', 'config.ini')):
    parser = configparser.ConfigParser()
    parser.read(config_txt, encoding='UTF-8')
    return parser


def save_config(sconfig: dict, config_txt=os.path.join(raw_dir, 'bin', 'config.ini')):
    parser = configparser.ConfigParser()
    parser.read_dict(sconfig)
    with open(config_txt, 'w', encoding='UTF-8') as f:
        parser.write(f)


def set_config(other_parameters: dict, config_txt=os.path.join(raw_dir, 'bin', 'config.ini'))\
        -> configparser.ConfigParser:
    global config
    config = load_config(config_txt)
    for key in other_parameters:
        have = False
        for section in config:
            if key in config[section]:
                if type(other_parameters[key]) == str:
                    config[section][key] = other_parameters[key]
                else:
                    config[section][key] = repr(other_parameters[key])
                have = True
                break
        if have:
            continue
        else:
            raise TypeError(f'{key} is an invalid keyword argument for pageant')
    save_config(config, config_txt)
    sleep(2)
    return config


def average_progress(all_progress=45):
    progresses = [10, 5, 10, 0, 10, 10, 5]
    return all_progress / sum([progresses[i] *
                               list(convert_dict(config['module']).values())[i] for i in range(len(progresses))])


def convert_dict(dict_like_obj: configparser.ConfigParser or dict or any, eval_trans: bool = True) -> dict:
    return {key: eval(dict_like_obj[key]) if eval_trans else dict_like_obj[key] for key in dict_like_obj}


def get_sys_font() -> str:
    if platform == 'Darwin':
        return "Baskerville"
    elif platform == 'Windows':
        return 'Cambria'
    elif platform == 'Linux':
        return '"Times New Roman"'
    else:
        raise OSError("Unsupported platform.")
