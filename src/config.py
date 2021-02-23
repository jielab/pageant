from collections import OrderedDict
import os
import pickle
from time import sleep


raw_dir = os.getcwd()
default_config = config = {
    'columns': {'columns_ind': OrderedDict({'code': 0, 'name': 1, 'sex': 3, 'Disease_description': 5}),
                'columns_quan': OrderedDict({'snp': 0, 'reference': 1, 'beta': 2}),
                'columns_qual': OrderedDict(
                    {'snp': 0, 'reference': 1, 'interpretation': 2, 'reverse-inter': 3}),
                'names_ind': {'name': 'Name', 'Description': 'Description', 'sex': 'Sex'}},
    'lan_lib': {'sex': {'male': ['male', 'Male', '男', '男性', '1', 1],
                        'female': ['female', 'Female', '女', '女性', '2', 2]},
                None: ['', '无', 'None', 'No', None], 'if': ['条件', 'if', 'IF'],
                'and': ['并列', 'and', 'AND', '&', 'And']},
    'parameter': {'file_type': 'vcf', 'logo_dir': f'{os.path.join("database", "logo")}',
                  'ref_structure': f'{os.path.abspath(os.path.join("reference", "ld_ref", "hapmap3.vcf.gz"))}',
                  'text_sep': '\t', 'encoding': 'UTF-8',
                  'SNP': 'SNP', 'EA': 'EA', 'P': 'P', 'BETA': 'BETA', 'sep': '\t', 'p_threshold': 1e-5,
                  'clump-p1': 1, 'clump-r2': 0.1, 'clump-kb': 250, 'clump-p2': 0.01,
                  'show_code': 1, 'description_suf': ".desc.txt", 'need_suf': '.snps.ref', 'gwas_suf': ".tsv",
                  'qc_maf': 0.01, 'qc_vmiss': 0.02, 'qc_smiss': 0.02, 'qc_hardy': 50,
                  'qc_het': 3, 'qc_male_F': 0.4, 'qc_female_F': 0.6, 'qc_pihat': 0.2,
                  'use_qc_ref': True,
                  },
    'module': {
        'sample_qc': True, 'ref_dist': True, 'ref_qc': True, 'vep': False,
        'query_database': True, 'pharmgkb': True, 'clinvar': True,
        'qr_code': True
    }
    }


def set_default_config(config_txt=os.path.join(raw_dir, 'bin', 'config.pydat')):
    save_config(default_config, config_txt)


def load_config(config_txt=os.path.join(raw_dir, 'bin', 'config.pydat')):
    with open(config_txt, 'rb') as f:
        return pickle.load(f)


def save_config(config: dict, config_txt=os.path.join(raw_dir, 'bin', 'config.pydat')):
    with open(config_txt, 'wb') as f:
        pickle.dump(config, f)


def set_config(other_parameters: dict, config_txt=os.path.join(raw_dir, 'bin', 'config.pydat')):
    lconfig = load_config(config_txt)
    for key in other_parameters:
        if key in lconfig['parameter']:
            lconfig['parameter'][key] = other_parameters[key]
        elif key in lconfig['module']:
            lconfig['module'][key] = other_parameters[key]
        else:
            raise TypeError(f'{key} is an invalid keyword argument for pageant')
    save_config(lconfig, config_txt)
    sleep(2)


def average_progress(all_progress=50):
    progresses = [10, 5, 10, 10, 0, 10, 10, 5]
    return all_progress / sum([progresses[i] * list(default_config['module'].values())[i] for i in range(8)])
