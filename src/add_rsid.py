import csv
import gzip
import io
import itertools
import os
import random
import math
import typing as ty
from contextlib import contextmanager
from pathlib import Path
from typing import List, Dict, Callable, Union, Iterator, Any

import wget
from boltons.fileutils import mkdir_p, AtomicSaver

data_dir = os.path.join(os.getcwd(), 'add_rsid')
dbsnp_version = '150'
genes_version = '19'
output_name = 'sites-rsids.tsv'


def make_basedir(path: Union[str, Path]) -> None:
    mkdir_p(os.path.dirname(path))


def get_generated_path(*path_parts: str) -> str:
    path = os.path.join(data_dir, *path_parts)
    make_basedir(path)
    return path


def get_filepath(kind: str, *, must_exist: bool = True) -> str:
    if kind not in _single_filepaths:
        raise Exception(f"Unknown kind of filepath: {repr(kind)}")
    filepath: str = _single_filepaths[kind]()
    if must_exist and not os.path.exists(filepath):
        raise Exception(f"Filepath {filepath} of kind {kind} was requested but doesn't exist")
    return filepath


def get_tmp_path(arg: Union[Path, str]) -> str:
    if isinstance(arg, Path): arg = str(arg)
    if arg.startswith(get_generated_path()):
        mkdir_p(get_generated_path('tmp'))
        tmp_basename = arg[len(get_generated_path()):].lstrip(os.path.sep).replace(os.path.sep, '-')
        ret = get_generated_path('tmp', tmp_basename)
    elif arg.startswith(os.path.sep):
        ret = arg + '.tmp'
    else:
        mkdir_p(get_generated_path('tmp'))
        ret = get_generated_path('tmp', arg)
    assert ret != arg, (ret, arg)
    while os.path.exists(ret):
        ret = '{}/{}-{}'.format(os.path.dirname(ret), random.choice('123456789'), os.path.basename(ret))
    return ret


def scientific_int(s: str) -> int:
    """like int(s) but accepts "1.23e2" == 123"""
    try:
        return int(s)
    except ValueError:
        x = float(s)
        if x.is_integer():
            return int(x)
        raise Exception(f"invalid scientific_int: {s}")


null_values = ['', '.', 'NA', 'N/A', 'n/a', 'nan', '-nan', 'NaN', '-NaN', 'null', 'NULL']
default_field = {
    'aliases': [],
    'required': False,
    'type': str,
    'nullable': False,
    'from_assoc_files': True,
}


def round_sig(x: float, digits: int) -> float:
    if x == 0:
        return 0
    elif abs(x) == math.inf or math.isnan(x):
        raise ValueError("Cannot round infinity or NaN")
    else:
        log = math.log10(abs(x))
        digits_above_zero = int(math.floor(log))
        return round(x, digits - 1 - digits_above_zero)


class Field:
    def __init__(self, d):
        self._d = d

    def parse(self, value):
        """parse from input file"""
        # nullable
        if self._d['nullable'] and value in null_values:
            return ''
        # type
        x = self._d['type'](value)
        # range
        if 'range' in self._d:
            assert self._d['range'][0] is None or x >= self._d['range'][0]
            assert self._d['range'][1] is None or x <= self._d['range'][1]
        if 'sigfigs' in self._d:
            x = round_sig(x, self._d['sigfigs'])
        if 'proportion_sigfigs' in self._d:
            if 0 <= x < 0.5:
                x = round_sig(x, self._d['proportion_sigfigs'])
            elif 0.5 <= x <= 1:
                x = 1 - round_sig(1-x, self._d['proportion_sigfigs'])
            else:
                raise ValueError('cannot use proportion_sigfigs on a number outside [0-1]')
        if 'decimals' in self._d:
            x = round(x, self._d['decimals'])
        return x

    def read(self, value):
        """read from internal file"""
        if self._d['nullable'] and value == '':
            return ''
        return self._d['type'](value)


_single_filepaths: Dict[str, Callable[[], str]] = {
    'rsids-hg19': (lambda: get_generated_path(f'resources/rsids-v{dbsnp_version}-hg19.tsv.gz')),
    'rsids-hg38': (lambda: get_generated_path(f'resources/rsids-v{dbsnp_version}-hg38.tsv.gz')),
    'unanno': (lambda: get_generated_path('sites-unannotated.tsv.gz')),
    'sites-rsids': (lambda: get_generated_path(f'{output_name}.gz')),
    'sites': (lambda: get_generated_path('sites.tsv.gz')),
}

chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ['X', 'Y', 'MT']
chrom_order = {chrom: index for index, chrom in enumerate(chrom_order_list)}
chrom_aliases = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}


# Note; key order in these dicts is the order of columns in VariantFileWriter
per_variant_fields: Dict[str, Dict[str, Any]] = {
    'chrom': {
        'aliases': ['#CHROM', 'chr', 'CHR'],
        'required': True,
        'tooltip_underscoretemplate': '<b><%= d.chrom %>:<%= d.pos.toLocaleString() %> <%= d.ref %> / <%= d.alt %></b><br>',
        'tooltip_lztemplate': False,
    },
    'pos': {
        'aliases': ['BEG', 'BEGIN', 'BP'],
        'required': True,
        'type': scientific_int,
        'range': [0, None],
        'tooltip_underscoretemplate': False,
        'tooltip_lztemplate': False,
    },
    'ref': {
        'aliases': ['reference', 'REF'],
        'required': True,
        'tooltip_underscoretemplate': False,
        'tooltip_lztemplate': False,
    },
    'alt': {
        'aliases': ['alternate', "ALT"],
        'required': True,
        'tooltip_underscoretemplate': False,
        'tooltip_lztemplate': False,
    },
    'rsids': {
        'from_assoc_files': False,
        'tooltip_underscoretemplate': '<% _.each(_.filter((d.rsids||"").split(",")), function(rsid) { %>rsid: <b><%= rsid %></b><br><% }) %>',
        'tooltip_lztemplate': {'condition': 'rsid', 'template': '<strong>{{rsid}}</strong><br>'},
    },
    'nearest_genes': {
        'from_assoc_files': False,
        'tooltip_underscoretemplate': 'nearest gene<%= _.contains(d.nearest_genes, ",")? "s":"" %>: <b><%= d.nearest_genes %></b><br>',
        'tooltip_lztemplate': False,
    },
    'consequence': {
        'from_assoc_files': False,
    },
}

per_assoc_fields: Dict[str, Dict[str, Any]] = {
    'pval': {
        'aliases': ['PVALUE', 'P', 'P.VALUE'],
        'required': True,
        'type': float,
        'nullable': True,
        'range': [0, 1],
        'sigfigs': 2,
        'tooltip_lztemplate': {
            'condition': False,
            'template': ('{{#if pvalue|is_numeric}}P-value: <strong>{{pvalue|scinotation}}</strong><br>{{/if}}\n' +
                         '{{#if pval|is_numeric}}P-value: <strong>{{pval|scinotation}}</strong><br>{{/if}}'),
        },
        'display': 'P-value',
    },
    'beta': {
        'type': float,
        'nullable': True,
        'sigfigs': 2,
        'tooltip_underscoretemplate': 'Beta: <b><%= d.beta %></b><% if(_.has(d, "sebeta")){ %> (se:<b><%= d.sebeta %></b>)<% } %><br>',
        'tooltip_lztemplate': 'Beta: <strong>{{beta}}</strong>{{#if sebeta|is_numeric}} (se:<strong>{{sebeta}}</strong>){{/if}}<br>',
        'display': 'Beta',
    },
    'sebeta': {
        'aliases': ['se'],
        'type': float,
        'nullable': True,
        'sigfigs': 2,
        'tooltip_underscoretemplate': False,
        'tooltip_lztemplate': False,
    },
    'or': {
        'type': float,
        'nullable': True,
        'range': [0, None],
        'sigfigs': 2,
        'display': 'Odds Ratio',
    },
    'maf': {
        'type': float,
        'range': [0, 0.5],
        'sigfigs': 2,
        'tooltip_lztemplate': {'transform': '|percent'},
        'display': 'MAF',
    },
    'af': {
        'aliases': ['A1FREQ', 'FRQ'],
        'type': float,
        'range': [0, 1],
        'proportion_sigfigs': 2,
        'tooltip_lztemplate': {'transform': '|percent'},
        'display': 'AF',
    },
    'case_af': {
        'aliases': ['af.cases'],
        'type': float,
        'range': [0, 1],
        'proportion_sigfigs': 2,
        'tooltip_lztemplate': {'transform': '|percent'},
        'display': 'AF among cases',
    },
    'control_af': {
        'aliases': ['af.controls'],
        'type': float,
        'range': [0, 1],
        'proportion_sigfigs': 2,
        'tooltip_lztemplate': {'transform': '|percent'},
        'display': 'AF among controls',
    },
    'ac': {
        'type': float,
        'range': [0, None],
        'decimals': 1,
        'display': 'AC',
    },
    'r2': {
        'type': float,
        'proportion_sigfigs': 2,
        'nullable': True,
        'display': 'R2',
    },
    'tstat': {
        'type': float,
        'sigfigs': 2,
        'nullable': True,
        'display': 'Tstat',
    },
}

per_pheno_fields: Dict[str, Dict[str, Any]] = {
    'num_cases': {
        'aliases': ['NS.CASE', 'N_cases'],
        'type': int,
        'nullable': True,
        'range': [0, None],
        'display': '#cases',
    },
    'num_controls': {
        'aliases': ['NS.CTRL', 'N_controls'],
        'type': int,
        'nullable': True,
        'range': [0, None],
        'display': '#controls',
    },
    'num_samples': {
        'aliases': ['NS', 'N'],
        'type': int,
        'nullable': True,
        'range': [0, None],
        'display': '#samples',
    },
    # TODO: phenocode, phenostring, category, &c?
    # TODO: include `assoc_files` with {never_send: True}?
}




fields: Dict[str, Dict[str, Any]] = {**per_variant_fields, **per_assoc_fields, **per_pheno_fields}
reader_for_field: ty.Dict[str, ty.Callable[[str], ty.Any]] = {}
parser_for_field: ty.Dict[str, ty.Callable[[str], ty.Any]] = {}
for field_name, field_dict in fields.items():
    for k,v in default_field.items():
        field_dict.setdefault(k, v)
for field_name, field_dict in fields.items():
    obj = Field(field_dict)
    parser_for_field[field_name] = obj.parse
    reader_for_field[field_name] = obj.read

csv.register_dialect(
    'pheweb-internal-dialect',
    delimiter='\t',
    doublequote=False,
    escapechar='\\',
    lineterminator='\n',
    quotechar='"',
    skipinitialspace=False,
    strict=True,
)


def get_rsids_for_build(hg_build_number: int) -> None:
    dest_filepath = Path(get_filepath(f'rsids-hg{hg_build_number}', must_exist=False))
    if dest_filepath.exists():
        return

    # Download from https://resources.pheweb.org/
    url = f'https://resources.pheweb.org/{dest_filepath.name}'
    print(f'Downloading {dest_filepath} from {url}')
    try:
        wget.download(url=url, out=str(os.path.dirname(dest_filepath)))
        print()
    except Exception as exc:
        raise Exception(f'Failed to download rsids from {url}.') from exc


def get_rsid_reader(rsids_f: Iterator[str], rsids_filepath: str) -> Iterator[Dict[str, Any]]:
    prev_chrom_idx = -1
    prev_pos = -1
    for line in rsids_f:
        if not line.startswith('##'):
            if line.startswith('#'):
                assert line.rstrip('\r\n').split('\t') == '#CHROM POS ID REF ALT QUAL FILTER INFO'.split(), repr(line)
            else:
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 5:
                    raise ValueError(f'Line has wrong number of fields: {line} - {fields}')
                chrom, pos, rsid, ref, alt_group = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
                if chrom not in chrom_order:
                    try:
                        chrom = chrom_aliases[chrom]
                    except KeyError:
                        raise Exception((
                                f'The rsids file, {rsids_filepath}, contains the unknown chromsome {chrom}.\n' +
                                f'The recognized chromosomes are: {list(chrom_order.keys())}.\n' +
                                f'Recognized aliases are: {list(chrom_aliases.keys())}.\n'))
                chrom_idx = chrom_order[chrom]
                if prev_chrom_idx > chrom_idx:
                    raise Exception((
                            f'The rsids file, {rsids_filepath}, contains chromosomes in the wrong order.' +
                            f'The order should be: {chrom_order_list}' +
                            f'but instead {chrom_order_list[prev_chrom_idx]} came before {chrom_order_list[chrom_idx]}'))
                if prev_chrom_idx == chrom_idx and prev_pos > pos:
                    raise Exception(f'The rsids file, {rsids_filepath}, on chromosome {chrom_order_list[chrom_idx]}, '
                                    f'has position {prev_pos} before {pos}.')
                assert rsid.startswith('rs')
                # Sometimes the reference contains `N`, and that's okay.
                assert all(base in 'ATCGN' for base in ref), (chrom, pos, ref, alt_group)
                for alt in alt_group.split(','):
                    # Alt can be a comma-separated list
                    if alt == '.':
                        continue
                    assert all(base in 'ATCGN' for base in alt), (chrom, pos, ref, alt)
                    yield {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'rsid': rsid}


def get_one_chr_pos_at_a_time(iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
    """
    Turns
    [{'chr':'1', 'pos':123, 'ref':'A', 'alt':'T'}, {'chr':'1', 'pos':123, 'ref':'A', 'alt':'GC'},
    {'chr':'1', 'pos':128, 'ref':'A', 'alt':'T'},...]
    into:
    [ [{'chr':'1', 'pos':123, 'ref':'A', 'alt':'T'}, {'chr':'1', 'pos':123, 'ref':'A', 'alt':'GC'}] ,
    [{'chr':'1', 'pos':128, 'ref':'A', 'alt':'T'}] ,...]
    where variants with the same position are in a list.
    """
    for k, g in itertools.groupby(iterator, key=lambda cpra: (cpra['chrom'], cpra['pos'])):
        yield list(g)


def are_match(seq1: str, seq2: str) -> bool:
    """Compares nucleotide sequences.  Eg, "A" == "A", "A" == "N", "A" != "AN"."""
    if seq1 == seq2: return True
    if len(seq1) == len(seq2) and 'N' in seq1 or 'N' in seq2:
        return all(b1 == b2 or b1 == 'N' or b2 == 'N' for b1, b2 in zip(seq1, seq2))
    return False


def mtime(filepath: Union[str, Path]) -> float:
    return os.stat(filepath).st_mtime


@contextmanager
def read_gzip(filepath):  # mypy doesn't like it
    # hopefully faster than `gzip.open(filepath, 'rt')` -- TODO: find out whether it is
    with gzip.GzipFile(filepath, 'rb') as f:  # leave in binary mode (default), let TextIOWrapper decode
        with io.BufferedReader(f, buffer_size=2 ** 18) as g:  # 256KB buffer
            with io.TextIOWrapper(g) as h:  # bytes -> unicode
                yield h


@contextmanager
def read_maybe_gzip(filepath: Union[str, Path]):
    if isinstance(filepath, Path): filepath = str(filepath)
    is_gzip = False
    with open(filepath, 'rb', buffering=0) as raw_f:  # no need for buffers
        if raw_f.read(3) == b'\x1f\x8b\x08':
            is_gzip = True
    if is_gzip:
        with read_gzip(filepath) as f:
            yield f
    else:
        with open(filepath, 'rt', buffering=2 ** 18) as f:  # 256KB buffer
            yield f


@contextmanager
def VariantFileReader(filepath: Union[str, Path], only_per_variant_fields: bool = False):
    """
    Reads variants (as dictionaries) from an internal file.  Iterable.  Exposes `.fields`.

        with VariantFileReader('a.tsv') as reader:
            print(reader.fields)
            for variant in reader:
                print(variant)
    """
    with read_maybe_gzip(filepath) as f:
        reader: Iterator[List[str]] = csv.reader(f, dialect='pheweb-internal-dialect')
        try:
            fields = next(reader)
        except StopIteration:
            raise Exception(f"It looks like the file {filepath} is empty")
        # This won't happen in normal use but it's convenient for temporary internal re-routing
        if fields[0].startswith('#'):
            fields[0] = fields[0][1:]
        for field in fields:
            assert field in per_variant_fields or field in per_assoc_fields, field
        if only_per_variant_fields:
            yield _vfr_only_per_variant_fields(fields, reader)
        else:
            yield _vfr(fields, reader)


class _vfr:
    def __init__(self, fields: List[str], reader: Iterator[List[str]]):
        self.fields = fields
        self._reader = reader

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return self._get_variants()

    def _get_variants(self) -> Iterator[Dict[str, Any]]:
        parsers: List[Callable[[str], Any]] = [reader_for_field[field] for field in self.fields]
        for unparsed_variant in self._reader:
            assert len(unparsed_variant) == len(self.fields), (unparsed_variant, self.fields)
            variant = {field: parser(value) for parser, field, value in zip(parsers, self.fields, unparsed_variant)}
            yield variant


class _vfr_only_per_variant_fields:
    def __init__(self, fields: List[str], reader: Iterator[List[str]]):
        self._all_fields = fields
        self._extractors = [(reader_for_field[field], field, colidx) for colidx, field in enumerate(fields) if
                            field in per_variant_fields]
        self.fields = [e[1] for e in self._extractors]
        self._reader = reader

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return self._get_variants()

    def _get_variants(self) -> Iterator[Dict[str, Any]]:
        for unparsed_variant in self._reader:
            assert len(unparsed_variant) == len(self._all_fields), (unparsed_variant, self._all_fields)
            variant = {field: parser(unparsed_variant[colidx]) for parser, field, colidx in self._extractors}
            yield variant


class _vfw:
    def __init__(self, f, allow_extra_fields: bool, filepath: str):
        self._f = f
        self._allow_extra_fields = allow_extra_fields
        self._filepath = filepath

    def write(self, variant: Dict[str, Any]) -> None:
        if not hasattr(self, '_writer'):
            fields: List[str] = []
            for field in fields:
                if field in variant:
                    fields.append(field)
            extra_fields = list(set(variant.keys()) - set(fields))
            if extra_fields:
                if not self._allow_extra_fields:
                    raise Exception(f"ERROR: found unexpected fields {extra_fields} among the expected fields "
                                    f"{fields} while writing {self._filepath}.")
                fields += extra_fields
            self._writer = csv.DictWriter(self._f, fieldnames=fields, dialect='pheweb-internal-dialect')
            self._writer.writeheader()
        self._writer.writerow(variant)

    def write_all(self, variants: Iterator[Dict[str, Any]]) -> None:
        for v in variants:
            self.write(v)


@contextmanager
def VariantFileWriter(filepath: str, allow_extra_fields: bool = False, use_gzip: bool = True):
    '''
    Writes variants (represented by dictionaries) to an internal file.

        with VariantFileWriter('a.tsv') as writer:
            writer.write({'chrom': '2', 'pos': 47, ...})

    Each variant/association/hit/loci written must have a subset of the keys of the first one.
    '''
    part_file = get_tmp_path(filepath)
    make_basedir(filepath)
    if use_gzip:
        with AtomicSaver(filepath, text_mode=False, part_file=part_file, overwrite_part=True,
                         rm_part_on_exc=False) as f:
            with gzip.open(f, 'wt', compresslevel=2) as f_gzip:
                yield _vfw(f_gzip, allow_extra_fields, filepath)
    else:
        with AtomicSaver(filepath, text_mode=True, part_file=part_file, overwrite_part=True, rm_part_on_exc=False) as f:
            yield _vfw(f, allow_extra_fields, filepath)


def run(raw_gwas: str) -> None:
    in_filepath = raw_gwas
    out_filepath = get_filepath('sites-rsids', must_exist=False)
    rsids_filepath = get_filepath(f'rsids-hg{genes_version}', must_exist=False)

    if not os.path.exists(rsids_filepath):
        print('Fetching rsids...')
        get_rsids_for_build(int(genes_version))

    if os.path.exists(out_filepath) and max(mtime(in_filepath), mtime(rsids_filepath)) <= mtime(out_filepath):
        print('rsid annotation is up-to-date!')
        return

    with VariantFileReader(in_filepath) as in_reader, \
            read_maybe_gzip(rsids_filepath) as rsids_f, \
            VariantFileWriter(out_filepath, True) as writer:

        rsid_group_reader = get_one_chr_pos_at_a_time(get_rsid_reader(rsids_f, rsids_filepath))
        cp_group_reader = get_one_chr_pos_at_a_time(in_reader)

        rsid_group = next(rsid_group_reader)
        for cp_group in cp_group_reader:

            # Advance rsid_group until it is up to/past cp_group
            while True:
                if rsid_group[0]['chrom'] == cp_group[0]['chrom']:
                    rsid_is_not_behind = rsid_group[0]['pos'] >= cp_group[0]['pos']
                else:
                    rsid_is_not_behind = chrom_order[rsid_group[0]['chrom']] >= chrom_order[cp_group[0]['chrom']]
                if rsid_is_not_behind:
                    break
                else:
                    try:
                        rsid_group = next(rsid_group_reader)
                    except StopIteration:
                        break

            if rsid_group[0]['chrom'] == cp_group[0]['chrom'] and rsid_group[0]['pos'] == cp_group[0]['pos']:
                # we have rsids at this position!  will they match on ref/alt?
                for cpra in cp_group:
                    rsids = [rsid['rsid'] for rsid in rsid_group if
                             (cpra['ref'] == rsid['ref'] and are_match(cpra['alt'], rsid['alt'])) or
                             (cpra['ref'] == rsid['alt'] and are_match(cpra['alt'], rsid['ref']))]
                    # if len(rsids) > 1:
                    #     print('WARNING: the variant {chrom}-{pos}-{ref}-{alt} has multiple rsids:
                    #     {rsids}'.format(**cpra, rsids=rsids))
                    cpra['rsids'] = ','.join(rsids)
                    writer.write(cpra)
            else:
                # No match, just print each cpra with an empty `rsids` column
                for cpra in cp_group:
                    cpra['rsids'] = ''
                    writer.write(cpra)


if __name__ == '__main__':
    data_dir = os.path.abspath('../add_rsid')
    run(os.path.abspath('../add_rsid/test.tsv'))
