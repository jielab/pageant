import os
from time import time
from tempfile import mkdtemp
from itertools import compress
from gzip import open as gzip_open
from collections import OrderedDict
from multiprocessing import Pool, cpu_count
from typing import Iterable, List, Tuple, Optional, IO


columns_name_l = OrderedDict({'CHR': ['#CHROM', 'chr', 'CHR', 'chrom'],
                              'BP': ['BEG', 'BEGIN', 'pos', 'BP']})
chrom_aliases = {'23': 'X', '24': 'Y', '25': 'MT', 'M': 'MT'}
sep_str_l = '\t'
new_columns = True


def open_(file: str, mode: str = 'rb', **kwargs) -> IO:
    if file.split('.')[-1] == 'gz':
        return gzip_open(file, mode=mode, **kwargs)
    else:
        return open(file, mode=mode, **kwargs)


def rm_dir(rmdir: str) -> None:
    if os.listdir(rmdir):
        for file in os.listdir(rmdir):
            os.remove(os.path.join(rmdir, file))
    os.rmdir(rmdir)


def select_list(ob_list: list, index) -> list:
    return [ob_list[i] if i < len(ob_list) else None for i in index]


def get_chain(file: str) -> Iterable[List[str]]:
    chain = []
    with open_(file) as f:
        for line in f:
            if line.strip():
                chain.append(line.strip().decode())
            else:
                yield chain
                chain = []


def get_ind(pos: List[int], pos_length: List[int]) -> List[int or None]:
    idx_rev = []
    for p in pos:
        if p is None:
            idx_rev.append(None)
        else:
            i = 0
            while p >= pos_length[i]:
                p -= pos_length[i]
                i += 1
            idx_rev.append(i)
    return [idx_rev.index(i) if i in idx_rev else None for i in range(len(pos_length))]


def get_columns_index(col_names: List[str], need_cols: List[List[str]]) -> List[int or None]:
    # col_names = [name.upper().replace(' ', '_') for name in col_names]
    needs = [need for need_col in need_cols for need in need_col]
    need_pos = [len(need_col) for need_col in need_cols]
    need_idx = [needs.index(name) if name in needs else None for name in col_names]
    return get_ind(need_idx, need_pos)


def file_div(file: str, n_div: Optional[int] = None, start: int = 0, min_step: int = 1) -> Iterable[Tuple[int, int]]:
    if not n_div:
        n_div = cpu_count()
    with open_(file) as f:
        end = 0
        for _ in f:
            end += 1
    print(f'Total line of {file}: {end}')
    step = int(end / n_div + 1)
    step = step if step > min_step else min_step
    while start < end:
        if start + step >= end:
            yield start, end
        else:
            yield start, start + step
        start += step


def read_file_part(file: str, start: int, end: Optional[int] = None) -> Iterable[Tuple[int, bytes]]:
    if not end:
        end = float('inf')
    with open_(file) as f:
        i = 0
        for line in f:
            if i >= start:
                if i < end:
                    yield i, line
                else:
                    break
            i += 1


def format_chrome(raw: str, chr_prefix: str = '') -> str:
    if chr_prefix:
        body = raw.lstrip(chr_prefix)
        try:
            body = int(body)
        except ValueError:
            return f'chr{body}'
        else:
            if body <= 22:
                return f'chr{body}'
            else:
                return f'chr{chrom_aliases[str(body)]}'
    return f'chr{raw}'


def reformat_chrom(raw: str, chr_prefix: str = '') -> str:
    return raw.replace('chr', chr_prefix, 1)


class Seq:
    def __init__(self, raw: List[str]):
        header = raw[0].split(' ')
        self.ref_chrome = header[2]
        self.ref_size = int(header[3])
        self.ref_strand = 1 if header[4] == '+' else -1
        self.ref_start = int(header[5])
        self.ref_end = int(header[6])
        self.query_chrome = header[7]
        self.query_size = int(header[8])
        self.query_strand = 1 if header[9] == '+' else -1
        self.query_start = int(header[10])
        self.query_end = int(header[11])
        self.id = int(header[12])
        self.alignments = Alignment(tuple(tuple(int(i) for i in alignment.split('\t')) for alignment in raw[1:]))
        if self.ref_strand == -1:
            self.alignments.__reversed__()
            self.ref_start, self.ref_end = self.ref_size - self.ref_end, self.ref_size - self.ref_start
            self.query_start, self.query_end = self.query_size - self.query_end, self.query_size - self.query_start

    @property
    def range(self):
        return self.ref_chrome, (self.ref_start, self.ref_end)

    def __contains__(self, item):
        if self.range[0] == item[0]:
            if self.range[1][0] <= item[1] <= self.range[1][1]:
                return True
        return False

    def __len__(self):
        return len(self.alignments)

    def __getitem__(self, i):
        return self.alignments[i]

    def __gt__(self, other):
        return self.id > other.id

    def __eq__(self, other):
        return self == other

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.id < other.id

    def __str__(self):
        return f'chain: {self.id}'


class Alignment:
    def __init__(self, raw: Tuple[Tuple[int]]):
        all_ = [j for i in raw for j in i]
        self._align_len = all_[0::3]
        self._ref_gap = all_[1::3]
        self._query_gap = all_[2::3]

    def __reversed__(self):
        self._align_len.reverse()
        self._ref_gap.reverse()
        self._query_gap.reverse()

    def __len__(self):
        return len(self.align_len)

    @property
    def align_len(self):
        return self._align_len

    @property
    def ref_gap(self):
        return self._ref_gap

    @property
    def query_gap(self):
        return self._query_gap

    def __getitem__(self, i):
        if i + 1 != self.__len__():
            return self._align_len[i], self._ref_gap[i], self._query_gap[i]
        else:
            return self.align_len[-1], 0, 0


def read_chain_file(file: str) -> Tuple[Seq]:
    chains = get_chain(file)
    return tuple(sorted(Seq(chain) for chain in chains))


class Trans:
    def __init__(self, seqs: Tuple[Seq]):
        self.seqs = seqs
        self.seq_idx = 0
        self.alignment_idx = 0
        self.ref_chrome: str = ''
        self.ref_size: int = 0
        self.ref_strand: int = 0
        self.ref_start: int = 0
        self.ref_end: int = 0
        self.query_chrome: str = ''
        self.query_size: int = 0
        self.query_strand: int = 0
        self.query_start: int = 0
        self.query_end: int = 0
        self.ref_pos: int = 0
        self.query_pos: int = 0
        self.align_len: int = 0
        self.align_range: Tuple[int, int] = (0, 0)
        self.get_chain_info(0)

    def get_trans(self, chrome: str, pos: int):
        if (chrome, pos) in self.seqs[self.seq_idx]:
            return self.get_trans_alignment(pos)
        for seq_idx in self.get_seq_iter():
            self.get_chain_info(seq_idx)
            if (chrome, pos) in self.seqs[self.seq_idx]:
                return self.get_trans_alignment(pos)
        return 'Fail', 'Fail'

    def __call__(self, chrom: str, pos: int):
        return self.get_trans(chrom, pos)

    def get_trans_alignment(self, pos: int):
        for _ in range(self.align_len):
            if self.align_range[0] <= pos <= self.align_range[1]:
                offset = pos - self.ref_pos
                return self.query_chrome, self.query_pos + offset * self.query_strand
            else:
                self.next_alignment()
        return self.query_chrome, 'Fail'

    def get_chain_info(self, idx: int):
        seq = self.seqs[idx]
        self.seq_idx = idx
        self.ref_chrome = seq.ref_chrome
        self.query_chrome = seq.query_chrome
        self.ref_start = seq.ref_start
        self.ref_strand = seq.ref_strand
        self.ref_size = seq.ref_size
        self.ref_end = seq.ref_end
        self.query_start = seq.query_start
        self.query_strand = seq.query_strand
        self.query_size = seq.query_size
        self.query_end = seq.query_end
        self.align_len = len(self.seqs[self.seq_idx].alignments)
        self.init_chain()

    def next_alignment(self):
        if self.alignment_idx + 1 < self.align_len:
            alignment = self.seqs[self.seq_idx].alignments[self.alignment_idx]
            self.ref_pos = self.ref_pos + alignment[0] + alignment[1]
            self.query_pos = self.query_pos + (alignment[0] + alignment[2]) * self.query_strand
            self.alignment_idx += 1
            self.align_range = self.ref_pos, self.ref_pos + self.seqs[self.seq_idx].alignments[self.alignment_idx][0]
        else:
            self.init_chain()

    def init_chain(self):
        self.ref_pos = self.ref_start
        self.query_pos = self.query_start
        self.alignment_idx = 0
        self.align_range = self.ref_start, self.ref_start + self.seqs[self.seq_idx].alignments[0][0]

    def get_seq_iter(self):
        seq_idxs = list(range(len(self.seqs)))
        seq_idxs.remove(self.seq_idx)
        return seq_idxs


def liftover_part(raw_file: str, start: int, end: int, idx: int, convertor: Trans, temp_dir: str,
                  columns_idx: List[int], sep_str: str, new_columns: bool, chr_prefix: str) -> None:
    with open(os.path.join(temp_dir, str(idx)), 'w', encoding='UTF-8') as f_write_part:
        for _, line in read_file_part(raw_file, start, end):
            line = line.decode().strip().split(sep_str)
            res = list(convertor(format_chrome(line[columns_idx[0]], chr_prefix), int(line[columns_idx[1]])))
            res[0] = reformat_chrom(res[0], chr_prefix)
            res[1] = str(res[1])
            if new_columns:
                line += res
            else:
                line[columns_idx[0]] = res[0]
                line[columns_idx[1]] = res[1]
            f_write_part.write(sep_str.join(line) + '\n')


def liftover(raw_file: str, chain_file: str, output_file: str, processes: Optional[int] = None, chr_prefix: str = '') \
        -> None:
    if not processes:
        processes = cpu_count()
    chains = read_chain_file(chain_file)
    convertor = Trans(chains)
    temp_dir = mkdtemp(suffix='liftover')
    time0 = time()
    print("Start...")
    try:
        with open_(raw_file) as f_raw, \
                open_(output_file, 'wb') as f_write:
            header = f_raw.readline()
            header = header.decode().strip().split(sep_str_l)
            columns_idx = get_columns_index(header, list(columns_name_l.values()))
            assert None not in columns_idx,\
                f"Cannot find column " \
                f"{list(compress(['CHR', 'POS'], map(lambda a: a is None, columns_idx)))} in file!"
            if new_columns:
                f_write.write(sep_str_l.join(header +
                                             list(map(lambda a: a+'_new', select_list(header, columns_idx)))).encode() +
                              b'\n')
            else:
                f_write.write(sep_str_l.join(header).encode() + b'\n')
        with Pool(processes=processes) as pool:
            idx = 0
            for start, end in file_div(raw_file, processes, 1, min_step=100):
                idx += 1
                pool.apply_async(liftover_part,
                                 args=(raw_file, start, end, idx, convertor,
                                       temp_dir, columns_idx, sep_str_l, new_columns, chr_prefix))
            pool.close()
            pool.join()
        with open_(output_file, 'ab') as f_write:
            for i in range(1, idx + 1):
                with open(os.path.join(temp_dir, str(i)), 'rb') as f_write_part:
                    for line in f_write_part:
                        f_write.write(line)
    finally:
        rm_dir(temp_dir)
    print(f'Finish! Used time: {time() - time0:.2f}s')


def liftover_run(raw_file: str, chain_file: str, output_file: str, processes: Optional[int] = None,
        sep: str = '\t', chr: str = 'CHR', pos: str = 'BP', new_column: bool = True, chr_prefix: str = '') -> None:
    global columns_name_l, sep_str_l, new_columns
    columns_name_l['CHR'].append(chr)
    columns_name_l['BP'].append(pos)
    sep_str_l = sep
    new_columns = new_column
    liftover(raw_file, chain_file, output_file, processes, chr_prefix)


if __name__ == '__main__':
    raw_file = '../add_rsid/test.tsv'
    output_file = '../add_rsid/test.tsv.gz'
    test_file = r'E:\Data\Database\LiftOver\hg19ToHg38.over.chain.gz'
    liftover_run(raw_file, test_file, output_file, processes=9, sep=' ', new_column=True, chr_prefix='CHR')
