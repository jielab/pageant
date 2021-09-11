import os
import rsa
import json
import pyDes
import qrcode
from re import compile
from zlib import decompress, compress
from pyzbar.pyzbar import decode, ZBarSymbol
from PIL import Image
from base64 import b85decode, b85encode
from typing import Tuple, List, Dict, Set, Optional


find_key = compile(r'(?<=-\n)[\S\n]*?(?=\n-)')
add_return = compile(r'\S{1,64}')


def trans_gt(gt: set or str, connector: str = '/') -> Set[str] or str:
    if type(gt) == set:
        return connector.join(list(gt)[0] * 2 if len(gt) == 1 else list(gt))  # Todo: not all biallelic
    elif type(gt) == str:
        for i in gt:
            if not i.isalpha() and i != '*':
                gt_n = gt.split(i)
                return set(gt_n)
        return set(gt)


def format_pub_key(key: str) -> bytes:
    return f'-----BEGIN RSA PUBLIC KEY-----\n{key}\n-----END RSA PUBLIC KEY-----\n'.encode()


def generate_key(directory: str, nbits: int = 1024, *args, **kwargs) -> None:
    public_key, private_key = rsa.newkeys(nbits, *args, **kwargs)
    public_key.save_pkcs1()
    if not os.path.isdir(directory):
        os.mkdir(directory)
    with open(os.path.join(directory, 'key'), 'wb') as f_key:
        f_key.write(private_key.save_pkcs1())
        f_key.write(public_key.save_pkcs1())


def load_key(key_file: str) -> Tuple[rsa.key.PublicKey, rsa.key.PrivateKey]:
    with open(key_file, 'rb') as f_key:
        key = f_key.read()
    return rsa.PublicKey.load_pkcs1(key), rsa.PrivateKey.load_pkcs1(key)


def extract_pub_key(key_file: str) -> str:
    public_key = load_key(key_file)[0]
    return find_key.search(public_key.save_pkcs1().decode())[0].replace('\n', '')


def make_qr_code(data: str, output: str, logo: Optional[str] = None, scale: int = 1):
    img_qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_M, box_size=12)
    img_qr.add_data(data)
    img_qr.make()
    img = img_qr.make_image().convert('RGB')
    if logo:
        logo = Image.open(logo)
        if scale != 1:
            logo.resize((logo.size[0] // scale, logo.size[1] // scale))
        pos = ((img.size[0] - logo.size[0]) // 2, (img.size[1] - logo.size[1]) // 2)
        img.paste(logo, pos)
    img.save(output)


def generate_send_qr(snp_list_txt: str, public_key: str, output: str, logo: Optional[str] = None) -> None:
    if os.path.isfile(snp_list_txt):
        with open(snp_list_txt) as f_snp:
            qr_code_snps = f_snp.read().strip()
    else:
        if ' ' in snp_list_txt:
            qr_code_snps = '\n'.join([i.strip() for i in snp_list_txt.split(' ') if i])
        else:
            qr_code_snps = snp_list_txt.strip()
    z_str = compress(qr_code_snps.encode())
    des = pyDes.des(public_key[-8:], pyDes.ECB, public_key[-8:], padmode=pyDes.PAD_PKCS5)
    secret_str = b85encode(des.encrypt(z_str)).decode()
    make_qr_code(f'{secret_str}\n{public_key}', output, logo=logo)


def decode_send_qr(qr_code: str) -> Tuple[List[str], rsa.key.PublicKey]:
    res = []
    i = 0
    while not res:
        if i > 10:
            raise Exception('No QR code was found')
        img = Image.open(qr_code)
        res = decode(img, symbols=[ZBarSymbol(64)])
        i += 1
    if len(res) > 1:
        print('More than one QR code were found')
    secret_str, pub_key_str = res[0].data.decode().split('\n')
    pub_key_str = '\n'.join(add_return.findall(pub_key_str))
    des = pyDes.des(pub_key_str[-8:], pyDes.ECB, pub_key_str[-8:], padmode=pyDes.PAD_PKCS5)
    pub_key = rsa.key.PublicKey.load_pkcs1(format_pub_key(pub_key_str))
    need_snp = decompress(des.decrypt(b85decode(secret_str.encode()))).decode()
    return need_snp.split('\n'), pub_key


def generate_give_qr(human: object, snp_list: List[str], pub_key: rsa.key.PublicKey, output: str,
                     logo: str or None = None) -> None:
    allele_text = ''
    for snps in snp_list:
        if snps in human.gt_data:
            allele_text += trans_gt(human.gt_data[snps], connector='')
        else:
            allele_text += 'NN'
    z_str = compress(allele_text.encode())
    limit_len = rsa.pkcs1.common.byte_size(pub_key.n) - 11
    secret_str = []
    for i in range(0, len(z_str), limit_len):
        secret_str.append(b85encode(rsa.encrypt(z_str, pub_key)).decode())
    make_qr_code('\n'.join(secret_str), output, logo, scale=3)


def decode_give_qr(give_qr_code: str, send_qr_code: str, private_key: rsa.key.PrivateKey) -> Dict[str, str]:
    res = []
    i = 0
    need_snp = decode_send_qr(send_qr_code)[0]
    while not res:
        if i > 10:
            raise Exception('No QR code was found')
        img = Image.open(give_qr_code)
        res = decode(img, symbols=[ZBarSymbol(64)])
        i += 1
    if len(res) > 1:
        print('More than one QR code were found')
    z_str = []
    try:
        for i in res[0].data.decode().split('\n'):
            z_str.append(rsa.decrypt(b85decode(i.encode()), private_key))
    except rsa.DecryptionError:
        raise rsa.DecryptionError('Error private key')
    genotype_text = ''.join([decompress(i).decode() for i in z_str])
    return {snp: genotype_text[idx * 2: (idx + 1) * 2] for idx, snp in enumerate(need_snp)}


def request(key_file: str, need_snp_list: str, save_img: str, logo: Optional[str] = None) -> None:
    # if not os.path.isdir(output):
    #     os.mkdir(output)
    # save_img = os.path.join(output, 'DR_QR_code.png')
    public_key_str = extract_pub_key(key_file)
    generate_send_qr(need_snp_list, public_key_str, save_img, logo)


def give(qr_code: str, human: object, output: str, logo: Optional[str] = None) -> None:
    snp_list, pub_key = decode_send_qr(qr_code)
    generate_give_qr(human, snp_list, pub_key, output, logo)


def obtain(qr_code: str, send_qr_code: str, key_file: str, output: str) -> None:
    private_key = load_key(key_file)[1]
    res = decode_give_qr(qr_code, send_qr_code, private_key)
    if not os.path.isdir(output):
        os.mkdir(output)
    with open(os.path.join(output, 'User_genotype.json'), 'w') as f:
        json.dump(res, f)
