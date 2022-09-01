from math import ceil, floor, log2
from copy import deepcopy
from typing import Any
from lattice_algebra import LatticeParameters, Polynomial, PolynomialVector

ALLOWABLE_SECURITY_PARAMETERS: list[int] = [512, 768, 1024]
PARAMS: dict[int, dict[str, Any]] = {
    512: {'n': 256, 'k': 2, 'q': 3329, 'eta_one': 3, 'eta_two': 2, 'd_u': 10, 'd_v': 4, 'log_delta': -139, 'lp': LatticeParameters(degree=256, length=2, modulus=3329)},
    768: {'n': 256, 'k': 3, 'q': 3329, 'eta_one': 2, 'eta_two': 2, 'd_u': 10, 'd_v': 4, 'log_delta': -164, 'lp': LatticeParameters(degree=256, length=3, modulus=3329)},
    1024: {'n': 256, 'k': 4, 'q': 3329, 'eta_one': 2, 'eta_two': 2, 'd_u': 11, 'd_v': 5, 'log_delta': -174, 'lp': LatticeParameters(degree=256, length=4, modulus=3329)}
}

ENCODED_CPA_PKE_SK_LEN: dict[int, int] = {
    security_parameter: 12 * PARAMS[security_parameter]['k'] * PARAMS[security_parameter]['n'] / 8 
    for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CPA_PKE_PK_LEN: dict[int, int] = {
    security_parameter: 12*PARAMS[security_parameter]['k']*PARAMS[security_parameter]['n']/8 + 32
    for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CPA_PKE_CIPHERTEXT_LEN: dict[int, int] = {
    security_parameter: (PARAMS[security_parameter]['d_u']*PARAMS[security_parameter]['k'] +
                         PARAMS[security_parameter]['d_v'])*PARAMS[security_parameter]['n']/8
    for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CCA_KEM_SK_LEN: dict[int, int] = {
    security_parameter: 24*PARAMS[security_parameter]['k']*PARAMS[security_parameter]['n']/8 + 96
    for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CCA_KEM_PK_LEN: dict[int, int] = deepcopy(ENCODED_CPA_PKE_PK_LEN)
ENCODED_CCA_KEM_CIPHERTEXT_LEN: dict[int, int] = {
    security_parameter: (PARAMS[security_parameter]['d_u']*PARAMS[security_parameter]['k'] +
                         PARAMS[security_parameter]['d_v'])*PARAMS[security_parameter]['n']/8
    for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}


def XOF(x: bytes) -> bytes:
    pass


def PRF(x: bytes) -> bytes:
    pass


def compress_int(q: int, d: int, x: int) -> int:
    tmp: float = 2**d * x / q
    result: int = floor(tmp) % 2**d
    if ceil(tmp) - tmp <= 0.5:
        result: int = ceil(tmp) % 2**d
    return result


def compress_polynomial(q: int, d: int, x: Polynomial) -> Polynomial:
    pass


def compress_vector(q: int, d: int, x: PolynomialVector) -> PolynomialVector:
    pass


def compress_matrix(q: int, d: int, x: list[PolynomialVector]) -> list[PolynomialVector]:
    pass


def decompress_int(q: int, d: int, x: int) -> int:
    tmp: float = q * x / 2**d
    result: int = floor(tmp)
    if ceil(tmp) - tmp <= 0.5:
        result: int = ceil(tmp)
    return result


def decompress_polynomial(q: int, d: int, x: Polynomial) -> Polynomial:
    pass


def decompress_vector(q: int, d: int, x: PolynomialVector) -> PolynomialVector:
    pass


def decompress_matrix(q: int, d: int, x: list[PolynomialVector]) -> list[PolynomialVector]:
    pass


def parse(x: bytes) -> PolynomialVector:
    pass


def cbd_eta(eta: int, x: bytes) -> Polynomial:
    pass


def encode_polynomial(j: int, x: Polynomial) -> bytes:
    pass


def encode_vector(j: int, x: PolynomialVector) -> bytes:
    pass


def encode_matrix(j: int, x: list[PolynomialVector]) -> bytes:
    pass


def decode_polynomial(j: int, x: bytes) -> Polynomial:
    pass


def decode_vector(j: int, x: bytes) -> PolynomialVector:
    pass


def decode_matrix(j: int, x: bytes) -> list[PolynomialVector]:
    pass


def cpa_pke_keygen(security_parameter: int) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(
            f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
            f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    pass


def cpa_pke_encrypt(security_parameter: int, pk: bytes, m: bytes, seed_r: bytes) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(
            f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
            f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    elif not isinstance(pk, bytes):
        raise ValueError(f'Input public key must be bytes but had type(pk)={type(pk)}.')
    elif len(pk) < ENCODED_CPA_PKE_PK_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CPA-PKE{security_parameter} public keys need to be at least ENCODED_CPA_PKE_PK_LEN=' +
            f'{ENCODED_CPA_PKE_PK_LEN[security_parameter]} bytes long, but had len(pk)={len(pk)}.')
    elif not isinstance(m, str):
        raise ValueError(f'Input message must be bytes but had type(m)={type(m)}.')
    elif len(m) < 32:
        raise ValueError(
            f'Encoded Kyber-CPA-PKE{security_parameter} messages must be at least 32 bytes, but had len(m)={len(m)}.')
    elif not isinstance(seed_r, bytes):
        raise ValueError(f'Input random coins must be a string but had type(r)={type(seed_r)}.')
    elif len(seed_r) < 32:
        raise ValueError(
            f'Encoded Kyber-CPA-PKE{security_parameter} random coins need to be at least 32 bytes long, but had' +
            f'len(r)={len(seed_r)}.')
    N: int = 0
    decoded_pk: PolynomialVector = decode_vector(j=12, x=pk[:-32])
    rho: bytes = pk[-32:]
    A: list[PolynomialVector] = []
    pass


def cpa_pke_decrypt(security_parameter: int, sk: bytes, c: bytes) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(
            f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS' +
            f'={ALLOWABLE_SECURITY_PARAMETERS}.')
    elif not isinstance(sk, bytes):
        raise ValueError(f'Input sk must be bytes but had type(sk)={type(sk)}.')
    elif len(sk) < ENCODED_CPA_PKE_SK_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CPA-PKE{security_parameter} private keys need to be at least ENCODED_CPA_PKE_SK_LEN=' +
            f'{ENCODED_CPA_PKE_SK_LEN[security_parameter]} bytes long but had len(sk)={len(sk)}.')
    elif not isinstance(c, bytes):
        raise ValueError(f'Input ciphertext must be bytes, but had type(c)={type(c)}.')
    elif len(c) < ENCODED_CPA_PKE_CIPHERTEXT_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CPA-PKE{security_parameter} ciphertexts need to be at least ' +
            f'ENCODED_CPA_PKE_CIPHERTEXT_LEN={ENCODED_CPA_PKE_CIPHERTEXT_LEN[security_parameter]} bytes long ' +
            f'but had len(c)={len(c)}.')
    pass


def cca_kem_keygen(security_parameter: int) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(f'Must have security_parameter={security_parameter} in ALLOWABLE_security_parameters=' +
                         f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    pass


def cca_kem_encapsulate(security_parameter: int, pk: bytes) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(f'Must have security_parameter={security_parameter} in ALLOWABLE_security_parameters=' +
                         f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    elif not isinstance(pk, bytes):
        raise ValueError(f'Input public key must be bytes but had type(pk)={type(pk)}.')
    elif len(pk) < ENCODED_CCA_KEM_PK_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CCA-KEM{security_parameter} public keys need to be at least ENCODED_CCA_KEM_PK_LEN=' +
            f'{ENCODED_CCA_KEM_PK_LEN[security_parameter]} bytes long, but had len(pk)={len(pk)}.')
    pass


def cca_kem_decapsulate(security_parameter: int, sk: bytes, c: bytes) -> bytes:
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(f'Must have security_parameter={security_parameter} in ALLOWABLE_security_parameters=' +
                         f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    elif not isinstance(sk, bytes):
        raise ValueError(f'Input secret key must be bytes but had type(pk)={type(sk)}.')
    elif len(sk) < ENCODED_CCA_KEM_SK_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CCA-KEM{security_parameter} secret keys need to be at least ENCODED_CCA_KEM_SK_LEN=' +
            f'{ENCODED_CCA_KEM_SK_LEN[security_parameter]} bytes long, but had len(sk)={len(sk)}.')
    elif not isinstance(c, str):
        raise ValueError(f'Input ciphertext must be bytes but had type(c)={type(c)}.')
    elif len(c) < ENCODED_CCA_KEM_CIPHERTEXT_LEN[security_parameter]:
        raise ValueError(
            f'Encoded Kyber-CCA-KEM{security_parameter} ciphertexts need to be at least ' +
            f'ENCODED_CCA_KEM_CIPHERTEXT_LEN={ENCODED_CCA_KEM_CIPHERTEXT_LEN[security_parameter]} bytes long ' +
            f'but had len(c)={len(c)}.')
    pass
