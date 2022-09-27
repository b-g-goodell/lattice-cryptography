from math import ceil, floor, log2
from copy import deepcopy
from typing import Any

# Order of variables: (x, y, z, w) <- (i, j, k, l) <- (u, v, r, s)
# Metasyntactic list: foo, bar, baz, qux, quux, corge, grault, garply, waldo, fred, plugh, xyzzy, thud

ALLOWABLE_SECURITY_PARAMETERS: list[int] = [512, 768, 1024]
DEGREE: int = 256
MODULUS: int = 3329
ETA_TWO: int = 2
ROU: int = 17
ROU_INV: int = 1175
PARAMS: dict[int, dict[str, Any]] = {
     512: {'k': 2, 'eta_one': 3, 'd_u': 10, 'd_v': 4, 'log_delta': -139},
     768: {'k': 3, 'eta_one': 2, 'd_u': 10, 'd_v': 4, 'log_delta': -164},
    1024: {'k': 4, 'eta_one': 2, 'd_u': 11, 'd_v': 5, 'log_delta': -174}
}
TWICE_DEGREE: int = 2*DEGREE
LOG_TWICE_DEGREE: int = ceil(log2(TWICE_DEGREE))
HALF_MODULUS: int = MODULUS//2
LOG_MODULUS: int = ceil(log2(MODULUS))
ZETAS: list[int] = [(ROU**k) % MODULUS for k in [TWICE_DEGREE//(2**(s+1)) for s in range(LOG_TWICE_DEGREE)]]
ZETA_INVS: list[int] = [(ROU_INV**k) % MODULUS for k in [TWICE_DEGREE//(2**(s+1)) for s in range(LOG_TWICE_DEGREE)]]
ENCODED_CPA_PKE_SK_LEN: dict[int, int] = {security_parameter: 12 * PARAMS[security_parameter]['k'] * DEGREE / 8 for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CPA_PKE_PK_LEN: dict[int, int] = {security_parameter: 12*PARAMS[security_parameter]['k'] * DEGREE/8 + 32 for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CPA_PKE_CIPHERTEXT_LEN: dict[int, int] = {security_parameter: (PARAMS[security_parameter]['d_u']*PARAMS[security_parameter]['k'] + PARAMS[security_parameter]['d_v'])*DEGREE/8 for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CCA_KEM_SK_LEN: dict[int, int] = {security_parameter: 24*PARAMS[security_parameter]['k'] * DEGREE/8 + 96 for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}
ENCODED_CCA_KEM_PK_LEN: dict[int, int] = deepcopy(ENCODED_CPA_PKE_PK_LEN)
ENCODED_CCA_KEM_CIPHERTEXT_LEN: dict[int, int] = {security_parameter: (PARAMS[security_parameter]['d_u']*PARAMS[security_parameter]['k'] + PARAMS[security_parameter]['d_v'])*DEGREE/8 for security_parameter in ALLOWABLE_SECURITY_PARAMETERS}


def bit_rev(num_bits: int, x: int) -> int:
    if isinstance(num_bits, int) and isinstance(x, int) and 0 <= x < 2**num_bits:
        return int(bin(x)[2:].zfill(num_bits)[::-1], 2)
    raise ValueError(f'Cannot reverse the bit string for x with num_bits={num_bits} unless 0 <= x < 2**num_bits, but x={x}.')


def is_pow_two(x: int) -> bool:
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    raise ValueError(f'Cannot test x as an integer power of 2 unless x > 0, but had x={x}.')


def bit_rev_cp(x: list[int], n: int) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(n, int) and ceil(log2(len(x))) == n and is_pow_two(len(x)):
        return [x[bit_rev(num_bits=n, x=i)] for i in range(len(x))]
    raise ValueError(f'Can only bit-reverse-copy arrays with integer power-of-two lengths but had len(x)={len(x)}.')


def cent_rem(x: int) -> int:
    if isinstance(x, int):
        y: int = x % MODULUS
        z: int = y - HALF_MODULUS - 1
        return y - (1 + (z >> LOG_MODULUS)) * MODULUS
    raise ValueError('Cannot compute cent_rem for a non-integer.')


def ntt(x: list[int], inv_flag: bool, const_time: bool = True) -> list[int]:
    if sum(int(y) for y in bin(len(x))[2:]) == 1:
        bit_rev_x: list[int] = bit_rev_cp(x=x, n=ceil(log2(len(x))))
        m: int = 1
        for s in range(1, LOG_TWICE_DEGREE + 1):
            m *= 2
            if inv_flag:
                this_zeta: int = ZETA_INVS[s - 1]
            else:
                this_zeta: int = ZETAS[s - 1]
            for k in range(0, TWICE_DEGREE, m):
                w: int = 1
                for j in range(m // 2):
                    t: int = w * bit_rev_x[k + j + m // 2]
                    u: int = bit_rev_x[k + j]
                    if const_time:
                        bit_rev_x[k + j]: int = cent_rem(x=u + t)
                        bit_rev_x[k + j + m // 2]: int = cent_rem(x=u - t)
                    else:
                        bit_rev_x[k + j]: int = (u + t) % MODULUS
                        if bit_rev_x[k + j] > HALF_MODULUS:
                            bit_rev_x[k + j] = bit_rev_x[k + j] - MODULUS
                        bit_rev_x[k + j + m // 2]: int = (u - t) % MODULUS
                        if bit_rev_x[k + j + m // 2] > HALF_MODULUS:
                            bit_rev_x[k + j + m // 2] = bit_rev_x[k + j + m // 2] - MODULUS
                    w *= this_zeta
        if inv_flag:
            n_inv: int = 1
            while (n_inv * TWICE_DEGREE) % MODULUS != 1:
                n_inv += 1
            if const_time:
                bit_rev_x: list[int] = [cent_rem(x=(n_inv * i)) for i in bit_rev_x]
            else:
                bit_rev_x: list[int] = [(n_inv * i) % MODULUS for i in bit_rev_x]
                bit_rev_x = [i if i <= HALF_MODULUS else i - MODULUS for i in bit_rev_x]
        return bit_rev_x
    raise ValueError("Can only NTT arrays with lengths that are powers of two.")


def encode(x: list[int], l: int) -> bytes:
    if len(x) == DEGREE and all(0 <= y < 2**l for y in x):
        return sum([bin(z % 2**l).zfill(l) for z in x]).encode()
    raise ValueError()


def decode(x: bytes, l: int):
    y: str = x.decode()
    return [int(y[i*l: (i+1)*l], 2) % 2**l for i in range(DEGREE)]


class PolynomialMatrix(object):
    modulus: int
    degree: int
    halfmod: int
    logmod: int
    coef_reps: list[list[list[int]]]
    ntt_reps: list[list[list[int]]]
    const_time: bool = True

    def __init__(self, modulus: int, degree: int, coef_reps: list[list[list[int]]]):
        self.modulus = modulus
        self.degree = degree
        self.halfmod = modulus//2
        self.logmod = ceil(log2(modulus))
        self.coef_reps = coef_reps
        self._coef_reps_to_ntt_reps()

    def _coef_reps_to_ntt_reps(self):
        # call this to set self.ntt_reps to the NTT representation of the coefficient representation
        self.ntt_reps = []
        for row in self.coef_reps:
            self.ntt_reps += [[]]
            for col in row:
                self.ntt_reps[-1] += [ntt(x=col, inv_flag=False, const_time=self.const_time)]

    def _ntt_reps_to_coef_reps(self):
        # call this to set self.coef_reps to the coefficient representation of the NTT representation
        self.coef_reps = []
        for row in self.ntt_reps:
            self.coef_reps += [[]]
            for col in row:
                self.coef_reps[-1] += [ntt(x=col, inv_flag=True, const_time=self.const_time)]

    def __add__(self, other):
        num_rows_in_self: int = len(self.ntt_reps)
        min_cols_in_self: int = min(len(x) for x in self.ntt_reps)
        max_cols_in_self: int = max(len(x) for x in self.ntt_reps)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        min_deg_in_self: int = min(len(x) for i in self.ntt_reps for x in i)
        max_deg_in_self: int = max(len(x) for i in self.ntt_reps for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        num_rows_in_other: int = len(other.ntt_reps)
        min_cols_in_other: int = min(len(x) for x in other.ntt_reps)
        max_cols_in_other: int = max(len(x) for x in other.ntt_reps)
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        min_deg_in_other: int = min(len(x) for i in other.ntt_reps for x in i)
        max_deg_in_other: int = max(len(x) for i in other.ntt_reps for x in i)
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        same_rows: bool = num_rows_in_self == num_rows_in_other
        same_cols: bool = max_cols_in_self == max_cols_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other
        if same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            result = deepcopy(self)
            for i, row in enumerate(other.ntt_reps):
                for j, col in enumerate(other.ntt_reps):
                    for k, x in enumerate(other.ntt_reps):
                        result.ntt_reps[i][j][k] = cent_rem(x=result.ntt_reps[i][j][k] + x)
            return result
        elif not consistent_cols_in_self or not consistent_deg_in_self or not consistent_cols_in_other or not consistent_deg_in_other:
            raise ValueError('Input matrices are not consistently sized.')
        else:
            raise ValueError('Dimension mismatch.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.ntt_reps = [[[-coef for coef in col] for col in row] for row in other.ntt_reps]
        return self.__add__(other=other)

    def __mul__(self, other):
        num_rows_in_self: int = len(self.ntt_reps)
        min_cols_in_self: int = min(len(x) for x in self.ntt_reps)
        max_cols_in_self: int = max(len(x) for x in self.ntt_reps)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        min_deg_in_self: int = min(len(x) for i in self.ntt_reps for x in i)
        max_deg_in_self: int = max(len(x) for i in self.ntt_reps for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        num_rows_in_other: int = len(other.ntt_reps)
        min_cols_in_other: int = min(len(x) for x in other.ntt_reps)
        max_cols_in_other: int = max(len(x) for x in other.ntt_reps)
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        min_deg_in_other: int = min(len(x) for i in other.ntt_reps for x in i)
        max_deg_in_other: int = max(len(x) for i in other.ntt_reps for x in i)
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other
        if num_rows_in_self == 1 and consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == 1 and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._scalar_mul(other=other)
        elif consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == num_rows_in_other and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._matrix_mul(other=other)
        else:
            raise ValueError('Can only multiply matrices A*B where A is 1x1 and both have consistent degrees, or where A is mxn, B is nxp, and both have consistent degrees.')

    def _scalar_mul(self, other):
        result = deepcopy(other)
        for i in range(len(other.ntt_reps)):
            for j in range(len(other.ntt_reps[0])):
                for k in range(len(self.ntt_reps[0][0])):
                    result.ntt_reps[i][j][k] = cent_rem(x=self.ntt_reps[0][0][k] * other.ntt_reps[i][j][k])
        return result

    def _matrix_mul(self, other):
        result = deepcopy(self)
        result.ntt_reps = [[[0 for k in range(len(self.ntt_reps[0][0]))] for j in range(len(other.ntt_reps[0]))] for
                           i in range(len(self.ntt_reps))]
        for i in range(len(self.ntt_reps)):
            for j in range(len(other.ntt_reps[0])):
                for k in range(len(self.ntt_reps[0][0])):
                    result.ntt_reps[i][j][k] = cent_rem(x=sum(
                        self.ntt_reps[i][l][k] * other.ntt_reps[l][j][k] for l in range(len(self.ntt_reps[0]))))
        return result

    def encode_coef_reps(self, l: int) -> bytes:
        return sum(encode(x=y, l=l) for i in self.coef_reps for y in i)

    def encode_ntt_reps(self, l: int) -> bytes:
        return sum(encode(x=y, l=l) for i in self.ntt_reps for y in i)


def decode_coef_reps(x: bytes, l: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialMatrix:
    if len(x) >= l*num_rows*num_cols*32:
        y = deepcopy(x)
        resulting_coef_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_coef_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32*l], y[32*l:]
                w = decode(x=z, l=l)
                resulting_coef_reps[-1] += [w]
        return PolynomialMatrix(modulus=modulus, degree=DEGREE, coef_reps=resulting_coef_reps)
    raise ValueError('Input does not have enough bytes to decode to a PolynomialMatrix.')


def decode_ntt_reps(x: bytes, l: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialMatrix:
    if len(x) >= l*num_rows*num_cols*32*2:
        # For NTT rep, we have double the degree
        y = deepcopy(x)
        resulting_ntt_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_ntt_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32*l*2], y[32*l*2:]
                w = decode(x=z, l=l)
                resulting_ntt_reps[-1] += [w]
        return PolynomialMatrix(modulus=modulus, degree=DEGREE, coef_reps=[[[0 for k in range(DEGREE)] for j in range(num_cols)] for i in range(num_rows)])
    raise ValueError('Input does not have enough bytes to decode to a PolynomialMatrix.')


def compress_int(d: int, x: int) -> int:
    y: float = x * (2**d / MODULUS)
    ceil_y: int = ceil(y)
    if ceil_y - y <= 0.5:
        return ceil_y % 2**d
    return floor(y) % 2**d


def decompress_int(d: int, x: int) -> int:
    y: float = x * (MODULUS/2**d)
    ceil_y: int = ceil(y)
    if ceil_y - y <= 0.5:
        return ceil_y
    return floor(y)


def XOF(x: bytes) -> bytes:
    # Replace with whatever your implementation of Kyber uses. NIST says SHAKE256 or AES.
    pass


def PRF(x: bytes) -> bytes:
    # Replace with whatever your implementation of Kyber uses. NIST says SHAKE256 or AES.
    pass


def parse(x: bytes):
    pass


def cbd(eta: int, x: bytes) -> int:
    z = x.decode()
    y: str = z[:len(z)//2]
    w: str = z[len(z)//2:]
    u: int = sum(int(i == '1') for i in y)
    v: int = sum(int(i == '1') for i in w)
    return u - v


def cpa_pke_keygen(security_parameter: int):
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS:
        pass
    raise ValueError(
        f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
        f'{ALLOWABLE_SECURITY_PARAMETERS}.')


def cpa_pke_encrypt(security_parameter: int, pk: bytes, m: bytes, seed_r: bytes):
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
    pass


def cpa_pke_decrypt(security_parameter: int, sk: bytes, c: bytes):
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


def cca_kem_keygen(security_parameter: int):
    if security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
        raise ValueError(f'Must have security_parameter={security_parameter} in ALLOWABLE_security_parameters=' +
                         f'{ALLOWABLE_SECURITY_PARAMETERS}.')
    pass


def cca_kem_encapsulate(security_parameter: int, pk: bytes):
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


def cca_kem_decapsulate(security_parameter: int, sk: bytes, c: bytes):
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
