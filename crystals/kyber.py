from math import ceil, floor, log2
from copy import deepcopy
from typing import Any
from secrets import randbits

# Order of variables: (x, y, z, w) <- (i, j, k, l) <- (u, v, r, s)
# Metasyntactic list: foo, bar, baz, qux, quux, corge, grault, garply, waldo, fred, plugh, xyzzy, thud

# Parameters from Kyber1024
N: int = 256
MODULUS: int = 3329
K: int = 4
ETA: int = 2
D_U: int = 11
D_V: int = 5
ROU: int = 17
ROU_INV: int = 1175
NEGATIVE_LOG_DELTA: int = 174  # in a sense, "target bits of security," although Kyber1024 claims only 128 bits of sec.
SEED_LEN: int = 32


# Convenient constants derived from the parameters above
TWICE_DEGREE: int = 2 * N
LOG_TWICE_DEGREE: int = ceil(log2(TWICE_DEGREE))
HALF_MODULUS: int = MODULUS // 2
LOG_MODULUS: int = ceil(log2(MODULUS))
LOG_K: int = ceil(log2(K))
ZETAS: list[int] = [(ROU ** i) % MODULUS for i in [TWICE_DEGREE // (2 ** (j + 1)) for j in range(LOG_TWICE_DEGREE)]]
ZETA_INVS: list[int] = [(ROU_INV ** i) % MODULUS for i in [TWICE_DEGREE // (2 ** (j + 1)) for j in range(LOG_TWICE_DEGREE)]]

CPA_PKE_SK_LEN: int = LOG_MODULUS * K * N // 8
CPA_PKE_PK_LEN: int = LOG_MODULUS * K * N // 8 + SEED_LEN
ENCODED_CPA_PKE_SK_LEN: int = CPA_PKE_SK_LEN // (K * N // 8)
ENCODED_CPA_PKE_PK_LEN: int = (CPA_PKE_PK_LEN - SEED_LEN) // (K * N // 8)
CPA_PKE_FIRST_CIPHERTEXT_LEN: int = D_U * K * N // 8
CPA_PKE_SECOND_CIPHERTEXT_LEN: int = D_V * N // 8
CPA_PKE_CIPHERTEXT_LEN: int = CPA_PKE_FIRST_CIPHERTEXT_LEN + CPA_PKE_SECOND_CIPHERTEXT_LEN

CCA_KEM_SK_LEN: int = 24 * K * N // 8 + 96
CCA_KEM_PK_LEN: int = deepcopy(CPA_PKE_PK_LEN)
CCA_KEM_CIPHERTEXT_LEN: int = D_U * K * N // 8 + D_V * N // 8


def int2bytes(x: int, length: int = LOG_K) -> bytes:
    return bin(x)[2:].encode().zfill(length)


def bytes2int(x: bytes) -> int:
    return int(x, 2)


# Utility functions for computing the NTT in constant time: bit_rev, is_pow_two, bit_rev_cp, reduce, ntt
def bit_rev(x: int, length: int) -> int:
    if isinstance(length, int) and length >= 1 and isinstance(x, int) and 0 <= x < 2 ** length:
        tmp: bytes = int2bytes(x=x, length=length)
        return bytes2int(x=tmp[::-1])  # reverse order
    raise ValueError(
        f'Cannot compute bit_rev for x, length unless length is an integer with length >= 1 and x is an length-bit integer, but had (x,length)={(x, length)}.')


def is_pow_two(x: int) -> bool:
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    return False


def bit_rev_cp(x: list[int], num_bits: int) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(num_bits,
                                                                                 int) and num_bits >= 1 and len(
            x) == 2 ** num_bits:
        return [x[bit_rev(x=i, length=num_bits)] for i in range(len(x))]
    raise ValueError(
        f'Cannot compute bit_reverse_cp for x, num_bits unless x is a list of integers, num_bits is an integer with num_bits >= 1, and the len(x)==2**num_bits, but had (x,num_bits)={(x, num_bits)}.')


def reduce(x: int) -> int:
    if isinstance(x, int):
        y: int = x % MODULUS
        z: int = y - HALF_MODULUS - 1
        return y - (1 + (z >> LOG_MODULUS)) * MODULUS
    raise ValueError(f'Cannot compute reduce for x unless x is an integer, but had x={x}.')


def round_without_bias(x: float) -> int:
    ceil_x: int = ceil(x)
    if ceil_x - x < 0.5:
        return ceil_x
    return floor(x)


def _parse_one(x: bytes) -> list[int]:
    # Our implementation does NOT parse according to specifications, because we only input a list of bytes, not a bytestream. See comments in the parse function.
    i: int = 0
    j: int = 0
    result: list[int] = []
    while j < N:
        d1 = x[i] + 256 * (x[i + 1] % 16)
        d2 = (x[i + 1] // 16) + 16 * x[i + 2]
        if d1 < MODULUS:
            result += [d1]
            j += 1
        if d2 < MODULUS and j < N:
            result += [d2]
            j += 1
        i += 3
    if len(result) < N:
        raise ValueError('Parsing failed!')
    return result


def _parse_many(x: bytes) -> list[list[list[int]]]:
    # Our implementation does NOT parse according to specifications, because we only input a list of bytes, not a bytestream. See comments in the parse function.
    result: list[list[list[int]]] = []
    for i in range(K):
        result += [[]]
        for j in range(K):
            next_bytes: bytes = x[(K*i+j)*((LOG_MODULUS + NEGATIVE_LOG_DELTA)//8): (K*i+j+1)*((LOG_MODULUS + NEGATIVE_LOG_DELTA)//8)]
            result[-1] += [[_parse_one(x=next_bytes)]]
    return result


def parse(x: bytes) -> list[int] | list[list[list[int]]]:
    # Our implementation does NOT parse according to specifications, because we only input a list of bytes, not a bytestream.
    # Our implementation is not compatible with the NIST standard: users will observe a different "A" matrix.
    # There is a small positive probability that parsing fails.
    if isinstance(x, bytes) and len(x) == N*(LOG_MODULUS + NEGATIVE_LOG_DELTA)//8:
        return _parse_one(x=x)
    elif isinstance(x, bytes) and len(x) >= K*K*N*(LOG_MODULUS + NEGATIVE_LOG_DELTA)//8:
        return _parse_many(x=x)
    raise ValueError(f'Cannot compute parse for x unless x is a bytes object of length {N*(LOG_MODULUS + NEGATIVE_LOG_DELTA)//8} or {K*K*N*(LOG_MODULUS + NEGATIVE_LOG_DELTA)//8} but had (type(x), len(x))={(type(x), len(x))}.')


def _cbd(x: bytes, eta: int) -> int:
    z = x.decode()
    y: str = z[:len(z) // 2]
    w: str = z[len(z) // 2:]
    u: int = sum(int(i == '1') for i in y)
    v: int = sum(int(i == '1') for i in w)
    return u - v


def cbd(x: bytes) -> int:
    if isinstance(x, bytes) and len(x) >= 64 * ETA:
        return _cbd(x=x, eta=ETA)
    raise ValueError(f'Cannot compute cbd for x unless x is a bytes object of length at least {64*ETA} but had len(x)={len(x)}.')


class PolyCoefs(object):
    modulus: int = MODULUS
    degree: int = N
    halfmod: int = HALF_MODULUS
    logmod: int = LOG_MODULUS
    vals: list[list[list[int]]] = []
    const_time: bool = True

    def __init__(self, vals: list[list[list[int]]] | None, modulus: int = MODULUS, degree: int = N):
        self.vals = vals
        self.modulus = modulus
        self.degree = degree

    def __add__(self, other):
        raise ValueError(
            'Addition with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __radd__(self, other):
        raise ValueError(
            'Addition with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __sub__(self, other):
        raise ValueError(
            'Subtraction with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __mul__(self, other):
        raise ValueError(
            'Multiplication with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def _scalar_mul(self, other):
        raise ValueError(
            'Scalar multiplication of PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def _matrix_mul(self, other):
        raise ValueError(
            'Matrix multiplication of PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')


class PolyNTT(object):
    modulus: int = MODULUS
    degree: int = N
    halfmod: int = HALF_MODULUS
    logmod: int = LOG_MODULUS
    vals: list[list[list[int]]] = []
    const_time: bool = True

    def __init__(self, vals: list[list[list[int]]] | None, modulus: int = MODULUS, degree: int = N):
        self.vals = vals
        self.modulus = modulus
        self.degree = degree

    def __add__(self, other):
        num_rows_in_self: int = len(self.vals)
        num_rows_in_other: int = len(other.vals)

        min_cols_in_self: int = min(len(x) for x in self.vals)
        min_cols_in_other: int = min(len(x) for x in other.vals)
        max_cols_in_self: int = max(len(x) for x in self.vals)
        max_cols_in_other: int = max(len(x) for x in other.vals)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other

        min_deg_in_self: int = min(len(x) for i in self.vals for x in i)
        min_deg_in_other: int = min(len(x) for i in other.vals for x in i)
        max_deg_in_self: int = max(len(x) for i in self.vals for x in i)
        max_deg_in_other: int = max(len(x) for i in other.vals for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other

        same_rows: bool = num_rows_in_self == num_rows_in_other
        same_cols: bool = max_cols_in_self == max_cols_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other

        if same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            result = deepcopy(self)
            for i, row in enumerate(other.vals):
                for j, col in enumerate(row):
                    for k, x in enumerate(col):
                        result.vals[i][j][k] = reduce(x=result.vals[i][j][k] + x)
            return result
        elif not consistent_cols_in_self or not consistent_deg_in_self or not consistent_cols_in_other or not consistent_deg_in_other:
            raise ValueError(
                f'Cannot compute PolynomialNTTMatrix.__add__ unless there are a consistent number of columns and degree in both self and other.')
        raise ValueError(
            f'Cannot compute PolynomialNTTMatrix.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.ntt_reps = [[[-coef for coef in col] for col in row] for row in other.vals]
        return self.__add__(other=other)

    def __mul__(self, other):
        num_rows_in_self: int = len(self.vals)
        num_rows_in_other: int = len(other.vals)

        min_cols_in_self: int = min(len(x) for x in self.vals)
        min_cols_in_other: int = min(len(x) for x in other.vals)
        max_cols_in_self: int = max(len(x) for x in self.vals)
        max_cols_in_other: int = max(len(x) for x in other.vals)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other

        min_deg_in_self: int = min(len(x) for i in self.vals for x in i)
        min_deg_in_other: int = min(len(x) for i in other.vals for x in i)
        max_deg_in_self: int = max(len(x) for i in self.vals for x in i)
        max_deg_in_other: int = max(len(x) for i in other.vals for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        
        # same_rows: bool = num_rows_in_self == num_rows_in_other
        # same_cols: bool = max_cols_in_self == max_cols_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other

        if num_rows_in_self == 1 and consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == 1 and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._scalar_mul(other=other)
        elif consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == num_rows_in_other and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._matrix_mul(other=other)
        raise ValueError(
            'Cannot compute PolynomialNTTMatrix.__mul__ with unless self is 1x1 and both have consistent degrees, or where self is mxn, other is nxp, and both have consistent degrees (dim mismatch).')

    def _scalar_mul(self, other):
        result = deepcopy(other)
        for i in range(len(other.vals)):
            for j in range(len(other.vals[0])):
                for k in range(len(self.vals[0][0])):
                    result.ntt_reps[i][j][k] = reduce(x=self.vals[0][0][k] * other.vals[i][j][k])
        return result

    def _matrix_mul(self, other):
        result = deepcopy(self)
        result.ntt_reps = [[[0 for k in range(len(self.vals[0][0]))] for j in range(len(other.vals[0]))] for i in
                           range(len(self.vals))]
        for i in range(len(self.vals)):
            for j in range(len(other.vals[0])):
                for k in range(len(self.vals[0][0])):
                    result.ntt_reps[i][j][k] = reduce(x=sum(
                        self.vals[i][l][k] * other.ntt_reps[l][j][k] for l in range(len(self.vals[0]))))
        return result


def _compress_one_int(x: int, d: int, q: int = MODULUS) -> int:
    return round_without_bias(x=x * 2 ** d / q) % 2 ** d


def _compress_list_of_ints(x: list[int], d: int, q: int = MODULUS) -> list[int]:
    return [_compress_one_int(x=y, d=d, q=q) for y in x]


def _compress_many_ints(x: list[list[list[int]]], d: int, q: int = MODULUS) -> list[list[list[int]]]:
    return [[_compress_list_of_ints(x=z, d=d, q=q) for z in y] for y in x]


def _compress_PolyCoefs(x: PolyCoefs, d: int, q: int = MODULUS) -> PolyCoefs:
    return PolyCoefs(vals=_compress_many_ints(x=x.vals, d=d, q=q), modulus=x.modulus, degree=x.degree)


def compress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, q: int = MODULUS) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
    if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _compress_one_int(x=x, d=d, q=q)
    elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _compress_list_of_ints(x=x, d=d, q=q)
    elif isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _compress_many_ints(x=x, d=d, q=q)
    elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _compress_PolyCoefs(x=x, d=d, q=q)
    raise ValueError(f'Cannot compute compress for x, d, q unless x is an integer, a list of integers, or a list of lists of lists of integers, or a  PolynomialCoefficientMatrix... and d is an integer >= 1 and q is an integer >= 2, but had (type(x), d, q)={(type(x), d, q)}.')


def _decompress_one_int(x: int, d: int, q: int = MODULUS) -> int:
    return round_without_bias(x=x * q / 2 ** d)


def _decompress_list_of_ints(x: list[int], d: int, q: int = MODULUS) -> list[int]:
    return [_decompress_one_int(x=y, d=d, q=q) for y in x]


def _decompress_many_ints(x: list[list[list[int]]], d: int, q: int = MODULUS) -> list[list[list[int]]]:
    return [[_decompress_list_of_ints(x=y, d=d, q=q) for z in y] for y in x]


def _decompress_PolyCoefs(x: PolyCoefs, d: int, q: int = MODULUS) -> PolyCoefs:
    return PolyCoefs(vals=_decompress_many_ints(x=x.vals, d=d, q=q), modulus=x.modulus, degree=x.degree)


def decompress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, q: int = MODULUS) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
    if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _decompress_one_int(x=x, d=d, q=q)
    elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _decompress_list_of_ints(x=x, d=d, q=q)
    elif isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _decompress_many_ints(x=x, d=d, q=q)
    elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(q, int) and q >= 2:
        return _decompress_PolyCoefs(x=x, d=d, q=q)
    raise ValueError(f'Cannot decompress with x, d, q unless x is an integer, or a list of integers, or a list of lists of lists of integers, or a PolynomialCoefficientMatrix... and d is an integer >= 1, and q is an integer >= 2, but had (type(x), d, q)={(type(x), d, q)}.')


def _encode_m_one_int(x: int, m: int) -> bytes:
    return int2bytes(x=x, length=m)  # bin(x)[2:].zfill(m).encode()


def _encode_m_list_of_ints(x: list[int], m: int) -> bytes:
    result = bytes(0)
    for y in x:
        result += _encode_m_one_int(x=y, m=m)
    return result


def _encode_m_many_ints(x: list[list[list[int]]], m: int) -> bytes:
    result = bytes(0)
    for y in x:
        for z in y:
            result += _encode_m_many_ints(x=z, m=m)
    return result


def _encode_m_matrix(x: PolyCoefs | PolyNTT, m: int) -> bytes:
    return _encode_m_many_ints(x=x.vals, m=m)


def encode_m(x: int | list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT, m: int) -> bytes:
    # We rename encode_l to encode_m to use m instead of l because l is an ambiguous character.
    if isinstance(m, int) and m >= 1 and isinstance(x, int) and 0 <= x < 2 ** m:
        return _encode_m_one_int(x=x, m=m)
    elif isinstance(m, int) and m >= 1 and isinstance(x, list) and all(isinstance(y, int) for y in x) and len(x) == N and all(0 <= y < 2 ** m for y in x):
        return _encode_m_list_of_ints(x=x, m=m)
    elif isinstance(m, int) and m >= 1 and isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and all(0 <= w < 2**m for y in x for z in y for w in z):
        return _encode_m_many_ints(x=x, m=m)
    elif isinstance(m, int) and m >= 1 and (isinstance(x, PolyCoefs) or isinstance(x, PolyNTT)):
        return _encode_m_matrix(x=x, m=m)
    raise ValueError(f'Cannot compute encode for (x, m) unless m >= 1 is an integer and x is an m-bit integer, or a list of m-bit integers, or x is a list of lists of lists of m-bit integers, or x is a PolynomialCoefficientMatrix, but had (type(x), m)={(type(x), m)}.')


def _decode_m_one_int(x: bytes, m: int) -> int:
    # light wrapper for bytes2int
    return bytes2int(x=x, length=m)


def _decode_m_list_of_ints(x: bytes, m: int) -> list[int]:
    return [_decode_m_one_int(x=y) for y in x]


def _decode_m_many(x: bytes, m: int) -> list[list[list[int]]]:
    result: list[list[list[int]]] = []
    y: list[bytes] = [x[i*N*m: (i+1)*N*m] for i in range(K)]
    for next_row in y:
        next_poly: list[int] = _decode_m_list_of_ints(x=next_row, m=m)
        result += [[next_poly]]
    return result


def _decode_m_matrix(x: bytes, m: int, ntt_matrix_flag: bool) -> PolyCoefs | PolyNTT:
    if ntt_matrix_flag:
        return PolyNTT(vals=_decode_m_many(x=x, m=m))
    return PolyCoefs(vals=_decode_m_many(x=x, m=m))


def decode_m(x: bytes, m: int, coef_matrix_flag: bool = False, ntt_matrix_flag: bool = False) -> list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT:
    if isinstance(x, bytes) and len(x) == N*m:
        return _decode_m_list_of_ints(x=x, m=m)
    elif isinstance(x, bytes) and len(x) >= K*N*m and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and not coef_matrix_flag and not ntt_matrix_flag:
        return _decode_m_many(x=x, m=m)
    elif isinstance(x, bytes) and len(x) >= K*N*m and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and (coef_matrix_flag ^ ntt_matrix_flag):
        return _decode_m_matrix(x=x, m=m, ntt_matrix_flag=ntt_matrix_flag)
    raise ValueError(f'Cannot compute decode_m for (x, m, coef_matrix_flag, ntt_matrix_flag) unless (i) x is a bytes object of length {N*m} or (ii) x is a bytes object with length at least {K*N*m} and both coef_matrix_flag and ntt_matrix_flag are booleans with at most one of them True, but had (type(x), len(x), type(coef_matrix_flag), coef_matrix_flag, type(ntt_matrix_flag), ntt_matrix_flag)={(type(x), len(x), type(coef_matrix_flag), coef_matrix_flag, type(ntt_matrix_flag), ntt_matrix_flag)}.')


def _ntt_one(x: list[int], inv_flag: bool, const_time: bool = True) -> list[int]:
    bit_rev_x: list[int] = bit_rev_cp(x=x, num_bits=ceil(log2(len(x))))
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
                    bit_rev_x[k + j]: int = reduce(x=u + t)
                    bit_rev_x[k + j + m // 2]: int = reduce(x=u - t)
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
            bit_rev_x: list[int] = [reduce(x=(n_inv * i)) for i in bit_rev_x]
        else:
            bit_rev_x: list[int] = [(n_inv * i) % MODULUS for i in bit_rev_x]
            bit_rev_x = [i if i <= HALF_MODULUS else i - MODULUS for i in bit_rev_x]
    return bit_rev_x


def _ntt_many(x: list[list[list[int]]], inv_flag: bool, const_time: bool = True) -> list[list[list[int]]]:
    return [[_ntt_one(x=z, inv_flag=inv_flag, const_time=const_time) for z in y] for y in x]


def _ntt_raw(x: list[int] | list[list[list[int]]], inv_flag: bool, const_time: bool = True) -> list[int] | list[list[list[int]]]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and is_pow_two(len(x)):
        return _ntt_one(x=x, inv_flag=inv_flag, const_time=const_time)
    elif isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and all(is_pow_two(len(z)) for y in x for z in y):
        return _ntt_many(x=x, inv_flag=inv_flag, const_time=const_time)
    raise ValueError(f'Cannot compute _ntt_raw for x, inv_flag, const_time unless x is a list of integers with power-of-two length or a list of lists of lists of integers (with the most internal of these lists having lengths that are powers-of-two), inv_flag and const_time are both bool, and x has power-of-two-length, but had (type(x), inv_flag, const_time)={(type(x), inv_flag, const_time)}.')


def ntt(x: PolyCoefs | PolyNTT) -> PolyCoefs | PolyNTT:
    if isinstance(x, PolyCoefs):
        return PolyNTT(vals=_ntt_raw(x.vals, inv_flag=False, const_time=True), modulus=x.modulus, degree=x.degree)
    elif isinstance(x, PolyNTT):
        return PolyCoefs(vals=_ntt_raw(x.vals, inv_flag=True, const_time=True), modulus=x.modulus, degree=x.degree)
    raise ValueError(f'Cannot compute NTT (or INTT) for x unless x is a PolynomialCoefficientMatrix or a PolynomialNTTMatrix, but had type(x)={x}.')


def transpose(x: PolyCoefs | PolyNTT) -> PolyCoefs | PolyNTT:
    if isinstance(x, PolyCoefs) or isinstance(x, PolyNTT):
        result = deepcopy(x)
        result.vals = [[x.vals[j][i] for j in range(len(x.vals))] for i in range(len(x.vals))]
        return result
    raise ValueError(f'Cannot transpose with x unless x is a PolynomialCoefficientMatrix or a PolynomialNTTMatrix, but had type(x)={type(x)}.')


def XOF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        # Pass input to whatever XOF your implementation requires.
        pass
    raise ValueError(f'Cannot compute XOF with x unless x is bytes, but had type(x)={type(x)}.')


def PRF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        # Pass input to whatever PRF your implementation requires.
        pass
    raise ValueError(f'Cannot compute PRF with x unless x is bytes, but had type(x)={type(x)}.')


def KDF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        # Pass input to whatever KDF your implementation requires.
        pass
    raise ValueError(f'Cannot compute KDF with x unless x is bytes, but had type(x)={type(x)}.')


def HashFunctionH(x: bytes) -> bytes:
    if isinstance(x, bytes):
        # Pass input to whatever hash function H your implementation requires.
        pass
    raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')


def HashFunctionG(x: bytes) -> bytes:
    if isinstance(x, bytes):
        # Pass input to whatever hash function G your implementation requires.
        pass
    raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')


def cpa_pke_keygen() -> bytes:
    d: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
    rho: bytes  # for typing
    sigma: bytes  # for typing
    rho_and_sigma: bytes = HashFunctionG(d)
    # assert len(rho_and_sigma) == 2 * SEED_LEN
    rho: bytes = rho_and_sigma[:len(rho_and_sigma) // 2]
    sigma: bytes = rho_and_sigma[len(rho_and_sigma) // 2:]
    index: int = 0

    # Make A_hat
    A_hat_vals: list[list[list[int]]] = []
    for i in range(K):
        A_hat_vals += [[]]
        for j in range(K):
            next_xof_input: bytes = rho + int2bytes(x=j) + int2bytes(x=i)  # first j, then i
            A_hat_vals[-1] += parse(XOF(next_xof_input))
    A_hat: PolyNTT = PolyNTT(vals=A_hat_vals)

    # Make s_hat
    s_vals: list[list[list[int]]] = [[] for _ in range(K)]
    for i in range(K):
        next_x: bytes = sigma + int2bytes(x=index)
        s_vals[i] = cbd(x=PRF(x=next_x))
        index += 1
    s_hat: PolyNTT = PolyNTT(vals=s_vals)

    # Make e_hat
    e_vals: list[list[list[int]]] = [[] for _ in range(K)]
    for i in range(K):
        next_x: bytes = sigma + int2bytes(x=index)
        e_vals[i] = cbd(x=next_x)
        index += 1
    e_hat: PolyNTT = PolyNTT(vals=e_vals)

    t_hat: PolyNTT = A_hat * s_hat + e_hat
    for i, row in enumerate(t_hat.vals):
        for j, col in enumerate(row):
            for k, coef in enumerate(col):
                t_hat.vals[i][j][k] = t_hat.vals[i][j][k] % MODULUS

    pk: bytes = encode_m(x=t_hat, m=ENCODED_CPA_PKE_PK_LEN) + rho
    sk: bytes = encode_m(x=s_hat, m=ENCODED_CPA_PKE_SK_LEN)
    return pk + sk


def _cpa_pke_enc(pk: bytes, plaintext: bytes, r: bytes) -> bytes:
    index: int = 0  # this is N in the specs, but we already used capital N for degree n
    encoded_t_hat: bytes = pk[:K]
    rho: bytes = pk[K:]
    t_hat: PolyNTT = decode_m(x=encoded_t_hat, m=ENCODED_CPA_PKE_PK_LEN)
    t_hat_transpose: PolyNTT = transpose(x=t_hat)

    A_hat_transpose_vals: list[list[list[int]]] = []
    for i in range(K):
        A_hat_transpose_vals += [[]]
        for j in range(K):
            next_xof_input: bytes = rho + int2bytes(x=i) + int2bytes(x=i)  # first i, then j
            A_hat_transpose_vals[-1] += [parse(XOF(next_xof_input))]
    A_hat_transpose: PolyNTT = PolyNTT(vals=A_hat_transpose_vals)

    r_vals: list[list[list[int]]] = []
    for i in range(K):
        next_prf_input: bytes = r + int2bytes(x=index)
        r_vals += [[cbd(x=PRF(next_prf_input))]]
        index += 1
    r: PolyCoefs = PolyCoefs(vals=r_vals)

    e_one_vals: list[list[list[int]]] = []
    for i in range(K):
        next_prf_input: bytes = r + int2bytes(x=index)
        e_one_vals += [[cbd(x=PRF(next_prf_input))]]
        index += 1
    e_one: PolyCoefs = PolyCoefs(vals=e_one_vals)

    e_two_vals: list[list[list[int]]] = []
    next_prf_input: bytes = r + int2bytes(x=index)
    e_two_vals += [[cbd(x=PRF(next_prf_input))]]
    e_two: PolyCoefs = PolyCoefs(vals=e_two_vals)

    r_hat: PolyNTT = ntt(x=r)
    partial_u_hat: PolyNTT = A_hat_transpose * r_hat
    u: PolyCoefs = ntt(x=partial_u_hat) + e_one
    partial_v_hat: PolyNTT = t_hat_transpose * r_hat
    partial_v: PolyCoefs = ntt(x=partial_v_hat)
    partial_v += e_two
    decoded_msg: PolyCoefs = decode_m(x=plaintext, m=1)
    decompressed_decoded_msg: PolyCoefs = decompress(x=decoded_msg, d=1)
    v: PolyCoefs = partial_v + decompressed_decoded_msg

    compressed_u: PolyCoefs = compress(x=u, d=D_U)
    encoded_compressed_u: bytes = encode_m(x=compressed_u, m=D_U)

    compressed_v: PolyCoefs = compress(x=v, d=D_V)
    encoded_compressed_v: bytes = encode_m(x=compressed_v, m=D_V)

    return encoded_compressed_u + encoded_compressed_v


def cpa_pke_encrypt(pk: bytes, plaintext: bytes, r: bytes) -> bytes:
    if isinstance(pk, bytes) and len(pk) >= CPA_PKE_PK_LEN and \
            isinstance(plaintext, bytes) and len(plaintext) >= SEED_LEN and \
            isinstance(r, bytes) and len(r) >= SEED_LEN:
        return _cpa_pke_enc(pk=pk, plaintext=plaintext, r=r)
    raise ValueError(
        f'Cannot cpa_pke_encrypt pk, plaintext, r unless pk is bytes with len(pk) >= {ENCODED_CPA_PKE_PK_LEN} ' +
        f'and plaintext is bytes with len(m) >= {SEED_LEN} and r is bytes len(r) >= {SEED_LEN}, but had ' +
        f'(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))=' +
        f'{(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))}.')


def _cpa_pke_dec(sk: bytes, ciphertext: bytes) -> bytes:
    encoded_c_one: bytes = ciphertext[:CPA_PKE_FIRST_CIPHERTEXT_LEN]
    encoded_c_two: bytes = ciphertext[CPA_PKE_FIRST_CIPHERTEXT_LEN:]
    c_one: PolyCoefs = decode_m(x=encoded_c_one, m=D_U)
    u: PolyCoefs = decompress(x=c_one, d=D_U)
    u_hat: PolyNTT = ntt(x=u)
    c_two: PolyCoefs = decode_m(x=encoded_c_two, m=D_V)
    v: PolyCoefs = decompress(x=c_two, d=D_V)
    s_hat: PolyNTT = decode_m(x=sk, m=ENCODED_CPA_PKE_SK_LEN)
    s_hat_dot_u_hat: PolyNTT = transpose(x=s_hat) * u_hat
    s_dot_u: PolyCoefs = ntt(x=s_hat_dot_u_hat)
    compressed_decoded_msg: PolyCoefs = v - s_dot_u
    decoded_msg: PolyCoefs = decompress(x=compressed_decoded_msg, d=1)
    msg: bytes = encode_m(x=decoded_msg, m=1)
    return msg


def cpa_pke_decrypt(sk: bytes, ciphertext: bytes) -> bytes:
    if isinstance(sk, bytes) and len(sk) >= CPA_PKE_SK_LEN and isinstance(ciphertext, bytes) and len(ciphertext) >= CPA_PKE_CIPHERTEXT_LEN:
        return _cpa_pke_dec(sk=sk, ciphertext=ciphertext)
    raise ValueError(f'Cannot cpa_pke_decrypt with sk, ciphertext unless sk is bytes with len(sk) >= ' +
                     f'{CPA_PKE_SK_LEN} and ciphertext is bytes with len(ciphertext) >= {CPA_PKE_CIPHERTEXT_LEN}' +
                     f' but had (type(sk), len(sk), type(ciphertext), len(ciphertext))=' +
                     f'{(type(sk), len(sk), type(ciphertext), len(ciphertext))}.')


def cca_kem_keygen() -> bytes:
    z: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
    cpa_pke_keys: bytes = cpa_pke_keygen()
    pk: bytes = cpa_pke_keys[:CPA_PKE_PK_LEN]
    sk_prime: bytes = cpa_pke_keys[CPA_PKE_PK_LEN:]
    hashed_pk: bytes = HashFunctionH(x=pk)
    sk: bytes = sk_prime + pk + hashed_pk + z
    return pk + sk


def _cca_kem_enc(pk: bytes) -> bytes:
    msg: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
    plaintext: bytes = HashFunctionH(x=msg)
    hashed_pk: bytes = HashFunctionH(x=pk)
    k_bar_and_r: bytes = HashFunctionG(plaintext + hashed_pk)
    k_bar: bytes = k_bar_and_r[:SEED_LEN]
    r: bytes = k_bar_and_r[SEED_LEN:]
    ciphertext: bytes = cpa_pke_encrypt(pk=pk, plaintext=plaintext, r=r)
    hashed_c: bytes = HashFunctionH(x=ciphertext)
    k: bytes = KDF(x=k_bar + hashed_c)
    return ciphertext + k


def cca_kem_encapsulate(pk: bytes) -> bytes:
    if isinstance(pk, bytes) and len(pk) >= CCA_KEM_PK_LEN:
        return _cca_kem_enc(pk=pk)
    raise ValueError(f'Cannot cca_kem_encapsulate for pk unless pk is bytes with length at least {CCA_KEM_PK_LEN} but had (type(pk), len(pk))={(type(pk), len(pk))}.')


def _cca_kem_dec(ciphertext: bytes, sk: bytes) -> bytes:
    pk: bytes = sk[CPA_PKE_SK_LEN:CPA_PKE_SK_LEN+CPA_PKE_PK_LEN]
    hashed_pk: bytes = sk[CPA_PKE_SK_LEN+CPA_PKE_PK_LEN:CPA_PKE_SK_LEN+CPA_PKE_PK_LEN+SEED_LEN]
    z: bytes = sk[CPA_PKE_SK_LEN+CPA_PKE_PK_LEN+SEED_LEN:]
    plaintext_prime: bytes = cpa_pke_decrypt(sk=sk, ciphertext=ciphertext)
    k_bar_prime_and_r_prime: bytes = HashFunctionG(x=plaintext_prime + hashed_pk)
    k_bar_prime: bytes = k_bar_prime_and_r_prime[:SEED_LEN]
    r_prime: bytes = k_bar_prime_and_r_prime[SEED_LEN:]
    ciphertext_prime: bytes = cpa_pke_encrypt(pk=pk, plaintext=plaintext_prime, r=r_prime)
    hashed_ciphertext: bytes = HashFunctionH(x=ciphertext)
    if ciphertext_prime == ciphertext:
        return KDF(x=k_bar_prime + hashed_ciphertext)
    return KDF(x=z + hashed_ciphertext)


def cca_kem_decapsulate(ciphertext: bytes, sk: bytes) -> bytes:
    if isinstance(ciphertext, bytes) and len(ciphertext) == CCA_KEM_CIPHERTEXT_LEN and isinstance(sk, bytes) and len(sk) == CCA_KEM_SK_LEN:
        return _cca_kem_dec(ciphertext=ciphertext, sk=sk)
    raise ValueError(f'Cannot cca_kem_decapsulate with ciphertext and sk unless ciphertext is bytes with length at least {CCA_KEM_CIPHERTEXT_LEN} and sk is bytes with length at least {CCA_KEM_SK_LEN}, but had (type(ciphertext), len(ciphertext), type(sk), len(sk))={(type(ciphertext), len(ciphertext), type(sk), len(sk))}.')