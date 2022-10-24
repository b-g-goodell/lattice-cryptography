from math import ceil, floor, log2
from copy import deepcopy
from typing import Any
from secrets import randbits

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


def bit_rev(x: int, length: int) -> int:
    if isinstance(length, int) and length >= 1 and isinstance(x, int) and 0 <= x < 2**length:
        return int(bin(x)[2:].zfill(length)[::-1], 2)
    raise ValueError(f'Cannot compute bit_rev for x, length unless length is an integer with length >= 1 and x is an length-bit integer, but had (x,length)={(x, length)}.')


def is_pow_two(x: int) -> bool:
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    return False


def bit_rev_cp(x: list[int], num_bits: int) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(num_bits, int) and num_bits >= 1 and len(x) == 2**num_bits:
        return [x[bit_rev(x=i, length=num_bits)] for i in range(len(x))]
    raise ValueError(f'Cannot compute bit_reverse_cp for x, num_bits unless x is a list of integers, num_bits is an integer with num_bits >= 1, and the len(x)==2**num_bits, but had (x,num_bits)={(x, num_bits)}.')


def cent_rem(x: int) -> int:
    if isinstance(x, int):
        y: int = x % MODULUS
        z: int = y - HALF_MODULUS - 1
        return y - (1 + (z >> LOG_MODULUS)) * MODULUS
    raise ValueError(f'Cannot compute cent_rem for x unless x is an integer, but had x={x}.')


def ntt(x: list[int], inv_flag: bool, const_time: bool = True) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and is_pow_two(len(x)):
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
    raise ValueError(f'Cannot compute ntt for x, inv_flag, const_time unless x is a list of integers, inv_flag and const_time are both boolean, and x has power-of-two-length, but had (x, inv_flag, const_time)={(x, inv_flag, const_time)}.')


def encode(x: list[int], num_bits: int) -> bytes:
    if isinstance(num_bits, int) and num_bits >= 1 and isinstance(x, list) and all(isinstance(y, int) for y in x) and len(x) == DEGREE and all(0 <= y < 2 ** num_bits for y in x):
        return sum(bin(z)[2:].zfill(num_bits).encode() for z in x)
    raise ValueError(f'Cannot compute encode for x, num_bits unless num_bits is an integer with num_bits >= 1, x is a list of integers, x has length {DEGREE}, and each entry in x is an num_bits-bit integer, but had (x,num_bits)={(x, num_bits)}.')


def decode(x: bytes, num_32_bytes: int):
    if isinstance(x, bytes) and isinstance(num_32_bytes, int) and num_32_bytes >= 1 and len(x) >= 32*num_32_bytes:
        y: str = x.decode()
        if len(y) == DEGREE * num_32_bytes:
            return [int(y[i * num_32_bytes: (i + 1) * num_32_bytes], 2) % 2 ** num_32_bytes for i in range(DEGREE)]
    raise ValueError(f'Cannot compute decode for x, num_32_bytes unless x is a bytes object, num_32_bytes is an integer with num_32_bytes >= 1, and x has at least 32*num_32_bytes bytes, but had (type(x), len(x), num_32_bytes)={(type(x), len(x), num_32_bytes)}.')


def compress_int(x: int, num_bits: int) -> int:
    if isinstance(x, int) and isinstance(num_bits, int) and num_bits >= 1:
        y: float = (x % MODULUS) * (2 ** num_bits / MODULUS)
        ceil_y: int = ceil(y)
        if ceil_y - y <= 0.5:
            return ceil_y % 2 ** num_bits
        return floor(y) % 2 ** num_bits
    raise ValueError(f'Cannot compute compress_int for x, num_bits unless x and num_bits are both integers with num_bits >= 1, but had (x, num_bits)={(x, num_bits)}.')


def decompress_int(x: int, num_bits: int) -> int:
    if isinstance(num_bits, int) and num_bits >= 1 and isinstance(x, int) and 0 <= x < 2**num_bits:
        y: float = x * (MODULUS / 2 ** num_bits)
        ceil_y: int = ceil(y)
        if ceil_y - y <= 0.5:
            return ceil_y
        return floor(y)
    raise ValueError(f'Cannot compute decompress_int for x, num_bits unless x and num_bits are both integers with num_bits >= 1 and x is a num_bits-bit integer, but had (x, num_bits)={(x, num_bits)}.')


class PolynomialMatrix(object):
    modulus: int
    degree: int
    halfmod: int
    logmod: int
    vals: list[list[list[int]]]
    const_time: bool = True

    def __init__(self, vals: list[list[list[int]]], modulus: int, degree: int):
        if isinstance(modulus, int) and modulus >= 1 and isinstance(degree, int) and degree >= 1 and isinstance(vals, list) and all(isinstance(y, list) for y in vals) and all(isinstance(y, list) for i in vals for y in i) and all(isinstance(y, int) for i in vals for u in i for v in u for y in v):
            self.modulus = modulus
            self.degree = degree
            self.halfmod = modulus//2
            self.logmod = ceil(log2(modulus))
            self.vals = vals
        raise ValueError(f'Cannot initialize a PolynomialMatrix object with vald, modulus, degree unless modulus and degree are both integers with modulus >= 1 and degree >= 1 and vals is a list of lists of lists of integers, but had (vals, modulus, degree)={(vals, modulus, degree)}.')


class PolynomialCoefficientMatrix(PolynomialMatrix):
    def __init__(self, vals: list[list[list[int]]] | None, modulus: int, degree: int):
        super().__init__(vals=vals, modulus=modulus, degree=degree)

    def ntt(self) -> PolynomialMatrix:
        # Note: inv_flag input to ntt is False
        ntt_reps = []
        for row in self.vals:
            ntt_reps += [[]]
            for col in row:
                ntt_reps[-1] += [ntt(x=col, inv_flag=False, const_time=self.const_time)]
        return PolynomialNTTMatrix(vals=ntt_reps, modulus=self.modulus, degree=self.degree)

    def __add__(self, other):
        raise ValueError('Addition with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __radd__(self, other):
        raise ValueError('Addition with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __sub__(self, other):
        raise ValueError('Subtraction with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def __mul__(self, other):
        raise ValueError('Multiplication with PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def _scalar_mul(self, other):
        raise ValueError('Scalar multiplication of PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def _matrix_mul(self, other):
        raise ValueError('Matrix multiplication of PolynomialCoefficientMatrix is not implemented. Perform arithmetic only with PolynomialNTTMatrix.')

    def encode(self, num_bits: int) -> bytes:
        if isinstance(num_bits, int) and num_bits >= 1:
            return sum(encode(x=y, num_bits=num_bits) for i in self.vals for y in i)
        raise ValueError(f'Cannot compute PolynomialCoefficientMatrix.encode with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')

    def compress(self, num_bits: int):
        if isinstance(num_bits, int) and num_bits >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = compress_int(x=coef, num_bits=num_bits)
        raise ValueError(f'Cannot compute PolynomialCoefficientMatrix.compress with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')

    def decompress(self, num_bits: int):
        if isinstance(num_bits, int) and num_bits >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = decompress_int(x=coef, num_bits=num_bits)
        raise ValueError(f'Cannot compute PolynomialCoefficientMatrix.decompress with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')


class PolynomialNTTMatrix(PolynomialMatrix):
    def __init__(self, vals: list[list[list[int]]] | None, modulus: int, degree: int):
        super().__init__(vals=vals, modulus=modulus, degree=degree)

    def inv_ntt(self) -> PolynomialMatrix:
        # Note: inv_flag input to ntt is True.
        result = []
        for row in self.vals:
            result += [[]]
            for col in row:
                result += [ntt(x=col, inv_flag=True, const_time=self.const_time)]
        return PolynomialCoefficientMatrix(vals=result, modulus=self.modulus, degree=self.degree)

    def __add__(self, other):
        num_rows_in_self: int = len(self.vals)
        min_cols_in_self: int = min(len(x) for x in self.vals)
        max_cols_in_self: int = max(len(x) for x in self.vals)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self

        min_deg_in_self: int = min(len(x) for i in self.vals for x in i)
        max_deg_in_self: int = max(len(x) for i in self.vals for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self

        num_rows_in_other: int = len(other.vals)
        min_cols_in_other: int = min(len(x) for x in other.vals)
        max_cols_in_other: int = max(len(x) for x in other.vals)
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        min_deg_in_other: int = min(len(x) for i in other.vals for x in i)
        max_deg_in_other: int = max(len(x) for i in other.vals for x in i)
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        same_rows: bool = num_rows_in_self == num_rows_in_other
        same_cols: bool = max_cols_in_self == max_cols_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other
        if same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            result = deepcopy(self)
            for i, row in enumerate(other.vals):
                for j, col in enumerate(row):
                    for k, x in enumerate(col):
                        result.vals[i][j][k] = cent_rem(x=result.vals[i][j][k] + x)
            return result
        elif not consistent_cols_in_self or not consistent_deg_in_self or not consistent_cols_in_other or not consistent_deg_in_other:
            raise ValueError(f'Cannot compute PolynomialNTTMatrix.__add__ unless there are a consistent number of columns and degree in both self and other.')
        raise ValueError(f'Cannot compute PolynomialNTTMatrix.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')

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
        min_cols_in_self: int = min(len(x) for x in self.vals)
        max_cols_in_self: int = max(len(x) for x in self.vals)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        min_deg_in_self: int = min(len(x) for i in self.vals for x in i)
        max_deg_in_self: int = max(len(x) for i in self.vals for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        num_rows_in_other: int = len(other.vals)
        min_cols_in_other: int = min(len(x) for x in other.vals)
        max_cols_in_other: int = max(len(x) for x in other.vals)
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        min_deg_in_other: int = min(len(x) for i in other.vals for x in i)
        max_deg_in_other: int = max(len(x) for i in other.vals for x in i)
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other
        if num_rows_in_self == 1 and consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == 1 and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._scalar_mul(other=other)
        elif consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == num_rows_in_other and consistent_deg_in_self and consistent_deg_in_other and same_deg:
            return self._matrix_mul(other=other)
        raise ValueError('Cannot compute PolynomialNTTMatrix.__mul__ with unless self is 1x1 and both have consistent degrees, or where self is mxn, other is nxp, and both have consistent degrees (dim mismatch).')

    def _scalar_mul(self, other):
        result = deepcopy(other)
        for i in range(len(other.vals)):
            for j in range(len(other.vals[0])):
                for k in range(len(self.vals[0][0])):
                    result.ntt_reps[i][j][k] = cent_rem(x=self.vals[0][0][k] * other.vals[i][j][k])
        return result

    def _matrix_mul(self, other):
        result = deepcopy(self)
        result.ntt_reps = [[[0 for k in range(len(self.vals[0][0]))] for j in range(len(other.vals[0]))] for i in range(len(self.vals))]
        for i in range(len(self.vals)):
            for j in range(len(other.vals[0])):
                for k in range(len(self.vals[0][0])):
                    result.ntt_reps[i][j][k] = cent_rem(x=sum(
                        self.vals[i][l][k] * other.ntt_reps[l][j][k] for l in range(len(self.vals[0]))))
        return result

    def encode(self, num_bits: int) -> bytes:
        if isinstance(num_bits, int) and num_bits >= 1:
            return sum(encode(x=y, num_bits=num_bits) for i in self.vals for y in i)
        raise ValueError(f'Cannot compute PolynomialNTTMatrix.encode with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')

    def compress(self, num_bits: int):
        if isinstance(num_bits, int) and num_bits >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = compress_int(x=coef, num_bits=num_bits)
        raise ValueError(f'Cannot compute PolynomialNTTMatrix.compress with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')

    def decompress(self, num_bits: int):
        if isinstance(num_bits, int) and num_bits >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = decompress_int(x=coef, num_bits=num_bits)
        raise ValueError(f'Cannot compute PolynomialNTTMatrix.decompress with num_bits unless num_bits is an integer with num_bits >= 1, but had num_bits={num_bits}.')


def decode_coefficient_matrix(x: bytes, num_32_bytes: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialCoefficientMatrix:
    if isinstance(x, bytes) and isinstance(num_32_bytes, int) and num_32_bytes >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1 and isinstance(modulus, int) and modulus >= 1 and len(x) >= num_32_bytes*num_rows*num_cols*32:
        y = deepcopy(x)
        resulting_coef_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_coef_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32 * num_32_bytes], y[32 * num_32_bytes:]
                w = decode(x=z, num_32_bytes=num_32_bytes)
                resulting_coef_reps[-1] += [w]
        return PolynomialCoefficientMatrix(vals=resulting_coef_reps, modulus=modulus, degree=DEGREE)
    raise ValueError(f'Cannot compute decode_coefficient_matrix with x, num_32_bytes, num_rows, num_cols, modulus unless x is bytes and num_32_bytes (and num_rows, num_cols, modulus, respectively) are integer with num_32_bytes >= 1 (and num_rows >= 1, num_cols >= 1, and modulus >= 1, respectively), and x consists of at least num_32_bytes*num_rows*num_cols*32 bytes, but had len(x)={len(x)}.')


def decode_ntt_matrix(x: bytes, num_32_bytes: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialNTTMatrix:
    if isinstance(x, bytes) and isinstance(num_32_bytes, int) and num_32_bytes >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1 and isinstance(modulus, int) and modulus >= 1 and len(x) >= num_32_bytes*num_rows*num_cols*32*2:
        # For NTT rep, we have double the degree
        y = deepcopy(x)
        resulting_ntt_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_ntt_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32 * num_32_bytes * 2], y[32 * num_32_bytes * 2:]
                w = decode(x=z, num_32_bytes=num_32_bytes)
                resulting_ntt_reps[-1] += [w]
        return PolynomialNTTMatrix(vals=resulting_ntt_reps, modulus=modulus, degree=DEGREE)
    raise ValueError(f'Cannot compute decode_ntt_matrix with x, num_32_bytes, num_rows, num_cols, modulus unless x is bytes and num_32_bytes (and num_rows, num_cols, modulus, respectively) are integers with num_32_bytes >= 1 (and num_rows >= 1, num_cols >= 1, and modulus >= 1, respectively), and x consists of at least num_32_bytes*num_rows*num_cols*32 bytes, but had len(x)={len(x)}.')


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


def parse(x: bytes) -> list[list[int]]:
    if isinstance(x, bytes):
        pass
    raise ValueError(f'Cannot compute parse with x unless x is bytes, but had type(x)={type(x)}.')


def cbd(x: bytes, eta: int) -> int:
    if isinstance(eta, int) and eta >= 1 and isinstance(x, bytes) and len(x) >= 64*eta:
        z = x.decode()
        y: str = z[:len(z)//2]
        w: str = z[len(z)//2:]
        u: int = sum(int(i == '1') for i in y)
        v: int = sum(int(i == '1') for i in w)
        return u - v
    raise ValueError(f'Cannot compute cbd with x, eta unless x is bytes, eta is an integer with eta >= 1, and len(x) >= 64*eta, but had (eta, len(x))={(eta, len(x))}.')


def cpa_pke_keygen(security_parameter: int) -> tuple[bytes, bytes]:
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS:
        d: bytes = bin(randbits(256))[2:].encode()
        rho_and_sigma: bytes = HashFunctionG(d)
        rho: bytes
        sigma: bytes
        rho, sigma = rho_and_sigma[:len(rho_and_sigma)//2], rho_and_sigma[len(rho_and_sigma)//2:]
        N: int = 0
        A_hat: list[list[list[int]]] = []
        for i in range(PARAMS[security_parameter]['k']):
            A_hat += [[]]
            for j in range(PARAMS[security_parameter]['k']):
                x: bytes = rho + bin(j)[2:].encode() + bin(i)[2:].encode()
                A_hat[-1] += parse(XOF(x))
        A_hat: PolynomialNTTMatrix = PolynomialNTTMatrix(vals=A_hat, modulus=MODULUS, degree=DEGREE)
        s: list[list[list[int]]] = [[] for i in range(PARAMS[security_parameter]['k'])]
        for i in range(PARAMS[security_parameter]['k']):
            next_x: bytes = sigma + bin(N)[2:].encode()
            s[i] = cbd(x=PRF(x=next_x), eta=PARAMS[security_parameter]['eta_one'])
            N += 1
        s_hat: PolynomialNTTMatrix = PolynomialNTTMatrix(vals=s, modulus=MODULUS, degree=DEGREE)
        for i, row in enumerate(s_hat.vals):
            for j, col in enumerate(row):
                for k, coef in enumerate(col):
                    s_hat.vals[i][j][k] = s_hat.vals[i][j][k] % MODULUS
        e: list[list[list[int]]] = [[] for i in range(PARAMS[security_parameter]['k'])]
        for i in range(PARAMS[security_parameter]['k']):
            next_x: bytes = sigma + bin(N)[2:].encode()
            e[i] = cbd(x=next_x, eta=PARAMS[security_parameter]['k'])
            N += 1
        e_hat: PolynomialNTTMatrix = PolynomialNTTMatrix(vals=e, modulus=MODULUS, degree=DEGREE)
        t_hat: PolynomialNTTMatrix = A_hat * s_hat + e_hat
        for i, row in enumerate(t_hat.vals):
            for j, col in enumerate(row):
                for k, coef in enumerate(col):
                    t_hat.vals[i][j][k] = t_hat.vals[i][j][k] % MODULUS
        pk: bytes = t_hat.encode(num_bits=12) + rho
        sk: bytes = s_hat.encode(num_bits=12)
        return pk, sk
    raise ValueError(
        f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
        f'{ALLOWABLE_SECURITY_PARAMETERS}.')


def cpa_pke_encrypt(security_parameter: int, pk: bytes, m: bytes, r: bytes) -> bytes:
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS:
        N: int = 0
        encoded_t_hat: bytes = pk[:ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']]
        encoded_rho: bytes = pk[ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']:]
        t_hat: PolynomialNTTMatrix = decode_ntt_matrix(x=encoded_t_hat, num_32_bytes=12, num_rows=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k'], num_cols=1, modulus=MODULUS)
        A_hat_transpose_vals: list[list[list[int]]] = []
        for i in range(ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']):
            A_hat_transpose_vals += [[]]
            for j in range(ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']):
                A_hat_transpose_vals[-1] = parse(XOF(encoded_rho, i, j))
        A_hat_transpose: PolynomialNTTMatrix = PolynomialNTTMatrix(vals=A_hat_transpose_vals, modulus=MODULUS, degree=DEGREE)
        r_vals: list[list[list[int]]] = []
        for i in range(ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']):
            r_vals += [[cbd(x=PRF(r, N), eta=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['eta_one'])]]
            N += 1
        r: PolynomialCoefficientMatrix = PolynomialCoefficientMatrix(vals=r_vals, modulus=MODULUS, degree=DEGREE)
        e_one_vals: list[list[list[int]]] = []
        for i in range(ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['k']):
            e_one_vals += [[cbd(x=PRF(r, N), eta=ETA_TWO)]]
            N += 1
        e_one: PolynomialCoefficientMatrix = PolynomialCoefficientMatrix(vals=e_one_vals, modulus=MODULUS, degree=DEGREE)
        e_two_vals: list[list[list[int]]] = []
        e_two_vals += [[cbd(x=PRF(r, N), eta=ETA_TWO)]]
        e_two: PolynomialCoefficientMatrix = PolynomialCoefficientMatrix(vals=e_two_vals, modulus=MODULUS, degree=DEGREE)

        r_hat: PolynomialNTTMatrix = r.ntt()

        tmp_one: PolynomialNTTMatrix = A_hat_transpose * r_hat
        u: PolynomialCoefficientMatrix = tmp.inv_ntt() + e_one

        tmp_two: PolynomialNTTMatrix = t_hat * r_hat
        tmp_three: PolynomialCoefficientMatrix = tmp_two.inv_ntt() + e_two
        decoded_msg: PolynomialCoefficientMatrix = decode(x=m, num_32_bytes=1)
        decompressed_decoded_msg: PolynomialCoefficientMatrix = decoded_msg.decompress(num_bits=1)
        v: PolynomialCoefficientMatrix = tmp_three + decompressed_decoded_msg

        compressed_u: PolynomialCoefficientMatrix = u.compress(num_bits=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['d_u'])
        encoded_compresssed_u: bytes = compressed_u.encode(num_bits=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['d_u'])

        compressed_v: PolynomialCoefficientMatrix = v.compress(num_bits=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['d_v'])
        encoded_compressed_v: bytes = compressed_v.encode(num_bits=ALLOWABLE_SECURITY_PARAMETERS[security_parameter]['d_v'])

        return encoded_compresssed_u + encoded_compressed_v


    raise ValueError(
        f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
        f'{ALLOWABLE_SECURITY_PARAMETERS}.')


def cpa_pke_decrypt(security_parameter: int, sk: bytes, c: bytes) -> tuple[bytes, bytes]:
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS and isinstance(sk, bytes) and len(sk) >= ENCODED_CPA_PKE_SK_LEN[security_parameter] and isinstance(c, bytes) and len(c) >= ENCODED_CPA_PKE_CIPHERTEXT_LEN[security_parameter]:
        pass
    elif security_parameter not in ALLOWABLE_SECURITY_PARAMETERS:
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
