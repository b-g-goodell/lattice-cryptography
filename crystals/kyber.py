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


def bit_rev(l: int, x: int) -> int:
    if isinstance(l, int) and l >= 1 and isinstance(x, int) and 0 <= x < 2**l:
        return int(bin(x)[2:].zfill(l)[::-1], 2)
    raise ValueError(f'Cannot reverse the bit string for x with num_bits={l} unless 0 <= x < 2**num_bits, but x={x}.')


def is_pow_two(x: int) -> bool:
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    raise False


def bit_rev_cp(x: list[int], n: int) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(n, int) and ceil(log2(len(x))) == n and is_pow_two(len(x)):
        return [x[bit_rev(l=n, x=i)] for i in range(len(x))]
    raise ValueError(f'Can only bit-reverse-copy arrays with integer power-of-two lengths but had len(x)={len(x)}.')


def cent_rem(x: int) -> int:
    if isinstance(x, int):
        y: int = x % MODULUS
        z: int = y - HALF_MODULUS - 1
        return y - (1 + (z >> LOG_MODULUS)) * MODULUS
    raise ValueError('Cannot compute cent_rem for a non-integer.')


def ntt(x: list[int], inv_flag: bool, const_time: bool = True) -> list[int]:
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and is_pow_two(len(x)):
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
    if isinstance(l, int) and l >= 1 and isinstance(x, list) and all(isinstance(y, int) for y in x) and len(x) == DEGREE and all(0 <= y < 2**l for y in x):
        return sum(bin(z)[2:].zfill(l).encode() for z in x)
    raise ValueError()


def decode(x: bytes, l: int):
    if isinstance(x, bytes) and isinstance(l, int) and l >= 1 and len(x) >= 32*l:
        y: str = x.decode()
        if len(y) == DEGREE * l:
            return [int(y[i*l: (i+1)*l], 2) % 2**l for i in range(DEGREE)]
    raise ValueError()


def compress_int(d: int, x: int) -> int:
    if isinstance(d, int) and d >= 1 and isinstance(x, int):
        y: float = (x % MODULUS) * (2**d / MODULUS)
        ceil_y: int = ceil(y)
        if ceil_y - y <= 0.5:
            return ceil_y % 2**d
        return floor(y) % 2**d
    raise ValueError()


def decompress_int(d: int, x: int) -> int:
    if isinstance(d, int) and d >= 1 and isinstance(x, int) and 0 <= x < 2**d:
        y: float = x * (MODULUS/2**d)
        ceil_y: int = ceil(y)
        if ceil_y - y <= 0.5:
            return ceil_y
        return floor(y)
    raise ValueError()


class PolynomialMatrix(object):
    modulus: int
    degree: int
    halfmod: int
    logmod: int
    vals: list[list[list[int]]]
    const_time: bool = True

    def __init__(self, modulus: int, degree: int, vals: list[list[list[int]]] | None):
        if isinstance(modulus, int) and modulus >= 1 and isinstance(degree, int) and degree >= 1 and isinstance(vals, list) and all(isinstance(y, list) for y in vals) and all(isinstance(y, list) for i in vals for y in i) and all(isinstance(y, list) for i in vals for u in i for y in u) and all(isinstance(y, int) for i in vals for u in i for v in u for y in v) and vals is not None:
            self.modulus = modulus
            self.degree = degree
            self.halfmod = modulus//2
            self.logmod = ceil(log2(modulus))
            self.vals = vals
        raise ValueError()


class PolynomialCoefficientMatrix(PolynomialMatrix):
    def __init__(self, modulus: int, degree: int, vals: list[list[list[int]]] | None):
        super().__init__(modulus=modulus, degree=degree, vals=vals)

    def ntt(self) -> PolynomialMatrix:
        # call this to set self.ntt_reps to the NTT representation of the coefficient representation
        ntt_reps = []
        for row in self.vals:
            ntt_reps += [[]]
            for col in row:
                ntt_reps[-1] += [ntt(x=col, inv_flag=False, const_time=self.const_time)]
        return PolynomialNTTMatrix(modulus=self.modulus, degree=self.degree, vals=ntt_reps)

    def __add__(self, other):
        # num_rows_in_self: int = len(self.ntt_reps)
        # min_cols_in_self: int = min(len(x) for x in self.ntt_reps)
        # max_cols_in_self: int = max(len(x) for x in self.ntt_reps)
        # consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        # min_deg_in_self: int = min(len(x) for i in self.ntt_reps for x in i)
        # max_deg_in_self: int = max(len(x) for i in self.ntt_reps for x in i)
        # consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        # num_rows_in_other: int = len(other.ntt_reps)
        # min_cols_in_other: int = min(len(x) for x in other.ntt_reps)
        # max_cols_in_other: int = max(len(x) for x in other.ntt_reps)
        # consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        # min_deg_in_other: int = min(len(x) for i in other.ntt_reps for x in i)
        # max_deg_in_other: int = max(len(x) for i in other.ntt_reps for x in i)
        # consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        # same_rows: bool = num_rows_in_self == num_rows_in_other
        # same_cols: bool = max_cols_in_self == max_cols_in_other
        # same_deg: bool = max_deg_in_self == max_deg_in_other
        # if same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg:
        #     result = deepcopy(self)
        #     for i, row in enumerate(other.ntt_reps):
        #         for j, col in enumerate(other.ntt_reps):
        #             for k, x in enumerate(other.ntt_reps):
        #                 result.vals[i][j][k] = cent_rem(x=result.vals[i][j][k] + x)
        #     return result
        # elif not consistent_cols_in_self or not consistent_deg_in_self or not consistent_cols_in_other or not consistent_deg_in_other:
        #     raise ValueError('Input matrices are not consistently sized.')
        # else:
        #     raise ValueError('Dimension mismatch.')
        raise ValueError('Warning: You are trying to do arithmetic with a PolynomialCoefficientMatrix, which is much slower than arithmetic with a PolynomialNTTMatrix. Consider creating the associated PolynomialNTTMatrix with self.ntt() before performing arithmetic.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.ntt_reps = [[[-coef for coef in col] for col in row] for row in other.ntt_reps]
        return self.__add__(other=other)

    def __mul__(self, other):
        # num_rows_in_self: int = len(self.ntt_reps)
        # min_cols_in_self: int = min(len(x) for x in self.ntt_reps)
        # max_cols_in_self: int = max(len(x) for x in self.ntt_reps)
        # consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        # min_deg_in_self: int = min(len(x) for i in self.ntt_reps for x in i)
        # max_deg_in_self: int = max(len(x) for i in self.ntt_reps for x in i)
        # consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        # num_rows_in_other: int = len(other.ntt_reps)
        # min_cols_in_other: int = min(len(x) for x in other.ntt_reps)
        # max_cols_in_other: int = max(len(x) for x in other.ntt_reps)
        # consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
        # min_deg_in_other: int = min(len(x) for i in other.ntt_reps for x in i)
        # max_deg_in_other: int = max(len(x) for i in other.ntt_reps for x in i)
        # consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
        # same_deg: bool = max_deg_in_self == max_deg_in_other
        # if num_rows_in_self == 1 and consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == 1 and consistent_deg_in_self and consistent_deg_in_other and same_deg:
        #     return self._scalar_mul(other=other)
        # elif consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == num_rows_in_other and consistent_deg_in_self and consistent_deg_in_other and same_deg:
        #     return self._matrix_mul(other=other)
        # else:
        #     raise ValueError('Can only multiply matrices A*B where A is 1x1 and both have consistent degrees, or where A is mxn, B is nxp, and both have consistent degrees.')
        raise ValueError('Warning: You are trying to do arithmetic with a PolynomialCoefficientMatrix, which is much slower than arithmetic with a PolynomialNTTMatrix. Consider creating the associated PolynomialNTTMatrix with self.ntt() before performing arithmetic.')

    # def _scalar_mul(self, other):
    #     result = deepcopy(other)
    #     for i in range(len(other.ntt_reps)):
    #         for j in range(len(other.ntt_reps[0])):
    #             for k in range(len(self.ntt_reps[0][0])):
    #                 result.ntt_reps[i][j][k] = cent_rem(x=self.ntt_reps[0][0][k] * other.ntt_reps[i][j][k])
    #     return result
    #
    # def _matrix_mul(self, other):
    #     result = deepcopy(self)
    #     result.ntt_reps = [[[0 for k in range(len(self.ntt_reps[0][0]))] for j in range(len(other.ntt_reps[0]))] for
    #                        i in range(len(self.ntt_reps))]
    #     for i in range(len(self.ntt_reps)):
    #         for j in range(len(other.ntt_reps[0])):
    #             for k in range(len(self.ntt_reps[0][0])):
    #                 result.ntt_reps[i][j][k] = cent_rem(x=sum(
    #                     self.ntt_reps[i][l][k] * other.ntt_reps[l][j][k] for l in range(len(self.ntt_reps[0]))))
    #     return result

    def encode(self, l: int) -> bytes:
        if isinstance(l, int) and l >= 1:
            return sum(encode(x=y, l=l) for i in self.vals for y in i)
        raise ValueError()

    # def encode_ntt_reps(self, l: int) -> bytes:
    #     if isinstance(l, int) and l >= 1:
    #         return sum(encode(x=y, l=l) for i in self.ntt_reps for y in i)
    #     raise ValueError()

    def compress(self, d: int):
        if isinstance(d, int) and d >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = compress_int(d=d, x=coef)
        raise ValueError('')

    def decompress(self, d: int):
        if isinstance(d, int) and d >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = decompress_int(d=d, x=coef)
        raise ValueError('')


class PolynomialNTTMatrix(PolynomialMatrix):
    def __init__(self, modulus: int, degree: int, vals: list[list[list[int]]] | None):
        super().__init__(modulus=modulus, degree=degree, vals=vals)

    def intt(self) -> PolynomialMatrix:
        # call this to set self.coef_reps to the coefficient representation of the NTT representation
        result = []
        for row in self.vals:
            result += [[]]
            for col in row:
                result += [ntt(x=col, inv_flag=True, const_time=self.const_time)]
        return PolynomialCoefficientMatrix(modulus=self.modulus, degree=self.degree, vals=result)

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
            raise ValueError('Input matrices are not consistently sized.')
        else:
            raise ValueError('Dimension mismatch.')

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
        else:
            raise ValueError('Can only multiply matrices A*B where A is 1x1 and both have consistent degrees, or where A is mxn, B is nxp, and both have consistent degrees.')

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

    def encode(self, l: int) -> bytes:
        if isinstance(l, int) and l >= 1:
            return sum(encode(x=y, l=l) for i in self.vals for y in i)
        raise ValueError()

    def compress(self, d: int):
        if isinstance(d, int) and d >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = compress_int(d=d, x=coef)
        raise ValueError('')

    def decompress(self, d: int):
        if isinstance(d, int) and d >= 1:
            for i, row in enumerate(self.vals):
                for j, col in enumerate(row):
                    for k, coef in enumerate(col):
                        self.vals[i][j][k] = decompress_int(d=d, x=coef)
        raise ValueError('')


def decode_coefficient_matrix(x: bytes, l: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialCoefficientMatrix:
    if isinstance(x, bytes) and isinstance(l, int) and l >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1 and isinstance(modulus, int) and modulus >= 1 and len(x) >= l*num_rows*num_cols*32:
        y = deepcopy(x)
        resulting_coef_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_coef_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32*l], y[32*l:]
                w = decode(x=z, l=l)
                resulting_coef_reps[-1] += [w]
        return PolynomialCoefficientMatrix(modulus=modulus, degree=DEGREE, vals=resulting_coef_reps)
    raise ValueError('Input does not have enough bytes to decode to a PolynomialMatrix.')


def decode_ntt_matrix(x: bytes, l: int, num_rows: int, num_cols: int, modulus: int) -> PolynomialNTTMatrix:
    if isinstance(x, bytes) and isinstance(l, int) and l >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1 and isinstance(modulus, int) and modulus >= 1 and len(x) >= l*num_rows*num_cols*32*2:
        # For NTT rep, we have double the degree
        y = deepcopy(x)
        resulting_ntt_reps: list[list[list[int]]] = []
        for i in range(num_rows):
            resulting_ntt_reps += [[]]
            for j in range(num_cols):
                z, y = y[:32*l*2], y[32*l*2:]
                w = decode(x=z, l=l)
                resulting_ntt_reps[-1] += [w]
        return PolynomialNTTMatrix(modulus=modulus, degree=DEGREE, vals=resulting_ntt_reps)
    raise ValueError('Input does not have enough bytes to decode to a PolynomialMatrix.')


def XOF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        pass
    raise ValueError('Must input bytes to XOF.')


def PRF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        pass
    raise ValueError('Must input bytes to PRF.')


def KDF(x: bytes) -> bytes:
    if isinstance(x, bytes):
        pass
    raise ValueError('')


def HashFunctionH(x: bytes) -> bytes:
    if isinstance(x, bytes):
        pass
    raise ValueError('')


def HashFunctionG(x: bytes) -> bytes:
    if isinstance(x, bytes):
        pass
    raise ValueError('')


def parse(x: bytes) -> list[list[int]]:
    if isinstance(x, bytes):
        pass
    raise ValueError()


def cbd(eta: int, x: bytes) -> int:
    if isinstance(eta, int) and eta >= 1 and isinstance(x, bytes) and len(x) >= 64*eta:
        z = x.decode()
        y: str = z[:len(z)//2]
        w: str = z[len(z)//2:]
        u: int = sum(int(i == '1') for i in y)
        v: int = sum(int(i == '1') for i in w)
        return u - v
    raise ValueError()


def cpa_pke_keygen(security_parameter: int) -> tuple[bytes, bytes]:
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS:
        d: bytes = bin(randbits(256))[2:].encode()
        rho_and_sigma: bytes = HashFunctionG(d)
        rho: bytes
        sigma: bytes
        rho, sigma = rho_and_sigma[:len(rho_and_sigma)//2], rho_and_sigma[len(rho_and_sigma)//2:]
        N: int = 0
        A_ntt_reps: list[list[list[int]]] = []
        for i in range(PARAMS[security_parameter]['k']):
            A_ntt_reps += [[]]
            for j in range(PARAMS[security_parameter]['k']):
                x: bytes = rho + bin(j)[2:].encode() + bin(i)[2:].encode()
                A_ntt_reps[-1] += parse(XOF(x))
        A: PolynomialNTTMatrix = PolynomialNTTMatrix(modulus=MODULUS, degree=DEGREE, vals=A_ntt_reps)
        s_coef_reps: list[list[list[int]]] = [[] for i in range(PARAMS[security_parameter]['k'])]
        for i in range(PARAMS[security_parameter]['k']):
            next_x: bytes = sigma + bin(N)[2:].encode()
            s_coef_reps[i] = cbd(eta=PARAMS[security_parameter]['eta_one'], x=PRF(x=next_x))
            N += 1
        s: PolynomialNTTMatrix = PolynomialNTTMatrix(modulus=MODULUS, degree=DEGREE, vals=s_coef_reps)
        for i, row in enumerate(s.vals):
            for j, col in enumerate(row):
                for k, coef in enumerate(col):
                    s.vals[i][j][k] = s.vals[i][j][k] % MODULUS
        e_coef_reps: list[list[list[int]]] = [[] for i in range(PARAMS[security_parameter]['k'])]
        for i in range(PARAMS[security_parameter]['k']):
            next_x: bytes = sigma + bin(N)[2:].encode()
            e_coef_reps[i] = cbd(eta=PARAMS[security_parameter]['k'], x=next_x)
            N += 1
        e: PolynomialNTTMatrix = PolynomialNTTMatrix(modulus=MODULUS, degree=DEGREE, vals=e_coef_reps)
        t: PolynomialNTTMatrix = A * s + e
        for i, row in enumerate(t.vals):
            for j, col in enumerate(row):
                for k, coef in enumerate(col):
                    t.vals[i][j][k] = t.vals[i][j][k] % MODULUS
        pk: bytes = t.encode(l=12) + rho
        sk: bytes = s.encode(l=12)
        return pk, sk
    raise ValueError(
        f'Must have security_parameter={security_parameter} in ALLOWABLE_SECURITY_PARAMETERS=' +
        f'{ALLOWABLE_SECURITY_PARAMETERS}.')


def cpa_pke_encrypt(security_parameter: int, pk: bytes, m: bytes) -> tuple[bytes, bytes]:
    if security_parameter in ALLOWABLE_SECURITY_PARAMETERS:
        pass
    raise ValueError()


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
