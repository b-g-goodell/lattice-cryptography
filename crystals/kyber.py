"""
A python implementation of the CRYSTALS-Kyber algorithm.

Documentation
-------------
Documentation is hosted at ???

Hosting
-------
The repository for this project can be found at ???

License
-------
Released under the MIT License, see LICENSE file for details. ???

Copyright
---------
Copyright (c) 2022 Geometry Labs, Inc.
Funded by The QRL Foundation.

Contributors
------------
Brandon Goodell (lead author), Mitchell Krawiec-Thayer.
"""
from math import ceil, floor, log2
from copy import deepcopy
from secrets import randbits

# Parameters from Kyber1024
N: int = 256
Q: int = 3329
K: int = 4
ETA: int = 2
D_U: int = 11
D_V: int = 5
ROU: int = 17
ROU_INVERSE: int = 1175
SECBITS: int = 174  # in a sense, "target bits of security," although Kyber1024 claims only 128 bits of sec.
SEED_LEN: int = 32


# Convenient constants derived from the parameters above
TWICE_N: int = 2 * N
LOG_TWICE_N: int = ceil(log2(TWICE_N))
HALF_Q: int = Q // 2
LOG_Q: int = ceil(log2(Q))
LOG_K: int = ceil(log2(K))
ZETAS: list[int] = [(ROU ** i) % Q for i in [TWICE_N // (2 ** (j + 1)) for j in range(LOG_TWICE_N)]]
ZETA_INVERSES: list[int] = [(ROU_INVERSE ** i) % Q for i in [TWICE_N // (2 ** (j + 1)) for j in range(LOG_TWICE_N)]]

CPA_PKE_SK_LEN: int = LOG_Q * K * N // 8
CPA_PKE_PK_LEN: int = LOG_Q * K * N // 8 + SEED_LEN
ENCODED_CPA_PKE_SK_LEN: int = CPA_PKE_SK_LEN // (K * N // 8)
ENCODED_CPA_PKE_PK_LEN: int = (CPA_PKE_PK_LEN - SEED_LEN) // (K * N // 8)
CPA_PKE_FIRST_CIPHERTEXT_LEN: int = D_U * K * N // 8
CPA_PKE_SECOND_CIPHERTEXT_LEN: int = D_V * N // 8
CPA_PKE_CIPHERTEXT_LEN: int = CPA_PKE_FIRST_CIPHERTEXT_LEN + CPA_PKE_SECOND_CIPHERTEXT_LEN
CCA_KEM_SK_LEN: int = CPA_PKE_SK_LEN + CPA_PKE_PK_LEN + 2*SEED_LEN
CCA_KEM_PK_LEN: int = CPA_PKE_PK_LEN
CCA_KEM_CIPHERTEXT_LEN: int = CPA_PKE_CIPHERTEXT_LEN


def _int2bytes(x: int, length: int) -> bytes:
    return bin(x)[2:].encode().zfill(length)


def int2bytes(x: int, length: int = LOG_K) -> bytes:
    """
    Essentially the encode function from the specs for a single integer: write the binary expansion of 'x', fill out
    to take up 'length' bits. NOTE: output of this function is a bytes object consisting of 'length//8' bytes.

    :param x: Input integer to be encoded as bytes.
    :type x: int
    :param length: Length
    :type length: int

    :return: Bytes representation of input x
    :rtype: bytes
    """
    if isinstance(x, int) and isinstance(length, int) and 0 <= x < 2**length:
        return _int2bytes(x=x, length=length)
    elif not isinstance(x, int):
        raise TypeError(f'Cannot int2bytes with x, length unless x is an integer, but had type(x)={type(x)}.')
    elif not isinstance(length, int):
        raise TypeError(f'Cannot int2bytes with x, length unless length is an integer, but had type(length)={type(length)}.')
    raise ValueError(f'Cannot int2bytes with x, length unless 0 <= x < {2**length}, but had x={x}.')


def _bytes2int(x: bytes) -> int:
    return int(x, 2)


def bytes2int(x: bytes) -> int:
    """
    Essentially the decode function from the specs for a single integer, interprets the input bytes object as the binary
    expansion of an integer in bits but written as bytes.

    :param x: Input integer to be encoded as bytes.
    :type x: int
    :return: Integer representation of input x
    :rtype: int
    """
    if isinstance(x, bytes):
        return _bytes2int(x=x)
    raise TypeError(f'Cannot bytes2int with x unless x is a bytes object, but had type(x)={type(x)}.')


def _bit_rev(x: int, length: int) -> int:
    return bytes2int(x=int2bytes(x=x, length=length)[::-1])


def bit_rev(x: int, length: int) -> int:
    """
    Reverse the bits in the binary expansion of x
    :param x: Input integer whose bits are to be reversed
    :type x: int
    :param length: Length of binary expansion
    :type length: int
    :return: Integer whose binary expansion is the reverse of the input x
    :rtype: int
    """
    if isinstance(length, int) and length >= 1 and isinstance(x, int) and 0 <= x < 2 ** length:
        return _bit_rev(x=x, length=length)
    elif not isinstance(x, int):
        raise TypeError(f'Cannot bit_rev with x, length unless x is an integer, but had type(x)={type(x)}.')
    elif not isinstance(length, int):
        raise TypeError(f'Cannot bit_rev with x, length unless length is an integer, but had type(length)={type(length)}.')
    elif length < 1:
        raise ValueError(f'Cannot bit_rev with x, length unless length >= 1, but had length={length}.')
    raise ValueError(f'Cannot bit_rev with x, length unless 0 <= x < {2**length} but had x={x}.')


def is_pow_two(x: int) -> bool:
    """
    Check whether input integer is a (positive) power-of-two.
    :param x: Input integer
    :type x: int
    :return: Boolean indicating whether x is a positive power of two.
    :rtype: bool
    """
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    return False


def _bit_rev_cp(x: list[int], num_bits: int) -> list[int]:
    return [x[bit_rev(x=i, length=num_bits)] for i in range(len(x))]


def bit_rev_cp(x: list[int], num_bits: int) -> list[int]:
    """
    Input a list of integers with power-of-two length and permute by reversing the digits in the binary expansions of
    the indices. For example: if x = [172, 31, 56, 7], the indices are 0, 1, 2, 3, whose binary expansions are
    00, 01, 10, 11. When reversed digit-by-digit, the binary expansions are 00, 10, 01, 11, or 0, 2, 1, 3. So bit_rev_cp
    will output [x[0], x[2], x[1], x[3]] = [172, 56, 31, 7]. For another example: if
    x = [172, 31, 56, 7, 202, 63, 17, 13], the indices are 0, 1, 2, 3, 4, 5, 6, 7, whose binary expansions are 000, 001,
    010, 011, 100, 101, 110, 111. When reversed, the binary expansions are 000, 100, 010, 110, 001, 101, 011, 111, or
    0, 4, 2, 6, 1, 5, 3, 7. So bit_rev_cp will output
    [x[0], x[4], x[2], x[6], x[1], x[5], x[3], x[7] = [172, 202, 56, 17, 31, 63, 56, 7, 13].
    :param x: List of integers
    :type x: list[int]
    :param num_bits: Number of bits required to describe the maximum index in x.
    :type x: int
    :return: List of integers with bit-reversed indices.
    :rtype: list[int]
    """
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(num_bits, int) and num_bits >= 1 and len(x) == 2 ** num_bits:
        return _bit_rev_cp(x=x, num_bits=num_bits)
    elif not isinstance(x, list):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless x is a list, but had type(x)={type(x)}.')
    elif not all(isinstance(y, int) for y in x):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless x is a list of integers, but had (type(y) for y in x)={(type(y) for y in x)}.')
    elif not isinstance(num_bits, int):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless num_bits is an integer, but had type(num_bits)={type(num_bits)}.')
    elif num_bits < 1:
        raise ValueError(f'Cannot bit_rev_cp with x, num_bits unless num_bits >= 1, but had num_bits={num_bits}.')
    raise ValueError(f'Cannot bit_rev_cp with x, num_bits unless len(x) == {2**num_bits} but had len(x)={len(x)}.')


def _reduce(x: int) -> int:
    y: int = x % Q
    z: int = y - HALF_Q - 1
    w: int = y - (1 + (z >> LOG_Q)) * Q
    return w


def reduce(x: int) -> int:
    """
    Compute some integer w such that -MODULUS//2 <= w <= MODULUS//2 and such that (x-w) % MODULUS == 0, in constant
    time.
    :param x: Input integer
    :type x: int
    :return: "Centered" representative of x % MODULUS.
    :rtype: int
    """
    if isinstance(x, int):
        return _reduce(x=x)
    raise TypeError(f'Cannot compute reduce for x unless x is an integer, but had type(x)={type(x)}.')


def _round_up(x: float | int) -> int:
    ceil_x: int = ceil(x)
    if ceil_x - x <= 0.5:
        return ceil_x
    return floor(x)


def round_up(x: float | int) -> int:
    """
    Round an input float (or integer) to the nearest integer, so that a float ending with the decimal .5 rounds up.

    :param x: Input number
    :type x: float or integer
    :return: Integer closest to x, (rounding up for floats ending in .5).
    :rtype: int
    """
    if isinstance(x, float) or isinstance(x, int):
        return _round_up(x=x)
    raise TypeError(f'Cannot round_up with x unless x is a float or an int, but had type(x)={type(x)}.')


def _parse_one(x: bytes) -> tuple[list[int], int]:
    """
    Parse an input bytes object into a 2-tuple, where the first entry is a list of N=256 integers in the list
    [0, 1, 2, 3, ..., MODULUS - 1] for MODULUS=3329, and the second entry is an integer describing the index of the
    first unused byte in x.
    TODO: Extend this to work for arbitrary N and arbitrary MODULUS.
    :param x: Input bytes
    :type x: bytes
    :return: A 2-tuple, where the first entry is a list of integers, and the second entry is an integer.
    :rtype: tuple[list[int], int]
    """
    i: int = 0
    j: int = 0
    result: list[int] = []
    while j < N and i+2 < len(x):
        d1 = x[i] + 256 * (x[i + 1] % 16)
        d2 = (x[i + 1] // 16) + 16 * x[i + 2]
        if d1 < Q:
            result += [d1]
            j += 1
        if d2 < Q and j < N:
            result += [d2]
            j += 1
        i += 3
    if len(result) < N:
        raise RuntimeError(f'Parsing failed! Did not have enough input bits to do the job.')
    return result, i


def _parse_many(x: bytes) -> list[list[list[int]]]:
    """
    Parse an input bytes object into a list of K=4 lists of K=4 lists of N=256 integers, each in the list
    [0, 1, 2, 3, ..., MODULUS-1] for MODULUS=3329, by calling _parse_one several times.
    TODO: Extend this to work for arbitrary K, arbitrary N, and arbitrary MODULUS.
    :param x: Input bytes
    :type x: bytes
    :return: A list of lists of lists of integers.
    :rtype: list[list[list[int]]]
    """
    result: list[list[list[int]]] = []
    starting_index: int = 0
    for i in range(K):
        result += [[]]
        for j in range(K):
            next_result, starting_index = _parse_one(x=x[starting_index:])
            result[-1] += [next_result]
    return result


def parse(x: bytes) -> list[int] | list[list[list[int]]]:
    """
    Parse an input bytes object into either a list of 256 integers, or a list of K=4 lists of K=4 lists of N=256
    integers, each in the list [0, 1, 2, 3, ..., MODULUS-1] for MODULUS=3329.

    Our implementation does NOT parse according to specifications, because we only input a list of bytes, not a bytestream.
    Our implementation is not compatible with the NIST standard: users will observe a different "A" matrix.
    There is a small positive probability that parsing fails.

    :param x: Input bytes
    :type x: bytes
    :return: A list of lists of lists of integers.
    :rtype: list[list[list[int]]]
    """
    if isinstance(x, bytes) and len(x) == N*(LOG_Q + SECBITS)//8:
        return _parse_one(x=x)[0]
    elif isinstance(x, bytes) and len(x) >= K*K*N*(LOG_Q + SECBITS)//8:
        return _parse_many(x=x)
    elif not isinstance(x, bytes):
        raise TypeError(f'Cannot parse with x unless x is a bytes object, but had type(x)={type(x)}.')
    raise ValueError(f'Cannot parse with x unless x is a bytes object with length {N * (LOG_Q + SECBITS) // 8} or length at least {K * K * N * (LOG_Q + SECBITS) // 8} but had len(x)={len(x)}.')


# def is_arithmetic_legal(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]]) -> bool:
#     num_rows_in_self: int = len(a_vals)
#     num_rows_in_other: int = len(b_vals)
#
#     min_cols_in_self: int = min(len(x) for x in a_vals)
#     min_cols_in_other: int = min(len(x) for x in b_vals)
#     max_cols_in_self: int = max(len(x) for x in a_vals)
#     max_cols_in_other: int = max(len(x) for x in b_vals)
#     consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
#     consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
#
#     min_deg_in_self: int = min(len(x) for i in a_vals for x in i)
#     min_deg_in_other: int = min(len(x) for i in b_vals for x in i)
#     max_deg_in_self: int = max(len(x) for i in a_vals for x in i)
#     max_deg_in_other: int = max(len(x) for i in b_vals for x in i)
#     consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
#     consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
#
#     same_rows: bool = num_rows_in_self == num_rows_in_other
#     same_cols: bool = max_cols_in_self == max_cols_in_other
#     same_deg: bool = max_deg_in_self == max_deg_in_other
#
#     return same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg
#
#
# def add(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]]) -> list[list[list[int]]]:
#     result = deepcopy(a_vals)
#     for i, row in enumerate(b_vals):
#         for j, col in enumerate(row):
#             for k, x in enumerate(col):
#                 result[i][j][k] = reduce(x=result[i][j][k] + x)
#     return result
#
#
# def mul(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]]) -> list[list[list[int]]]:
#     result = deepcopy(a_vals)
#     for i, row in enumerate(b_vals):
#         for j, col in enumerate(row):
#             for k, x in enumerate(col):
#                 result[i][j][k] = reduce(x=result[i][j][k] * x)
#     return result
#
#
# class PolyCoefs(object):
#     """
#     Class for coefficient representations of matrices of polynomials. This class does not support arithmetic.
#     # TODO: Implement addition
#     # TODO: Add a multiplication method but throw a NotImplementedError
#
#     Attributes
#     ----------
#         modulus: int
#             Integer modulus for all coefficients (0 <= coefficient < modulus)
#         degree: int
#             Integer degree for all polynomials
#         halfmod: int
#             Half of modulus, rounded down.
#         logmod: int
#             Number of bits to describe a coefficient.
#         vals: list[list[list[int]]]
#             2-dimensional matrix with list[int] entries, each entry is a coefficient representation of a polynomial.
#         const_time: bool
#             Flag for describing whether arithmetic is constant-time - in prod, always keep at True.
#
#     Methods
#     -------
#         __init__(self, vals, modulus, degree)
#     """
#     modulus: int = Q
#     degree: int = N
#     halfmod: int = HALF_MODULUS
#     logmod: int = LOG_MODULUS
#     vals: list[list[list[int]]] = []
#     const_time: bool = True
#
#     def __init__(self, vals: list[list[list[int]]] | None, modulus: int = Q, degree: int = N, const_time: bool = True):
#         self.modulus = modulus
#         self.degree = degree
#         self.halfmod = modulus//2
#         self.logmod = ceil(log2(modulus))
#         self.vals = vals
#         self.const_time = const_time
#
#     def __add__(self, other):
#         if is_arithmetic_legal(a_vals=self.vals, b_vals=other.vals):
#             result = deepcopy(self)
#             result.vals = add(a_vals=self.vals, b_vals=other.vals)
#             return result
#         raise ValueError(f'Cannot compute PolyCoefs.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')
#
#     def __radd__(self, other):
#         if other == 0:
#             return self
#         return self.__add__(other=other)
#
#     def __sub__(self, other):
#         negative_other = deepcopy(other)
#         negative_other.vals = [[[-coef for coef in col] for col in row] for row in other.vals]
#         return self.__add__(other=other)
#
#     def __mul__(self, other):
#         raise NotImplementedError(f'Multiplication of PolyCoefs is not implemented. Apply the NTT and multiply the resulting PolyNTT objects instead.')
#
#
# class PolyNTT(object):
#     """
#     Class for coefficient representations of matrices of polynomials. This class does not support arithmetic.
#
#     Attributes
#     ----------
#         modulus: int
#             Integer modulus for all coefficients (0 <= coefficient < modulus)
#         degree: int
#             Integer degree for all polynomials
#         halfmod: int
#             Half of modulus, rounded down.
#         logmod: int
#             Number of bits to describe a coefficient.
#         vals: list[list[list[int]]]
#             2-dimensional matrix with list[int] entries, each entry is a coefficient representation of a polynomial.
#         const_time: bool
#             Flag for describing whether arithmetic is constant-time - in prod, always keep at True.
#
#     Methods
#     -------
#         __init__(self, vals, modulus, degree, const_time)
#             Set self.modulus to modulus, self.degree to degree, compute self.halfmod and self.logmod, set self.vals to vals, and set self.const_time to const_time
#         __add__(self, other)
#             Checks if dimensions of self and other are compatible, and then return the sum of self and other.
#         __radd__(self, other)
#             If other == 0, return self, other call __add__. This way, we can use sum().
#         __sub__(self, other)
#             Deepcopy other, negate each entry in vals, and then call __add__.
#         __mul__(self, other)
#             Check if self is 1-by-1 or if self and other have dimensions match. Then, if self is 1-by-1, call _scalar_mul, otherwise call _matrix_mul.
#         _scalar_mul(self, other)
#             Multiply each polynomial in other by self.
#         _matrix_mul(self, other)
#             Multiply the matrices self and other coordinate-wise.
#     """
#     modulus: int = Q
#     degree: int = N
#     halfmod: int = HALF_MODULUS
#     logmod: int = LOG_MODULUS
#     vals: list[list[list[int]]] = []
#     const_time: bool = True
#
#     def __init__(self, vals: list[list[list[int]]] | None, modulus: int = Q, degree: int = N, const_time: bool = True):
#         self.modulus = modulus
#         self.degree = degree
#         self.halfmod = modulus//2
#         self.logmod = ceil(log2(modulus))
#         self.vals = vals
#         self.const_time = const_time
#
#     def __add__(self, other):
#         if is_arithmetic_legal(a_vals=self.vals, b_vals=other.vals):
#             result = deepcopy(self)
#             result.vals = add(a_vals=self.vals, b_vals=other.vals)
#             return result
#         raise ValueError(f'Cannot compute PolyNTT.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')
#
#     def __radd__(self, other):
#         if other == 0:
#             return self
#         return self.__add__(other=other)
#
#     def __sub__(self, other):
#         negative_other = deepcopy(other)
#         negative_other.vals = [[[-coef for coef in col] for col in row] for row in other.vals]
#         return self.__add__(other=other)
#
#     def __mul__(self, other):
#         num_rows_in_self: int = len(self.vals)
#         num_rows_in_other: int = len(other.vals)
#
#         min_cols_in_self: int = min(len(x) for x in self.vals)
#         min_cols_in_other: int = min(len(x) for x in other.vals)
#         max_cols_in_self: int = max(len(x) for x in self.vals)
#         max_cols_in_other: int = max(len(x) for x in other.vals)
#         consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
#         consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other
#
#         min_deg_in_self: int = min(len(x) for i in self.vals for x in i)
#         min_deg_in_other: int = min(len(x) for i in other.vals for x in i)
#         max_deg_in_self: int = max(len(x) for i in self.vals for x in i)
#         max_deg_in_other: int = max(len(x) for i in other.vals for x in i)
#         consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
#         consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other
#
#         # same_rows: bool = num_rows_in_self == num_rows_in_other
#         # same_cols: bool = max_cols_in_self == max_cols_in_other
#         same_deg: bool = max_deg_in_self == max_deg_in_other
#
#         if num_rows_in_self == 1 and consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == 1 and consistent_deg_in_self and consistent_deg_in_other and same_deg:
#             return self._scalar_mul(other=other)
#         elif consistent_cols_in_self and consistent_cols_in_other and max_cols_in_self == num_rows_in_other and consistent_deg_in_self and consistent_deg_in_other and same_deg:
#             return self._matrix_mul(other=other)
#         raise ValueError(
#             'Cannot compute PolynomialNTTMatrix.__mul__ with unless self is 1x1 and both have consistent degrees, or where self is mxn, other is nxp, and both have consistent degrees (dim mismatch).')
#
#     def _scalar_mul(self, other):
#         result = deepcopy(other)
#         for i in range(len(other.vals)):
#             for j in range(len(other.vals[0])):
#                 for k in range(len(self.vals[0][0])):
#                     result.ntt_reps[i][j][k] = reduce(x=self.vals[0][0][k] * other.vals[i][j][k])
#         return result
#
#     def _matrix_mul(self, other):
#         result = deepcopy(self)
#         result.ntt_reps = [[[0]*len(self.vals[0][0])]*len(other.vals[0])]*len(self.vals)
#         for i in range(len(self.vals)):
#             for j in range(len(other.vals[0])):
#                 for k in range(len(self.vals[0][0])):
#                     result.ntt_reps[i][j][k] = reduce(x=sum(self.vals[i][l][k] * other.ntt_reps[l][j][k] for l in range(len(self.vals[0]))))
#         return result
#
#
# def _cbd_eta(x: bytes, eta: int = ETA) -> list[int]:
#     x_as_bits: list[int] = [int(y == '1') for y in x.decode()]
#     result: list[int] = [0 for _ in range(N)]
#     for i in range(N):
#         a: int = sum([x_as_bits[2*i*eta + j] for j in range(eta)])
#         b: int = sum([x_as_bits[2*i*eta + eta + j] for j in range(eta)])
#         result[i] = a - b
#     return result
#
#
# def cbd_eta(x: bytes, eta: int = ETA) -> list[int]:
#     """
#     Sample a list of integers with length N from input bytes x, such that the integers are sampled from the centered
#     binomial distribution (CBD) with support [-eta, -eta + 1, ..., eta - 1, eta]. Works by looking at x as a bitstring,
#     and taking the difference of sums of sections of the bit string.
#
#     :param x: Input bytes
#     :type x: bytes
#     :param eta: Input bound.
#     :type eta: int
#
#     :return: List of integers.
#     :rtype: list[int]
#     """
#     if isinstance(x, bytes) and len(x) >= 2*N*eta//8 and isinstance(eta, int):
#         return _cbd_eta(x=x, eta=eta)
#     elif not isinstance(x, bytes):
#         raise TypeError(f'Cannot cbd_eta with x, eta unless x is a bytes object, but type(x)={type(x)}.')
#     elif len(x) < 2*N*eta//8:
#         raise ValueError(f'Cannot cbd_eta with x, eta unless len(x)>={2*N*eta//8} but had len(x)={len(x)}.')
#     raise TypeError(f'Cannot cbd_eta with x, eta unless eta is an integer, but had type(x)={type(x)}.')
#
#
# def _cbd_polycoefs(x: bytes, eta: int = ETA, num_rows: int = K, num_cols: int = 1) -> PolyCoefs:
#     vals: list[list[list[int]]] = []
#     for i in range(num_rows):
#         next_x: bytes = x[i*num_cols*2*N*eta//8: (i+1)*num_cols*2*N*eta//8]
#         vals += [[]]
#         for j in range(num_cols):
#             vals[-1] += [cbd_eta(x=next_x[j * (2 * N * eta // 8): (j + 1) * (2 * N * eta // 8)], eta=eta)]
#     return PolyCoefs(vals=vals)
#
#
# def cbd_polycoefs(x: bytes, eta: int = ETA, num_rows: int = K, num_cols: int = 1) -> PolyCoefs:
#     """
#     Sample a matrix of polynomials with num_rows rows, num_cols columns, degrees N, and such that each coefficient is in
#     the list [-eta, -eta + 1, ..., eta - 1, eta]. Works by calling cbd_eta for each polynomial in the matrix.
#
#     :param x: Input bytes
#     :type x: bytes
#     :param eta: Integer bound eta
#     :type eta: int
#     :param num_rows: Number of rows.
#     :type num_rows: int
#     :param num_cols: Number of columns.
#     :type num_cols: int
#
#     :return: PolyCoefs object
#     :rtype: PolyCoefs
#     """
#     if isinstance(x, bytes) and len(x) >= num_rows*num_cols*N*2*eta//8 and isinstance(eta, int) and eta >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1:
#         return _cbd_polycoefs(x=x, eta=eta, num_rows=num_rows, num_cols=num_cols)
#     elif not isinstance(x, bytes):
#         raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless x is a bytes object, but had type(x)={type(x)}.')
#     elif len(x) < num_rows*num_cols*N*2*eta//8:
#         raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless len(x) >= {num_rows*num_cols*N*2*eta//8} but had len(x)={len(x)}.')
#     elif not isinstance(eta, int):
#         raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless eta is an integer, but had type(eta)={type(eta)}.')
#     elif eta < 1:
#         raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless eta >= 1, but had eta={eta}.')
#     elif not isinstance(num_rows, int):
#         raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_rows is an integer, but had type(num_rows)={type(num_rows)}.')
#     elif num_rows < 1:
#         raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_rows >= 1, but had num_rows={num_rows}.')
#     elif not isinstance(num_cols, int):
#         raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_cols is an integer, but had type(num_cols)={type(num_cols)}.')
#     raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_cols >= 1, but had num_cols={num_cols}.')
#
#
# def _compress_one_int(x: int, d: int, p: int = Q) -> int:
#     # We rename q to p to avoid ambiguity with the constant Q
#     return round_up(x=x * 2 ** d / p) % 2 ** d
#
#
# def _compress_list_of_ints(x: list[int], d: int, p: int = Q) -> list[int]:
#     return [_compress_one_int(x=y, d=d, p=p) for y in x]
#
#
# def _compress_many_ints(x: list[list[list[int]]], d: int, p: int = Q) -> list[list[list[int]]]:
#     return [[_compress_list_of_ints(x=z, d=d, p=p) for z in y] for y in x]
#
#
# def _compress_polycoefs(x: PolyCoefs, d: int, p: int = Q) -> PolyCoefs:
#     return PolyCoefs(vals=_compress_many_ints(x=x.vals, d=d, p=p), modulus=x.modulus, degree=x.degree)
#
#
# def _should_compress_many(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> bool:
#     return isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2
#
#
# def compress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
#     """
#     Compresses an integer modulo p (or a list of them, or a list of lists of lists of them) or the integers modulo p in
#     a PolyCoefs object to a d-bit integer as specified. Works by, essentially, looking at x/p as a ratio, and going that
#     far along the list [0, 1, 2, ..., 2**d - 1].
#
#     :param x: Input data
#     :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
#     :param d: Input number of bits for compressed data
#     :type d: int
#     :param p: Input modulus
#     :type p: int
#
#     :return: Compressed data
#     :rtype: int | list[int] | list[list[list[int]]] | PolyCoefs
#
#     """
#     if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _compress_one_int(x=x, d=d, p=p)
#     elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _compress_list_of_ints(x=x, d=d, p=p)
#     elif _should_compress_many(x=x, d=d, p=p):
#         return _compress_many_ints(x=x, d=d, p=p)
#     elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _compress_polycoefs(x=x, d=d, p=p)
#     raise ValueError(f'Cannot compute compress for x, d, p unless x is an integer, a list of integers, or a list of lists of lists of integers, or a PolyCoefs... and d is an integer >= 1 and p is an integer >= 2, but had (type(x), d, p)={(type(x), d, p)}.')
#
#
# def _decompress_one_int(x: int, d: int, p: int = Q) -> int:
#     return round_up(x=x * p / 2 ** d)
#
#
# def _decompress_list_of_ints(x: list[int], d: int, p: int = Q) -> list[int]:
#     return [_decompress_one_int(x=y, d=d, p=p) for y in x]
#
#
# def _decompress_many_ints(x: list[list[list[int]]], d: int, p: int = Q) -> list[list[list[int]]]:
#     return [[_decompress_list_of_ints(x=z, d=d, p=p) for z in y] for y in x]
#
#
# def _decompress_polycoefs(x: PolyCoefs, d: int, p: int = Q) -> PolyCoefs:
#     return PolyCoefs(vals=_decompress_many_ints(x=x.vals, d=d, p=p), modulus=x.modulus, degree=x.degree)
#
#
# def decompress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
#     """
#     Decompresses an integer modulo 2**d (or a list of them, or a list of lists of lists of them) or the integers modulo
#     2**d in a PolyCoefs object to an integer modulo p as specified. Works by, essentially, looking at x/2**d as a ratio,
#     and going that far along the list [0, 1, 2, ..., p-1].
#
#     :param x: Input data
#     :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
#     :param d: Input number of bits for compressed data
#     :type d: int
#     :param p: Input modulus
#     :type p: int
#
#     :return: Decompressed data
#     :rtype: int | list[int] | list[list[list[int]]] | PolyCoefs
#     """
#     if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _decompress_one_int(x=x, d=d, p=p)
#     elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _decompress_list_of_ints(x=x, d=d, p=p)
#     elif _should_compress_many(x=x, d=d, p=p):
#         return _decompress_many_ints(x=x, d=d, p=p)
#     elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
#         return _decompress_polycoefs(x=x, d=d, p=p)
#     raise ValueError(f'Cannot decompress with x, d, p unless x is an integer, or a list of integers, or a list of lists of lists of integers, or a PolyCoefs... and d is an integer >= 1, and p is an integer >= 2, but had (type(x), d, p)={(type(x), d, p)}.')
#
#
# def _encode_m_one_int(x: int, m: int) -> bytes:
#     return int2bytes(x=x, length=m)  # bin(x)[2:].zfill(m).encode()
#
#
# def _encode_m_list_of_ints(x: list[int], m: int) -> bytes:
#     result = bytes(0)
#     for y in x:
#         result += _encode_m_one_int(x=y, m=m)
#     return result
#
#
# def _encode_m_many_ints(x: list[list[list[int]]], m: int) -> bytes:
#     result = bytes(0)
#     for y in x:
#         for z in y:
#             result += _encode_m_list_of_ints(x=z, m=m)
#     return result
#
#
# def _encode_m_matrix(x: PolyCoefs | PolyNTT, m: int) -> bytes:
#     return _encode_m_many_ints(x=x.vals, m=m)
#
#
# def encode_m(x: int | list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT, m: int) -> bytes:
#     """
#     We rename encode_l to encode_m to use m instead of l because l is an ambiguous character. Encodes an integer (or a
#     list of integers, or a list of lists of lists of integers) or the integers in a PolyCoefs or PolyNTT object to bytes. Works
#     with the usual serialization... each piece of input data is an integer modulo p, and we just pad out the binary
#     expansion of the input data to m bytes.
#
#     :param x: Input data
#     :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
#     :param m: Number of bytes
#     :type m: int
#
#     :return: Encoded data
#     :rtype: bytes
#     """
#     if isinstance(m, int) and m >= 1 and isinstance(x, int) and 0 <= x < 2 ** m:
#         return _encode_m_one_int(x=x, m=m)
#     elif isinstance(m, int) and m >= 1 and isinstance(x, list) and all(isinstance(y, int) for y in x) and len(x) == N and all(0 <= y < 2 ** m for y in x):
#         return _encode_m_list_of_ints(x=x, m=m)
#     elif isinstance(m, int) and m >= 1 and isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and all(0 <= w < 2**m for y in x for z in y for w in z):
#         return _encode_m_many_ints(x=x, m=m)
#     elif isinstance(m, int) and m >= 1 and (isinstance(x, PolyCoefs) or isinstance(x, PolyNTT)):
#         return _encode_m_matrix(x=x, m=m)
#     raise ValueError(f'Cannot compute encode for (x, m) unless m >= 1 is an integer and x is an m-bit integer, or a list of m-bit integers, or x is a list of lists of lists of m-bit integers, or x is a PolyCoefs, but had (type(x), m)={(type(x), m)}.')
#
#
# def _decode_m_one_int(x: bytes) -> int:
#     # TODO: WRONG
#     # light wrapper for bytes2int
#     return bytes2int(x=x)
#
#
# def _decode_m_list_of_ints(x: bytes, m: int) -> list[int]:
#     # TODO: WRONG
#     return [_decode_m_one_int(x=x[i * m: (i + 1) * m]) for i in range(len(x) // m)]
#
#
# def _decode_m_many(x: bytes, m: int) -> list[list[list[int]]]:
#     # TODO: WRONG
#     result: list[list[list[int]]] = []
#     y: list[bytes] = [x[i*N*m: (i+1)*N*m] for i in range(K)]
#     for next_row in y:
#         next_poly: list[int] = _decode_m_list_of_ints(x=next_row, m=m)
#         result += [[next_poly]]
#     return result
#
#
# def _decode_m_matrix(x: bytes, m: int, ntt_matrix_flag: bool) -> PolyCoefs | PolyNTT:
#     # TODO: WRONG
#     if ntt_matrix_flag:
#         return PolyNTT(vals=_decode_m_many(x=x, m=m))
#     return PolyCoefs(vals=_decode_m_many(x=x, m=m))
#
#
# def decode_m(x: bytes, m: int, coef_matrix_flag: bool = False, ntt_matrix_flag: bool = False) -> list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT:
#     """
#     TODO: WRONG
#     We rename decode_l to decode_m to use m instead of l because l is an ambiguous character. Decodes a bytes object to
#     an integer (or a list of integers, or a list of lists of lists of integers) or a PolyCoefs or PolyNTT object. Works
#     with the usual serialization... interpret each chunk of m bytes as determined by the binary expansion of an integer,
#     written as bytes.
#
#     :param x: Input data
#     :type x: bytes
#     :param m: Number of bytes per integer
#     :type m: int
#     :param coef_matrix_flag: Flag indicating whether input x should be decoded to a PolyCoefs object.
#     :type coef_matrix_flag: bool
#     :param ntt_matrix_flag: Flag indicating whether input x should be decoded to a PolyNTT object.
#     :type ntt_matrix_flag: bool
#
#     :return: Encoded data
#     :rtype: bytes
#     """
#     if isinstance(x, bytes) and len(x) == N*m:
#         return _decode_m_list_of_ints(x=x, m=m)
#     elif isinstance(x, bytes) and len(x) >= K*N*m and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and not coef_matrix_flag and not ntt_matrix_flag:
#         return _decode_m_many(x=x, m=m)
#     elif isinstance(x, bytes) and len(x) >= K*N*m and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and (coef_matrix_flag ^ ntt_matrix_flag):
#         return _decode_m_matrix(x=x, m=m, ntt_matrix_flag=ntt_matrix_flag)
#     raise ValueError(f'Cannot compute decode_m for (x, m, coef_matrix_flag, ntt_matrix_flag) unless (i) x is a bytes object of length {N*m} or (ii) x is a bytes object with length at least {K*N*m} and both coef_matrix_flag and ntt_matrix_flag are booleans with at most one of them True, but had (type(x), len(x), type(coef_matrix_flag), coef_matrix_flag, type(ntt_matrix_flag), ntt_matrix_flag)={(type(x), len(x), type(coef_matrix_flag), coef_matrix_flag, type(ntt_matrix_flag), ntt_matrix_flag)}.')
#
#
# def _ntt_one(x: list[int], inv_flag: bool, const_time: bool = True) -> list[int]:
#     bit_rev_x: list[int] = bit_rev_cp(x=x, num_bits=ceil(log2(len(x))))
#     m: int = 1
#     for s in range(1, LOG_TWICE_DEGREE + 1):
#         m *= 2
#         if inv_flag:
#             this_zeta: int = ZETA_INVERSES[s - 1]
#         else:
#             this_zeta: int = ZETAS[s - 1]
#         for k in range(0, TWICE_DEGREE, m):
#             w: int = 1
#             for j in range(m // 2):
#                 t: int = w * bit_rev_x[k + j + m // 2]
#                 u: int = bit_rev_x[k + j]
#                 if const_time:
#                     bit_rev_x[k + j]: int = reduce(x=u + t)
#                     bit_rev_x[k + j + m // 2]: int = reduce(x=u - t)
#                 else:
#                     bit_rev_x[k + j]: int = (u + t) % Q
#                     if bit_rev_x[k + j] > HALF_MODULUS:
#                         bit_rev_x[k + j] = bit_rev_x[k + j] - Q
#                     bit_rev_x[k + j + m // 2]: int = (u - t) % Q
#                     if bit_rev_x[k + j + m // 2] > HALF_MODULUS:
#                         bit_rev_x[k + j + m // 2] = bit_rev_x[k + j + m // 2] - Q
#                 w *= this_zeta
#     if inv_flag:
#         n_inv: int = 1
#         while (n_inv * TWICE_DEGREE) % Q != 1:
#             n_inv += 1
#         if const_time:
#             bit_rev_x: list[int] = [reduce(x=(n_inv * i)) for i in bit_rev_x]
#         else:
#             bit_rev_x: list[int] = [(n_inv * i) % Q for i in bit_rev_x]
#             bit_rev_x = [i if i <= HALF_MODULUS else i - Q for i in bit_rev_x]
#     return bit_rev_x
#
#
# def _ntt_many(x: list[list[list[int]]], inv_flag: bool, const_time: bool = True) -> list[list[list[int]]]:
#     return [[_ntt_one(x=z, inv_flag=inv_flag, const_time=const_time) for z in y] for y in x]
#
#
# def _ntt_raw(x: list[int] | list[list[list[int]]], inv_flag: bool, const_time: bool = True) -> list[int] | list[list[list[int]]]:
#     if isinstance(x, list) and is_pow_two(len(x)) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and all(isinstance(y, int) for y in x):
#         return _ntt_one(x=x, inv_flag=inv_flag, const_time=const_time)
#     elif isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(inv_flag, bool) and isinstance(const_time, bool) and all(is_pow_two(len(z)) for y in x for z in y):
#         return _ntt_many(x=x, inv_flag=inv_flag, const_time=const_time)
#     raise ValueError(f'Cannot compute _ntt_raw for x, inv_flag, const_time unless x is a list of integers with power-of-two length or a list of lists of lists of integers (with the most internal of these lists having lengths that are powers-of-two), inv_flag and const_time are both bool, and x has power-of-two-length, but had (type(x), inv_flag, const_time)={(type(x), inv_flag, const_time)}.')
#
#
# def ntt(x: PolyCoefs | PolyNTT) -> PolyCoefs | PolyNTT:
#     """
#     Performs the NTT and the inverse NTT, depending on input. If a PolyCoefs object is input, then the NTT of the input
#     data is output. If a PolyNTT object is input, then the inverse NTT of the input data is output.
#
#     See []
#     :params x: Input polynomial matrix representation.
#     :type x: PolyCoefs | PolyNTT
#
#     :return: Transformed data
#     :rtype: PolyCoefs | PolyNTT
#     """
#     if isinstance(x, PolyCoefs):
#         return PolyNTT(vals=_ntt_raw(x.vals, inv_flag=False, const_time=True), modulus=x.modulus, degree=x.degree)
#     elif isinstance(x, PolyNTT):
#         return PolyCoefs(vals=_ntt_raw(x.vals, inv_flag=True, const_time=True), modulus=x.modulus, degree=x.degree)
#     raise ValueError(f'Cannot compute NTT (or inverse NTT) for x unless x is a PolyCoefs or a PolynomialNTTMatrix, but had type(x)={x}.')
#
#
# def transpose(x: PolyCoefs | PolyNTT) -> PolyCoefs | PolyNTT:
#     """
#     Transposes the input data.
#
#     :params x: Input polynomial matrix representation.
#     :type x: PolyCoefs | PolyNTT
#
#     :return: Transposed data
#     :rtype: PolyCoefs | PolyNTT
#     """
#     if isinstance(x, PolyCoefs) or isinstance(x, PolyNTT):
#         result = deepcopy(x)
#         result.vals = [[x.vals[j][i] for j in range(len(x.vals))] for i in range(len(x.vals))]
#         return result
#     raise ValueError(f'Cannot transpose with x unless x is a PolyCoefs or a PolynomialNTTMatrix, but had type(x)={type(x)}.')
#
#
# def xof(x: bytes) -> bytes:
#     """
#     Pass input to whatever XOF your implementation requires.
#
#     :params x: Input data
#     :type x: bytes
#
#     :return: XOF of input data
#     :type x: bytes
#     """
#     if isinstance(x, bytes):
#         return bytes(0)
#     raise ValueError(f'Cannot compute XOF with x unless x is bytes, but had type(x)={type(x)}.')
#
#
# def prf(x: bytes) -> bytes:
#     """
#     Pass input to whatever PRF your implementation requires.
#
#     :params x: Input data
#     :type x: bytes
#
#     :return: PRF of input data
#     :type x: bytes
#     """
#     if isinstance(x, bytes):
#         return bytes(0)
#     raise ValueError(f'Cannot compute PRF with x unless x is bytes, but had type(x)={type(x)}.')
#
#
# def kdf(x: bytes) -> bytes:
#     """
#     Pass input to whatever KDF your implementation requires.
#
#     :params x: Input data
#     :type x: bytes
#
#     :return: KDF of input data
#     :type x: bytes
#     """
#     if isinstance(x, bytes):
#         return bytes(0)
#     raise ValueError(f'Cannot compute KDF with x unless x is bytes, but had type(x)={type(x)}.')
#
#
# def hash_h(x: bytes) -> bytes:
#     """
#     Pass input to whatever hash function H your implementation requires.
#
#     :params x: Input data
#     :type x: bytes
#
#     :return: HashFunctionH of input data
#     :type x: bytes
#     """
#     if isinstance(x, bytes):
#         return bytes(0)
#     raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')
#
#
# def hash_g(x: bytes) -> bytes:
#     """
#     Pass input to whatever hash function G your implementation requires.
#
#     :params x: Input data
#     :type x: bytes
#
#     :return: HashFunctionG of input data
#     :type x: bytes
#     """
#     if isinstance(x, bytes):
#         return bytes(0)
#     raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')
#
#
# def cpa_pke_keygen() -> bytes:
#     """
#     Key generation for the Kyber CPA PKE scheme.
#
#     Works by:
#         1) Sample a random seed d, then hashing d to two new seeds, rho and sigma.
#         2) Feed rho into the XOF to expand out to the matrix A_hat.
#         3) Sample s_hat and e_hat with the CBD_eta function
#         4) Set t_hat = A_hat * s_hat + e_hat
#         5) The secret key is s_hat (encoded) and the public key is rho and t_hat (encoded).
#     """
#     d: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
#     rho: bytes  # for typing
#     sigma: bytes  # for typing
#     rho_and_sigma: bytes = hash_g(d)
#     # assert len(rho_and_sigma) == 2 * SEED_LEN
#     rho: bytes = rho_and_sigma[:len(rho_and_sigma) // 2]
#     sigma: bytes = rho_and_sigma[len(rho_and_sigma) // 2:]
#     index: int = 0
#
#     # Make A_hat
#     a_hat_vals: list[list[list[int]]] = []
#     for i in range(K):
#         a_hat_vals += [[]]
#         for j in range(K):
#             next_xof_input: bytes = rho + int2bytes(x=j) + int2bytes(x=i)  # first j, then i
#             a_hat_vals[-1] += parse(xof(next_xof_input))
#     a_hat: PolyNTT = PolyNTT(vals=a_hat_vals)
#
#     # Make s_hat
#     s_vals: list[list[list[int]]] = []
#     for _ in range(K):
#         s_vals += [[]]
#         next_x: bytes = sigma + int2bytes(x=index)
#         s_vals[-1] += [cbd_eta(x=prf(x=next_x))]
#         index += 1
#     s_hat: PolyNTT = PolyNTT(vals=s_vals)
#
#     # Make e_hat
#     e_vals: list[list[list[int]]] = []
#     for _ in range(K):
#         e_vals += [[]]
#         next_x: bytes = sigma + int2bytes(x=index)
#         e_vals[-1] += [cbd_eta(x=next_x)]
#         index += 1
#     e_hat: PolyNTT = PolyNTT(vals=e_vals)
#
#     # Make t_hat
#     t_hat: PolyNTT = a_hat * s_hat + e_hat
#
#     pk: bytes = encode_m(x=t_hat, m=ENCODED_CPA_PKE_PK_LEN) + rho
#     sk: bytes = encode_m(x=s_hat, m=ENCODED_CPA_PKE_SK_LEN)
#     return pk + sk
#
#
# def _cpa_pke_enc(pk: bytes, plaintext: bytes, randomness: bytes) -> bytes:
#     index: int = 0  # this is N in the specs, but we already used capital N for degree n
#     encoded_t_hat: bytes = pk[:K]
#     rho: bytes = pk[K:]
#     t_hat: PolyNTT = decode_m(x=encoded_t_hat, m=ENCODED_CPA_PKE_PK_LEN)
#     t_hat_transpose: PolyNTT = transpose(x=t_hat)
#
#     a_hat_transpose_vals: list[list[list[int]]] = []
#     for i in range(K):
#         a_hat_transpose_vals += [[]]
#         for j in range(K):
#             next_xof_input: bytes = rho + int2bytes(x=i) + int2bytes(x=i)  # first i, then j
#             a_hat_transpose_vals[-1] += [parse(xof(next_xof_input))]
#     a_hat_transpose: PolyNTT = PolyNTT(vals=a_hat_transpose_vals)
#
#     r_vals: list[list[list[int]]] = []
#     for _ in range(K):
#         r_vals += [[]]
#         next_prf_input: bytes = randomness + int2bytes(x=index)
#         r_vals[-1] += [cbd_eta(x=prf(next_prf_input))]
#         index += 1
#     r: PolyCoefs = PolyCoefs(vals=r_vals)
#
#     e_one_vals: list[list[list[int]]] = []
#     for _ in range(K):
#         e_one_vals += [[]]
#         next_prf_input: bytes = randomness + int2bytes(x=index)
#         e_one_vals[-1] += [cbd_eta(x=prf(next_prf_input))]
#         index += 1
#     e_one: PolyCoefs = PolyCoefs(vals=e_one_vals)
#
#     e_two_vals: list[list[list[int]]] = []
#     next_prf_input: bytes = randomness + int2bytes(x=index)
#     e_two_vals += [[cbd_eta(x=prf(next_prf_input))]]
#     e_two: PolyCoefs = PolyCoefs(vals=e_two_vals)
#
#     r_hat: PolyNTT = ntt(x=r)
#     partial_u_hat: PolyNTT = a_hat_transpose * r_hat
#     u: PolyCoefs = ntt(x=partial_u_hat) + e_one
#     partial_v_hat: PolyNTT = t_hat_transpose * r_hat
#     partial_v: PolyCoefs = ntt(x=partial_v_hat)
#     partial_v += e_two
#     decoded_msg: PolyCoefs = decode_m(x=plaintext, m=1)
#     decompressed_decoded_msg: PolyCoefs = decompress(x=decoded_msg, d=1)
#     v: PolyCoefs = partial_v + decompressed_decoded_msg
#
#     compressed_u: PolyCoefs = compress(x=u, d=D_U)
#     encoded_compressed_u: bytes = encode_m(x=compressed_u, m=D_U)
#
#     compressed_v: PolyCoefs = compress(x=v, d=D_V)
#     encoded_compressed_v: bytes = encode_m(x=compressed_v, m=D_V)
#
#     return encoded_compressed_u + encoded_compressed_v
#
#
# def cpa_pke_encrypt(pk: bytes, plaintext: bytes, r: bytes) -> bytes:
#     """
#     Encryption for the Kyber CPA PKE scheme.
#
#     Works by:
#         1) Parsing the input pk to get t_hat and rho, and transposing t_hat
#         2) Expanding rho to A_hat and then computing the transpose A_hat_transpose
#         3) Sampling a random r and e_one with the CBD_eta function
#         4) Sampling a polynomial e_two with the CBD_eta function.
#         5) Computing u = A_transpose * r + e_1, then compressing and encoding u
#         6) Compute v = t_hat_transpose * r + e_2 + Decompress(Decode(plaintext)), then compressing and encoding v.
#         7) The ciphertext is encode(compress(u)) + encode(compress(v))
#
#     :params pk: Input encoded public key.
#     :type pk: bytes
#     :params plaintext: Input plaintext.
#     :type plaintext: bytes
#     :params: r: Input randomness
#     :type r: bytes
#
#     :return: Ciphertext
#     :rtype: bytes
#     """
#     if isinstance(pk, bytes) and len(pk) >= CPA_PKE_PK_LEN and \
#             isinstance(plaintext, bytes) and len(plaintext) >= SEED_LEN and \
#             isinstance(r, bytes) and len(r) >= SEED_LEN:
#         return _cpa_pke_enc(pk=pk, plaintext=plaintext, randomness=r)
#     raise ValueError(
#         f'Cannot cpa_pke_encrypt pk, plaintext, r unless pk is bytes with len(pk) >= {ENCODED_CPA_PKE_PK_LEN} ' +
#         f'and plaintext is bytes with len(m) >= {SEED_LEN} and r is bytes len(r) >= {SEED_LEN}, but had ' +
#         f'(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))=' +
#         f'{(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))}.')
#
#
# def _cpa_pke_dec(sk: bytes, ciphertext: bytes) -> bytes:
#     encoded_c_one: bytes = ciphertext[:CPA_PKE_FIRST_CIPHERTEXT_LEN]
#     encoded_c_two: bytes = ciphertext[CPA_PKE_FIRST_CIPHERTEXT_LEN:]
#     c_one: PolyCoefs = decode_m(x=encoded_c_one, m=D_U)
#     u: PolyCoefs = decompress(x=c_one, d=D_U)
#     u_hat: PolyNTT = ntt(x=u)
#     c_two: PolyCoefs = decode_m(x=encoded_c_two, m=D_V)
#     v: PolyCoefs = decompress(x=c_two, d=D_V)
#     s_hat: PolyNTT = decode_m(x=sk, m=ENCODED_CPA_PKE_SK_LEN)
#     s_hat_dot_u_hat: PolyNTT = transpose(x=s_hat) * u_hat
#     s_dot_u: PolyCoefs = ntt(x=s_hat_dot_u_hat)
#     compressed_decoded_msg: PolyCoefs = v - s_dot_u
#     decoded_msg: PolyCoefs = decompress(x=compressed_decoded_msg, d=1)
#     msg: bytes = encode_m(x=decoded_msg, m=1)
#     return msg
#
#
# def cpa_pke_decrypt(sk: bytes, ciphertext: bytes) -> bytes:
#     """
#     Decryption for the Kyber CPA PKE scheme.
#
#     Works by:
#         1) Parse the input ciphertext to get the encoded and compressed u and v
#         2) Decode then decompress to get u and decode then decompress to get v
#         3) Parse sk to get s_hat
#         4) Compute v - s_transpose * u as the decoded and decompressed plaintext
#         5) Compress then encode v - s_transpose * u to get the plaintext.
#
#     :params sk: Input secret key
#     :type sk: bytes
#     :params ciphertext: Input ciphertext
#     :type ciphertext: bytes
#
#     :return: Plaintext
#     :rtype: bytes
#     """
#     if isinstance(sk, bytes) and len(sk) >= CPA_PKE_SK_LEN and isinstance(ciphertext, bytes) and len(ciphertext) >= CPA_PKE_CIPHERTEXT_LEN:
#         return _cpa_pke_dec(sk=sk, ciphertext=ciphertext)
#     raise ValueError(f'Cannot cpa_pke_decrypt with sk, ciphertext unless sk is bytes with len(sk) >= ' +
#                      f'{CPA_PKE_SK_LEN} and ciphertext is bytes with len(ciphertext) >= {CPA_PKE_CIPHERTEXT_LEN}' +
#                      f' but had (type(sk), len(sk), type(ciphertext), len(ciphertext))=' +
#                      f'{(type(sk), len(sk), type(ciphertext), len(ciphertext))}.')
#
#
# def cca_kem_keygen() -> bytes:
#     """
#     Key generation for the Kyber CCA KEM scheme.
#
#     Works by:
#         1) Sampling a random seed z
#         2) Running Kyber CPA PKE Key Generation to get a pk and sk_prime, and hashing pk with H to get hashed_pk
#         3) The public key is pk and the secret key is the concatenated sk_prime + pk + hashed_pk + z
#     """
#     z: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
#     cpa_pke_keys: bytes = cpa_pke_keygen()
#     pk: bytes = cpa_pke_keys[:CPA_PKE_PK_LEN]
#     sk_prime: bytes = cpa_pke_keys[CPA_PKE_PK_LEN:]
#     hashed_pk: bytes = hash_h(x=pk)
#     sk: bytes = sk_prime + pk + hashed_pk + z
#     return pk + sk
#
#
# def _cca_kem_enc(pk: bytes) -> bytes:
#     msg: bytes = int2bytes(x=randbits(SEED_LEN*8), length=SEED_LEN)
#     plaintext: bytes = hash_h(x=msg)
#     hashed_pk: bytes = hash_h(x=pk)
#     k_bar_and_r: bytes = hash_g(plaintext + hashed_pk)
#     k_bar: bytes = k_bar_and_r[:SEED_LEN]
#     r: bytes = k_bar_and_r[SEED_LEN:]
#     ciphertext: bytes = cpa_pke_encrypt(pk=pk, plaintext=plaintext, r=r)
#     hashed_c: bytes = hash_h(x=ciphertext)
#     k: bytes = kdf(x=k_bar + hashed_c)
#     return ciphertext + k
#
#
# def cca_kem_encapsulate(pk: bytes) -> bytes:
#     """
#     Encapsulation for the Kyber CCA KEM scheme.
#
#     Works by:
#         1) Sampling a random seed message msg
#         2) Hash the msg with H to conceal system randomness, and hash the public key to get hashed_pk
#         3) Concatenate hashed_msg + hashed_pk, and then hashing the result with G to get k_bar and r
#         4) Compute the Kyber CPA PKE ciphertext of hashed_msg for pk using randomness r
#         5) Hashing the ciphertext to get hashed_ciphertext
#         6) Concatenate k_bar with hashed_ciphertext and then feeding this through KDF to get the key K
#         7) Output the ciphertext and the key K
#
#     :params pk: Input public key.
#     :type pk: bytes
#
#     :return: Encapsulated key and ciphertext
#     :rtype: bytes
#     """
#     if isinstance(pk, bytes) and len(pk) >= CCA_KEM_PK_LEN:
#         return _cca_kem_enc(pk=pk)
#     raise ValueError(f'Cannot cca_kem_encapsulate for pk unless pk is bytes with length at least {CCA_KEM_PK_LEN} but had (type(pk), len(pk))={(type(pk), len(pk))}.')
#
#
# def _cca_kem_dec(ciphertext: bytes, sk: bytes) -> bytes:
#     pk: bytes = sk[CPA_PKE_SK_LEN:CPA_PKE_SK_LEN+CPA_PKE_PK_LEN]
#     hashed_pk: bytes = sk[CPA_PKE_SK_LEN+CPA_PKE_PK_LEN:CPA_PKE_SK_LEN+CPA_PKE_PK_LEN+SEED_LEN]
#     z: bytes = sk[CPA_PKE_SK_LEN+CPA_PKE_PK_LEN+SEED_LEN:]
#     plaintext_prime: bytes = cpa_pke_decrypt(sk=sk, ciphertext=ciphertext)
#     k_bar_prime_and_r_prime: bytes = hash_g(x=plaintext_prime + hashed_pk)
#     k_bar_prime: bytes = k_bar_prime_and_r_prime[:SEED_LEN]
#     r_prime: bytes = k_bar_prime_and_r_prime[SEED_LEN:]
#     ciphertext_prime: bytes = cpa_pke_encrypt(pk=pk, plaintext=plaintext_prime, r=r_prime)
#     hashed_ciphertext: bytes = hash_h(x=ciphertext)
#     if ciphertext_prime == ciphertext:
#         return kdf(x=k_bar_prime + hashed_ciphertext)
#     return kdf(x=z + hashed_ciphertext)
#
#
# def cca_kem_decapsulate(ciphertext: bytes, sk: bytes) -> bytes:
#     """
#     Decapsulation for the Kyber CCA KEM scheme.
#
#     Works by:
#         1) Parse the input to get pk, the hashed_pk, the random z
#         2) Use Kyber CPA PKE Decryption to get a purported plaintext
#         3) Concatenate the purported plaintext and the hashed_pk, then feed this through the hash function G to get k_bar and r
#         4) Re-encrypt to a ciphertext
#         5) If the ciphertexts match, concatenate k_bar and hashed_ciphertext and feed the result through KDF to get K
#         6) Otherwise, concatenate z and hashed_ciphertext and feed the result through KDF to get K.
#
#     :params ciphertext: Input ciphertext
#     :type ciphertext: bytes
#     :params sk: Input secret key sk
#     :type sk: bytes
#
#     :return: Decapsulated key
#     :rtype: bytes
#     """
#     if isinstance(ciphertext, bytes) and len(ciphertext) == CCA_KEM_CIPHERTEXT_LEN and isinstance(sk, bytes) and len(sk) == CCA_KEM_SK_LEN:
#         return _cca_kem_dec(ciphertext=ciphertext, sk=sk)
#     raise ValueError(f'Cannot cca_kem_decapsulate with ciphertext and sk unless ciphertext is bytes with length at least {CCA_KEM_CIPHERTEXT_LEN} and sk is bytes with length at least {CCA_KEM_SK_LEN}, but had (type(ciphertext), len(ciphertext), type(sk), len(sk))={(type(ciphertext), len(ciphertext), type(sk), len(sk))}.')
