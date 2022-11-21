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
from math import ceil, floor, log2, sqrt
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
SEED_LEN_IN_BYTES: int = 32


# Convenient constants derived from the parameters above
LOG_N: int = ceil(log2(N))
HALF_Q: int = Q // 2
LOG_Q: int = ceil(log2(Q))
LOG_Q_BITS_IN_BYTES: int = ceil(LOG_Q / 8)
LOG_K: int = ceil(log2(K))
ZETAS: list[int] = [(ROU ** i) % Q for i in [N // (2 ** (j + 1)) for j in range(LOG_N)]]
ZETA_INVERSES: list[int] = [(ROU_INVERSE ** i) % Q for i in [N // (2 ** (j + 1)) for j in range(LOG_N)]]

CPA_PKE_SK_LEN: int = ceil(K * N * LOG_Q / 8)
CPA_PKE_PK_LEN: int = CPA_PKE_SK_LEN + SEED_LEN_IN_BYTES
CPA_PKE_FIRST_CIPHERTEXT_LEN: int = ceil(K * N * D_U / 8)
CPA_PKE_SECOND_CIPHERTEXT_LEN: int = ceil(N * D_V / 8)
CPA_PKE_CIPHERTEXT_LEN: int = CPA_PKE_FIRST_CIPHERTEXT_LEN + CPA_PKE_SECOND_CIPHERTEXT_LEN
CCA_KEM_SK_LEN: int = CPA_PKE_SK_LEN + CPA_PKE_PK_LEN + 2 * SEED_LEN_IN_BYTES
CCA_KEM_PK_LEN: int = CPA_PKE_PK_LEN
CCA_KEM_CIPHERTEXT_LEN: int = CPA_PKE_CIPHERTEXT_LEN


def _bit_rev(x: int, length_in_bits: int) -> int:
    return int(bin(x)[2:].zfill(length_in_bits)[::-1], 2)


def bit_rev(x: int, length_in_bits: int) -> int:
    """
    Reverse the bits in the binary expansion of x
    :param x: Input integer whose bits are to be reversed
    :type x: int
    :param length_in_bits: Length of binary expansion
    :type length_in_bits: int
    :return: Integer whose binary expansion is the reverse of the input x
    :rtype: int
    """
    if isinstance(length_in_bits, int) and length_in_bits >= 1 and isinstance(x, int) and 0 <= x < 2 ** length_in_bits:
        return _bit_rev(x=x, length_in_bits=length_in_bits)
    elif not isinstance(x, int):
        raise TypeError(f'Cannot bit_rev with x, length unless x is an integer, but had type(x)={type(x)}.')
    elif not isinstance(length_in_bits, int):
        raise TypeError(f'Cannot bit_rev with x, length unless length is an integer, but had type(length)={type(length_in_bits)}.')
    elif length_in_bits < 1:
        raise ValueError(f'Cannot bit_rev with x, length unless length >= 1, but had length={length_in_bits}.')
    raise ValueError(f'Cannot bit_rev with x, length unless 0 <= x < {2 ** length_in_bits} but had x={x}.')


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


def _bit_rev_cp(x: list[int], length_in_bits: int) -> list[int]:
    return [x[bit_rev(x=i, length_in_bits=length_in_bits)] for i in range(len(x))]


def bit_rev_cp(x: list[int], length_in_bits: int) -> list[int]:
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
    :param length_in_bits: Number of bits required to describe the maximum index in x.
    :type x: int
    :return: List of integers with bit-reversed indices.
    :rtype: list[int]
    """
    if isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(length_in_bits, int) and length_in_bits >= 1 and len(x) == 2 ** length_in_bits:
        return _bit_rev_cp(x=x, length_in_bits=length_in_bits)
    elif not isinstance(x, list):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless x is a list, but had type(x)={type(x)}.')
    elif not all(isinstance(y, int) for y in x):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless x is a list of integers, but had (type(y) for y in x)={(type(y) for y in x)}.')
    elif not isinstance(length_in_bits, int):
        raise TypeError(f'Cannot bit_rev_cp with x, num_bits unless num_bits is an integer, but had type(num_bits)={type(length_in_bits)}.')
    elif length_in_bits < 1:
        raise ValueError(f'Cannot bit_rev_cp with x, num_bits unless num_bits >= 1, but had num_bits={length_in_bits}.')
    raise ValueError(f'Cannot bit_rev_cp with x, num_bits unless len(x) == {2 ** length_in_bits} but had len(x)={len(x)}.')


def _reduce(x: int, q: int) -> int:
    y: int = x % q
    z: int = y - q//2 - 1
    w: int = y - (1 + (z >> ceil(log2(q)))) * q
    return w


def reduce(x: int, q: int = Q) -> int:
    """
    Compute some integer w such that -Q//2 <= w <= Q//2 and such that (x-w) % Q == 0, in constant
    time.
    :param x: Input integer
    :type x: int
    :param q: Input modulus
    :type q: int
    :return: "Centered" representative of x % Q.
    :rtype: int
    """
    if isinstance(x, int) and isinstance(q, int) and q >= 2:
        return _reduce(x=x, q=q)
    elif not isinstance(x, int):
        raise TypeError(f'Cannot compute reduce for x unless x is an integer, but had type(x)={type(x)}.')
    elif not isinstance(q, int):
        raise TypeError(f'Cannot compute reduce with q unless q is an integer, but had type(q)={type(q)}.')
    elif q < 2:
        raise ValueError(f'Cannot compute reduce with q unless q >= 2, but had q={q}.')


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


def _our_parse_one(x: bytes) -> tuple[list[int], int]:
    """
    Parse an input bytes object into a 2-tuple, where the first entry is a list of N=256 integers in the list
    [0, 1, 2, 3, ..., Q - 1] for Q=3329, and the second entry is an integer describing the index of the
    first unused byte in x.
    TODO: Extend this to work for arbitrary N and arbitrary Q.
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


def _our_parse_many(x: bytes) -> list[list[list[int]]]:
    """
    Parse an input bytes object into a list of K=4 lists of K=4 lists of N=256 integers, each in the list
    [0, 1, 2, 3, ..., Q-1] for Q=3329, by calling _parse_one several times.
    TODO: Extend this to work for arbitrary K, arbitrary N, and arbitrary Q.
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
            next_result, starting_index = _our_parse_one(x=x[starting_index:])
            result[-1] += [next_result]
    return result


def parse(x: bytes) -> list[int] | list[list[list[int]]]:
    """
    Parse an input bytes object into either a list of 256 integers, or a list of K=4 lists of K=4 lists of N=256
    integers, each in the list [0, 1, 2, 3, ..., Q-1] for Q=3329.

    Our implementation does NOT parse according to specifications, because we only input a list of bytes, not a
    bytestream. Therefore, our implementation is not compatible with the NIST standard. The standard has no chance of
    failing, but our implementation has a small positive probability of failing. However, if parsing does not
    fail, then the resulting 'A' matrix should match the standard.

    :param x: Input bytes
    :type x: bytes
    :return: A list of lists of lists of integers.
    :rtype: list[list[list[int]]]
    """
    if isinstance(x, bytes) and len(x) == N*(LOG_Q + SECBITS)//8:
        return _our_parse_one(x=x)[0]
    elif isinstance(x, bytes) and len(x) >= K*K*N*(LOG_Q + SECBITS)//8:
        return _our_parse_many(x=x)
    elif not isinstance(x, bytes):
        raise TypeError(f'Cannot parse with x unless x is a bytes object, but had type(x)={type(x)}.')
    raise ValueError(f'Cannot parse with x unless x is a bytes object with length {N * (LOG_Q + SECBITS) // 8} or length at least {K * K * N * (LOG_Q + SECBITS) // 8} but had len(x)={len(x)}.')


def is_arithmetic_legal(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]], q_a: int, q_b: int) -> bool:
    """
    Return true if both inputs are consistently sized (all rows have the same number of columns and all columns have the
    same number of degrees) and both inputs have the same number of rows, columns, and degrees.

    :param x: Input bytes
    :type x: bytes
    :return: A list of lists of lists of integers.
    :rtype: list[list[list[int]]]
    """
    if q_a == q_b >= 1:
        num_rows_in_self: int = len(a_vals)
        num_rows_in_other: int = len(b_vals)

        min_cols_in_self: int = min(len(x) for x in a_vals)
        min_cols_in_other: int = min(len(x) for x in b_vals)
        max_cols_in_self: int = max(len(x) for x in a_vals)
        max_cols_in_other: int = max(len(x) for x in b_vals)
        consistent_cols_in_self: bool = max_cols_in_self == min_cols_in_self
        consistent_cols_in_other: bool = max_cols_in_other == min_cols_in_other

        min_deg_in_self: int = min(len(x) for i in a_vals for x in i)
        min_deg_in_other: int = min(len(x) for i in b_vals for x in i)
        max_deg_in_self: int = max(len(x) for i in a_vals for x in i)
        max_deg_in_other: int = max(len(x) for i in b_vals for x in i)
        consistent_deg_in_self: bool = max_deg_in_self == min_deg_in_self
        consistent_deg_in_other: bool = max_deg_in_other == min_deg_in_other

        same_rows: bool = num_rows_in_self == num_rows_in_other
        same_cols: bool = max_cols_in_self == max_cols_in_other
        same_deg: bool = max_deg_in_self == max_deg_in_other

        return same_rows and consistent_cols_in_self and consistent_cols_in_other and same_cols and consistent_deg_in_self and consistent_deg_in_other and same_deg
    return False


def add(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]], q: int) -> list[list[list[int]]]:
    """
    Input two lists of lists of lists of integers such that is_arithmetic_legal returns true, and simply adds their
    entries coordinate-wise before reducing in constant time.

    :param a_vals: Input first matrix
    :type a_vals: list[list[list[int]]]
    :param b_vals: Input second matrix
    :type b_vals: list[list[list[int]]]
    :return: Coordinate-wise sum, reduced
    :rtype: list[list[list[int]]]
    """
    result = deepcopy(a_vals)
    for i, row in enumerate(b_vals):
        for j, col in enumerate(row):
            for k, x in enumerate(col):
                tmp = result[i][j][k] + x
                result[i][j][k] = reduce(x=tmp, q=q)
    return result


def mul(a_vals: list[list[list[int]]], b_vals: list[list[list[int]]], q: int) -> list[list[list[int]]]:
    """
    Input two lists of lists of lists of integers such that is_arithmetic_legal returns true, and simply multiplies
    their entries coordinate-wise before reducing in constant time.

    :param a_vals: Input first matrix
    :type a_vals: list[list[list[int]]]
    :param b_vals: Input second matrix
    :type b_vals: list[list[list[int]]]
    :return: Coordinate-wise sum, reduced
    :rtype: list[list[list[int]]]
    """
    result = deepcopy(a_vals)
    for i, row in enumerate(b_vals):
        for j, col in enumerate(row):
            for k, x in enumerate(col):
                result[i][j][k] = reduce(x=result[i][j][k] * x, q=q)
    return result


class PolyCoefs(object):
    """
    Class for coefficient representations of matrices of polynomials. Supports addition but not multiplication.

    TODO: Removed checks for q-1 % 2n == 0

    Attributes
    ----------
        q: int
            Integer modulus for all coefficients (0 <= coefficient < modulus)
        n: int
            Integer degree for all polynomials
        k1: int
            Number of rows
        k2: int
            Number of cols
        half_q: int
            Half of modulus, rounded down.
        log_q: int
            Number of bits to describe a coefficient.
        zetas: list[int]
            Powers of the root of unity
        zeta_inverses: list[int]
            Powers of the inverse of the root of unity
        vals: list[list[list[int]]]
            2-dimensional matrix with list[int] entries, each entry is a coefficient representation of a polynomial.
        const_time_flag: bool
            Flag for describing whether arithmetic is constant-time - in prod, always keep at True.

    Methods
    -------
        __init__(self, vals, modulus, degree)
    """
    q: int = Q
    n: int = N
    k1: int = K
    k2: int = 1
    half_q: int = HALF_Q
    log_q: int = LOG_Q
    zetas: list[int] = []
    zeta_inverses: list[int] = []
    vals: list[list[list[int]]] = []
    const_time_flag: bool = True

    def __init__(self, vals: list[list[list[int]]] | None, q: int = Q, n: int = N, k1: int = K, k2: int = 1, zetas: list[int] = ZETAS, zeta_inverses: list[int] = ZETA_INVERSES, const_time_flag: bool = True):
        if isinstance(q, int) and q >= 2 and isinstance(n, int) and n >= 1 and is_pow_two(x=n) and isinstance(k1, int) and k1 >= 1 and isinstance(k2, int) and k2 >= 1 and isinstance(const_time_flag, bool) and isinstance(vals, list) and all(isinstance(x, list) for x in vals) and all(isinstance(y, list) for x in vals for y in x) and all(isinstance(z, int) for x in vals for y in x for z in y) and len(vals) == k1 and all(len(x)==k2 for x in vals) and all(len(y) == n for x in vals for y in x) and (all(0 <= z < q for x in vals for y in x for z in y) or all(-q//2 <= z <= q//2 for x in vals for y in x for z in y)):
            self.q = q
            self.n = n
            self.k1 = k1
            self.k2 = k2
            self.half_q = q // 2
            self.log_q = ceil(log2(q))
            self.zetas, self.zeta_inverses = make_zetas_and_invs(q=q, n=n, lgn=ceil(log2(n)))
            self.vals = vals
            self.const_time_flag = const_time_flag
        elif not isinstance(q, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer modulus q, but had type(q)={type(q)}.')
        elif q < 2:
            raise ValueError(f'Cannot instantiate PolyCoefs object with modulus q < 2 but had q={q}.')
        elif not isinstance(n, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer degree n, but had type(n)={type(n)}.')
        elif n < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with degree n < 1 but had n={n}.')
        elif not is_pow_two(x=n):
            raise ValueError(f'Cannot instantiate PolyCoefs object with non-power-of-two degree n, but had n={n}.')
        elif not isinstance(k1, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer number of rows k1, but had type(k1)={type(k1)}.')
        elif k1 < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with number of rows k1 < 1, but had k1={k1}.')
        elif not isinstance(k2, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer number of columns k2, but had type(k2)={type(k2)}.')
        elif k2 < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with number of columns k2 < 1, but had k2={k2}.')
        elif not isinstance(const_time_flag, bool):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-boolean const_time_flag, but had type(const_time_flag)={const_time_flag}.')
        elif not isinstance(vals, list):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-list vals, but had type(vals)={type(vals)}.')
        elif len(vals) != k1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals without k1 rows, but had (len(vals), k1)={(len(vals), k1)}. Check that you specify k1 for non-standard sized data.')
        elif not all(isinstance(x, list) for x in vals):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every row in vals is a list.')
        elif not all(len(x) == k2 for x in vals):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row in vals has k2={k2} columns. Check that you specify k1 for non-standard sized data.')
        elif not all(isinstance(y, list) for x in vals for y in x):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every column in every row is a list.')
        elif not all(len(y) == n for x in vals for y in x):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row and every column is a list with degree n={n}.')
        elif not all(isinstance(z, int) for x in vals for y in x for z in y):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every row, column, and degree is an integer.')
        elif not all(0 <= z < q for x in vals for y in x for z in y) and not all(-q//2 <= z < q//2 for x in vals for y in x for z in y):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row, column, and degree is an integer modulo q.')
        else:
            raise RuntimeError(f'Some unspecified error occurred while instantiating a PolyCoefs object. Please contact the developers with a description of the input to the __init__ function.')

    def __eq__(self, other):
        return self.q == other.q and \
               self.k1 == other.k1 and \
               self.k2 == other.k2 and \
               all((a-b) % self.q == 0 for x, y in zip(self.vals, other.vals) for z, w in zip(x, y) for a, b in zip(z, w))

    def __add__(self, other):
        if is_arithmetic_legal(a_vals=self.vals, b_vals=other.vals, q_a=self.q, q_b=other.q):
            result = deepcopy(self)
            result.vals = add(a_vals=self.vals, b_vals=other.vals, q=self.q)
            return result
        raise ValueError(f'Cannot compute PolyCoefs.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [[[-coef for coef in col] for col in row] for row in other.vals]
        return self.__add__(other=other)

    def __mul__(self, other):
        raise NotImplementedError(f'Multiplication of PolyCoefs is not implemented. Apply the NTT and multiply the resulting PolyNTT objects instead.')

    def __mod__(self, other):
        # Computes modulus in NON-CONSTANT TIME. Only use this if the input data is public information.
        result = deepcopy(self)
        result.vals = [[[z % other for z in y] for y in x] for x in result.vals]
        return result


class PolyNTT(object):
    """
    Class for coefficient representations of matrices of polynomials. This class does not support arithmetic.

    Attributes
    ----------
        q: int
            Integer modulus for all coefficients (0 <= coefficient < modulus)
        n: int
            Integer degree for all polynomials
        k1: int
            Number of rows
        k2: int
            Number of cols
        half_q: int
            Half of modulus, rounded down.
        log_q: int
            Number of bits to describe a coefficient.
        zetas: list[int]
            Powers of the root of unity
        zeta_inverses: list[int]
            Powers of the inverse of the root of unity
        vals: list[list[list[int]]]
            2-dimensional matrix with list[int] entries, each entry is a coefficient representation of a polynomial.
        const_time: bool
            Flag for describing whether arithmetic is constant-time - in prod, always keep at True.

    Methods
    -------
        __init__(self, vals, modulus, degree, const_time)
            Set self.modulus to modulus, self.degree to degree, compute self.halfmod and self.logmod, set self.vals to vals, and set self.const_time to const_time
        __add__(self, other)
            Checks if dimensions of self and other are compatible, and then return the sum of self and other.
        __radd__(self, other)
            If other == 0, return self, other call __add__. This way, we can use sum().
        __sub__(self, other)
            Deepcopy other, negate each entry in vals, and then call __add__.
        __mul__(self, other)
            Check if self is 1-by-1 or if self and other have dimensions match. Then, if self is 1-by-1, call _scalar_mul, otherwise call _matrix_mul.
        _scalar_mul(self, other)
            Multiply each polynomial in other by self.
        _matrix_mul(self, other)
            Multiply the matrices self and other coordinate-wise.
    """
    q: int = Q
    n: int = N
    k1: int = K
    k2: int = 1
    half_q: int = HALF_Q
    log_q: int = LOG_Q
    zetas: list[int] = []
    zeta_inverses: list[int] = []
    vals: list[list[list[int]]] = []
    const_time_flag: bool = True

    def __init__(self, vals: list[list[list[int]]] | None, q: int = Q, n: int = N, k1: int = K, k2: int = 1, zetas: list[int]=ZETAS, zeta_inverses: list[int]=ZETA_INVERSES, const_time_flag: bool = True):
        if isinstance(q, int) and q >= 2 and isinstance(n, int) and n >= 1 and is_pow_two(x=n) and isinstance(k1, int) and k1 >= 1 and isinstance(k2, int) and k2 >= 1 and isinstance(const_time_flag, bool) and isinstance(vals, list) and all(isinstance(x, list) for x in vals) and all(isinstance(y, list) for x in vals for y in x) and all(isinstance(z, int) for x in vals for y in x for z in y) and len(vals) == k1 and all(len(x)==k2 for x in vals) and all(len(y) == n for x in vals for y in x) and (all(0 <= z < q for x in vals for y in x for z in y) or all(-q//2 <= z <= q//2 for x in vals for y in x for z in y)):
            self.q = q
            self.n = n
            self.k1 = k1
            self.k2 = k2
            self.half_q = q // 2
            self.log_q = ceil(log2(q))
            self.zetas, self.zeta_inverses = make_zetas_and_invs(q=q, n=n, lgn=ceil(log2(n)))
            self.vals = vals
            self.const_time_flag = const_time_flag
        elif not isinstance(q, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer modulus q, but had type(q)={type(q)}.')
        elif q < 2:
            raise ValueError(f'Cannot instantiate PolyCoefs object with modulus q < 2 but had q={q}.')
        elif not isinstance(n, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer degree n, but had type(n)={type(n)}.')
        elif n < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with degree n < 1 but had n={n}.')
        elif not is_pow_two(x=n):
            raise ValueError(f'Cannot instantiate PolyCoefs object with non-power-of-two degree n, but had n={n}.')
        elif not isinstance(k1, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer number of rows k1, but had type(k1)={type(k1)}.')
        elif k1 < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with number of rows k1 < 1, but had k1={k1}.')
        elif not isinstance(k2, int):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-integer number of columns k2, but had type(k2)={type(k2)}.')
        elif k2 < 1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with number of columns k2 < 1, but had k2={k2}.')
        elif not isinstance(const_time_flag, bool):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-boolean const_time_flag, but had type(const_time_flag)={const_time_flag}.')
        elif not isinstance(vals, list):
            raise TypeError(f'Cannot instantiate PolyCoefs object with non-list vals, but had type(vals)={type(vals)}.')
        elif len(vals) != k1:
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals without k1 rows, but had (len(vals), k1)={(len(vals), k1)}.')
        elif not all(isinstance(x, list) for x in vals):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every row in vals is a list.')
        elif not all(len(x) == k2 for x in vals):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row in vals has k2={k2} columns.')
        elif not all(isinstance(y, list) for x in vals for y in x):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every column in every row is a list.')
        elif not all(len(y) == n for x in vals for y in x):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row and every column is a list with degree n={n}.')
        elif not all(isinstance(z, int) for x in vals for y in x for z in y):
            raise TypeError(f'Cannot instantiate PolyCoefs object with vals unless every row, column, and degree is an integer.')
        elif not all(0 <= z < q for x in vals for y in x for z in y) and not all(-q//2 <= z < q//2 for x in vals for y in x for z in y):
            raise ValueError(f'Cannot instantiate PolyCoefs object with vals unless every row, column, and degree is an integer modulo q.')
        else:
            raise RuntimeError(f'Some unspecified error occurred while instantiating a PolyCoefs object. Please contact the developers with a description of the input to the __init__ function.')

    def __eq__(self, other):
        return self.q == other.q and \
               self.k1 == other.k1 and \
               self.k2 == other.k2 and \
               all(0 == (a-b) % self.q for x, y in zip(self.vals, other.vals) for z, w in zip(x, y) for a, b in zip(z, w))

    def __add__(self, other):
        if is_arithmetic_legal(a_vals=self.vals, b_vals=other.vals, q_a=self.q, q_b=other.q):
            result = deepcopy(self)
            result.vals = add(a_vals=self.vals, b_vals=other.vals, q=self.q)
            return result
        raise ValueError(f'Cannot compute PolyNTT.__add__ unless dimensions of both matrices match (dim mismatch). Check if number of rows, number of columns, and degrees all match.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [[[-coef for coef in col] for col in row] for row in other.vals]
        return self.__add__(other=other)

    def __mul__(self, other):
        if isinstance(other, PolyNTT) and isinstance(other, PolyNTT):
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
                'Cannot compute PolynomialNTTMatrix.__mul__ unless self is 1x1 and both have consistent degrees, or where self is mxn, other is nxp, and both have consistent degrees (dim mismatch).')
        raise NotImplementedError(f'Cannot compute PolyNTT.__mul__ unless both self and other are PolyNTT.')

    def _scalar_mul(self, other):
        result = deepcopy(other)
        for i in range(len(other.vals)):
            for j in range(len(other.vals[0])):
                for k in range(len(self.vals[0][0])):
                    result.vals[i][j][k] = (self.vals[0][0][k] * other.vals[i][j][k]) % self.q
        return result

    def _matrix_mul(self, other):
        result = deepcopy(self)
        result.k1 = self.k1
        result.k2 = other.k2
        result.vals = [[[0 for k in range(self.n)] for j in range(other.k2)] for i in range(self.k1)]
        for i in range(self.k1):
            for j in range(other.k2):
                tmp: list[int] = [0 for _ in range(self.n)]
                for l in range(self.k2):
                    for k in range(self.n):
                        tmp[k] += self.vals[i][l][k] * other.vals[l][j][k]
                result.vals[i][j] = [coef % self.q for coef in tmp]  # [reduce(x=sum(self.vals[i][l][k] * other.vals[l][j][k] for l in range(self.k2)), q=self.q) for k in range(self.n)]
        return result

    def __mod__(self, other):
        # Computes modulus in NON-CONSTANT TIME. Only use this if the input data is public information.
        result = deepcopy(self)
        result.vals = [[[z % other for z in y] for y in x] for x in result.vals]
        return result


def _cbd_eta(x: bytes, eta: int = ETA) -> list[int]:
    x_as_bits: str = ''
    for next_bytes in x:
        next_bytes_as_bits = bin(next_bytes)[2:]
        x_as_bits += bin(next_bytes)[2:].zfill(8*ceil(len(next_bytes_as_bits)/8))
    result: list[int] = [0] * N
    for i in range(N):
        a: int = sum([int('1' == x_as_bits[2*i*eta + j]) for j in range(eta)])
        b: int = sum([int('1' == x_as_bits[2*i*eta + eta + j]) for j in range(eta)])
        result[i] = a - b
    return result


def cbd_eta(x: bytes, eta: int = ETA) -> list[int]:
    """
    Sample a list of integers with length N from input bytes x, such that the integers are sampled from the centered
    binomial distribution (CBD) with support [-eta, -eta + 1, ..., eta - 1, eta]. Works by looking at x as a bitstring,
    and taking the difference of sums of sections of the bit string.

    :param x: Input bytes
    :type x: bytes
    :param eta: Input bound.
    :type eta: int

    :return: List of integers.
    :rtype: list[int]
    """
    if isinstance(x, bytes) and len(x) >= 2*N*eta//8 and isinstance(eta, int) and eta >= 1:
        return _cbd_eta(x=x, eta=eta)
    elif not isinstance(x, bytes):
        raise TypeError(f'Cannot cbd_eta with x, eta unless x is a bytes object, but type(x)={type(x)}.')
    elif len(x) < 2*N*eta//8:
        raise ValueError(f'Cannot cbd_eta with x, eta unless len(x)>={2*N*eta//8} but had len(x)={len(x)}.')
    elif not isinstance(eta, int):
        raise TypeError(f'Cannot cbd_eta with x, eta unless eta is an integer, but had type(eta)={type(eta)}.')
    elif eta < 1:
        raise ValueError(f'Cannot cbd_eta with x, eta unless eta >= 1 but had eta={eta}.')
    raise TypeError(f'Cannot cbd_eta with x, eta unless eta is an integer, but had type(x)={type(x)}.')


def _cbd_polycoefs(x: bytes, eta: int = ETA, num_rows: int = K, num_cols: int = 1) -> PolyCoefs:
    vals: list[list[list[int]]] = []
    for i in range(num_rows):
        next_x: bytes = x[i*num_cols*(2*N*eta//8): (i+1)*num_cols*(2*N*eta)//8]
        vals += [[]]
        for j in range(num_cols):
            next_bytes: bytes = bytes(next_x[j * (2 * N * eta)//8: (j + 1) * (2 * N * eta)//8])
            vals[-1] += [cbd_eta(x=next_bytes, eta=eta)]
    return PolyCoefs(vals=vals, k1=num_rows, k2=num_cols)


def cbd_polycoefs(x: bytes, eta: int = ETA, num_rows: int = K, num_cols: int = 1) -> PolyCoefs:
    """
    Sample a matrix of polynomials with num_rows rows, num_cols columns, degrees N, and such that each coefficient is in
    the list [-eta, -eta + 1, ..., eta - 1, eta]. Works by calling cbd_eta for each polynomial in the matrix.

    :param x: Input bytes
    :type x: bytes
    :param eta: Integer bound eta
    :type eta: int
    :param num_rows: Number of rows.
    :type num_rows: int
    :param num_cols: Number of columns.
    :type num_cols: int

    :return: PolyCoefs object
    :rtype: PolyCoefs
    """
    if isinstance(x, bytes) and len(x) >= num_rows*num_cols*N*2*eta//8 and isinstance(eta, int) and eta >= 1 and isinstance(num_rows, int) and num_rows >= 1 and isinstance(num_cols, int) and num_cols >= 1:
        return _cbd_polycoefs(x=x, eta=eta, num_rows=num_rows, num_cols=num_cols)
    elif not isinstance(x, bytes):
        raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless x is a bytes object, but had type(x)={type(x)}.')
    elif len(x) < num_rows*num_cols*N*2*eta//8:
        raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless len(x) >= {num_rows*num_cols*N*2*eta//8} but had len(x)={len(x)}.')
    elif not isinstance(eta, int):
        raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless eta is an integer, but had type(eta)={type(eta)}.')
    elif eta < 1:
        raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless eta >= 1, but had eta={eta}.')
    elif not isinstance(num_rows, int):
        raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_rows is an integer, but had type(num_rows)={type(num_rows)}.')
    elif num_rows < 1:
        raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_rows >= 1, but had num_rows={num_rows}.')
    elif not isinstance(num_cols, int):
        raise TypeError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_cols is an integer, but had type(num_cols)={type(num_cols)}.')
    raise ValueError(f'Cannot cbd_PolyCoefs with x, eta, num_rows, num_cols unless num_cols >= 1, but had num_cols={num_cols}.')


def _compress_one_int(x: int, d: int, p: int = Q) -> int:
    # We rename q to p to avoid ambiguity with the constant Q
    return round_up(x=x * 2 ** d / p) % 2 ** d


def _compress_list_of_ints(x: list[int], d: int, p: int = Q) -> list[int]:
    return [_compress_one_int(x=y, d=d, p=p) for y in x]


def _compress_many_ints(x: list[list[list[int]]], d: int, p: int = Q) -> list[list[list[int]]]:
    return [[_compress_list_of_ints(x=z, d=d, p=p) for z in y] for y in x]


def _compress_polycoefs(x: PolyCoefs, d: int, p: int = Q) -> PolyCoefs:
    return PolyCoefs(vals=_compress_many_ints(x=x.vals, d=d, p=p), q=x.q, n=x.n, k1=x.k1, k2=x.k2, const_time_flag=x.const_time_flag)


def _should_compress_many(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> bool:
    return isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2


def compress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
    """
    Compresses an integer modulo p (or a list of them, or a list of lists of lists of them) or the integers modulo p in
    a PolyCoefs object to a d-bit integer as specified. Works by, essentially, looking at x/p as a ratio, and going that
    far along the list [0, 1, 2, ..., 2**d - 1].

    :param x: Input data
    :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
    :param d: Input number of bits for compressed data
    :type d: int
    :param p: Input modulus
    :type p: int

    :return: Compressed data
    :rtype: int | list[int] | list[list[list[int]]] | PolyCoefs

    """
    if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _compress_one_int(x=x, d=d, p=p)
    elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _compress_list_of_ints(x=x, d=d, p=p)
    elif _should_compress_many(x=x, d=d, p=p):
        return _compress_many_ints(x=x, d=d, p=p)
    elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _compress_polycoefs(x=x, d=d, p=p)
    raise ValueError(f'Cannot compute compress for x, d, p unless x is an integer, a list of integers, or a list of lists of lists of integers, or a PolyCoefs... and d is an integer >= 1 and p is an integer >= 2, but had (type(x), d, p)={(type(x), d, p)}.')


def _decompress_one_int(x: int, d: int, p: int = Q) -> int:
    unrounded_result = x * p / 2 ** d
    rounded_result = round_up(x=unrounded_result) % p
    return rounded_result


def _decompress_list_of_ints(x: list[int], d: int, p: int = Q) -> list[int]:
    decompressed_ints: list[int] = [_decompress_one_int(x=y, d=d, p=p) for y in x]
    return decompressed_ints


def _decompress_many_ints(x: list[list[list[int]]], d: int, p: int = Q) -> list[list[list[int]]]:
    return [[_decompress_list_of_ints(x=z, d=d, p=p) for z in y] for y in x]


def _decompress_polycoefs(x: PolyCoefs, d: int, p: int = Q) -> PolyCoefs:
    decompressed_vals: list[list[list[int]]] = _decompress_many_ints(x=x.vals, d=d, p=p)
    return PolyCoefs(vals=decompressed_vals, q=x.q, n=x.n, k1=x.k1, k2=x.k2)


def decompress(x: int | list[int] | list[list[list[int]]] | PolyCoefs, d: int, p: int = Q) -> int | list[int] | list[list[list[int]]] | PolyCoefs:
    """
    Decompresses an integer modulo 2**d (or a list of them, or a list of lists of lists of them) or the integers modulo
    2**d in a PolyCoefs object to an integer modulo p as specified. Works by, essentially, looking at x/2**d as a ratio,
    and going that far along the list [0, 1, 2, ..., p-1].

    :param x: Input data
    :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
    :param d: Input number of bits for compressed data
    :type d: int
    :param p: Input modulus
    :type p: int

    :return: Decompressed data
    :rtype: int | list[int] | list[list[list[int]]] | PolyCoefs
    """
    if isinstance(x, int) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _decompress_one_int(x=x, d=d, p=p)
    elif isinstance(x, list) and all(isinstance(y, int) for y in x) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _decompress_list_of_ints(x=x, d=d, p=p)
    elif _should_compress_many(x=x, d=d, p=p):
        return _decompress_many_ints(x=x, d=d, p=p)
    elif isinstance(x, PolyCoefs) and isinstance(d, int) and d >= 1 and isinstance(p, int) and p >= 2:
        return _decompress_polycoefs(x=x, d=d, p=p)
    raise ValueError(f'Cannot decompress with x, d, p unless x is an integer, or a list of integers, or a list of lists of lists of integers, or a PolyCoefs... and d is an integer >= 1, and p is an integer >= 2, but had (type(x), d, p)={(type(x), d, p)}.')


def _encode_m_list_of_ints_to_bytes(x: list[int], bits_per_int: int) -> bytes:
    resulting_ints: list[int] = []
    resulting_bitstring: str = ''
    for y in x:
        y_in_bin: str = bin(y)[2:].zfill(bits_per_int)  # big endian
        nib_ni_y: str = y_in_bin[::-1]  # lil endian
        resulting_bitstring += nib_ni_y
    while len(resulting_bitstring) % 8 != 0:
        resulting_bitstring += '0'  # pad zeros
    next_byte: str
    while len(resulting_bitstring) > 0:
        next_byte, resulting_bitstring = resulting_bitstring[:8], resulting_bitstring[8:]
        resulting_ints += [int(next_byte, 2)]
    return bytes(resulting_ints)


def _encode_m_many(x: list[list[list[int]]], bits_per_int: int) -> bytes:
    result = bytes(0)
    for y in x:
        for z in y:
            result += _encode_m_list_of_ints_to_bytes(x=z, bits_per_int=bits_per_int)
    return result


def _encode_m_matrix(x: PolyCoefs | PolyNTT, bits_per_int: int) -> tuple[bytes, int, int, int, int, bool]:
    flag: bool = isinstance(x, PolyCoefs) and not isinstance(x, PolyNTT)
    return (_encode_m_many(x=x.vals, bits_per_int=bits_per_int), x.q, x.n, x.k1, x.k2, flag)


def encode_m(x: int | list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT, bits_per_int: int) -> bytes | tuple[bytes, int, int, int, int, bool]:
    """
    We rename encode_l to encode_m to use m instead of l because l is an ambiguous character. Encodes an integer (or a
    list of integers, or a list of lists of lists of integers) or the integers in a PolyCoefs or PolyNTT object to
    bytes. Works with the usual serialization... each piece of input data is an integer modulo p, and we just pad out
    the binary expansion of the input data to m bytes.

    :param x: Input data
    :type x: int | list[int] | list[list[list[int]]] | PolyCoefs
    :param bits_per_int: Number of bytes
    :type bits_per_int: int

    :return: Encoded data
    :rtype: bytes
    """
    if isinstance(bits_per_int, int) and bits_per_int >= 1 and isinstance(x, list) and all(isinstance(y, int) for y in x) and all(0 <= y < 2 ** bits_per_int for y in x):
        return _encode_m_list_of_ints_to_bytes(x=x, bits_per_int=bits_per_int)
    elif isinstance(bits_per_int, int) and bits_per_int >= 1 and isinstance(x, list) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and all(isinstance(w, int) for y in x for z in y for w in z) and all(0 <= w < 2 ** bits_per_int for y in x for z in y for w in z):
        return _encode_m_many(x=x, bits_per_int=bits_per_int)
    elif isinstance(bits_per_int, int) and bits_per_int >= 1 and not isinstance(x, list) and (isinstance(x, PolyCoefs) ^ isinstance(x, PolyNTT)):
        return _encode_m_matrix(x=x, bits_per_int=bits_per_int)
    elif not isinstance(bits_per_int, int):
        raise TypeError(f'Cannot compute encode with bits_per_int unless bits_per_int is an integer, but had type(bits_per_int)={type(bits_per_int)}.')
    elif bits_per_int < 1:
        raise ValueError(f'Cannot compute decode with negative or zero bits_per_int, but had bits_per_int={bits_per_int}.')
    elif not isinstance(x, list) and not (isinstance(x, PolyCoefs) ^ isinstance(x, PolyNTT)):
        raise TypeError(f'Cannot compute decode with x unless x is a list, a PolyCoefs object, or a PolyNTT object, but had type(x)={type(x)}.')
    elif isinstance(x, list) and not isinstance(x, PolyCoefs) and not isinstance(x, PolyNTT) and all(isinstance(y, list) for y in x) and all(isinstance(z, list) for y in x for z in y) and not all(isinstance(w, int) for y in x for z in y for w in z):
        raise TypeError(f'Cannot compute decode with x when x is not PolyCoefs, not PolyNTT, and not list[int], but instead is a list[list[list[_]]] for some _ other than integers, but had (type(w) for y in x for z in y for w in z)={(type(w) for y in x for z in y for w in z)}.')
    elif isinstance(x, list) and not isinstance(x, PolyCoefs) and not isinstance(x, PolyNTT) and not all(isinstance(y, list) for y in x) and not all(isinstance(y, int) for y in x):
        raise TypeError(f'Cannot compute decode with x when x is not PolyCoefs, not PolyNTT, and not list[list[list[int]]] and is a list[_] for some _ other than integers, but had (type(y) for y in x)={(type(y) for y in x)}.')
    raise ValueError(f'Cannot compute encode with x unless x is a list[int] or a list[list[list[int]]] where all the integers are appropriately sized.')


def _decode_m_list_of_ints_from_bytes(x: bytes, bits_per_int: int) -> list[int]:
    resulting_ints: list[int] = []
    resulting_bitstring: str = ''
    tni_txen: str
    next_int_in_bin: str
    for y in x:
        resulting_bitstring += bin(y)[2:].zfill(8)
    while len(resulting_bitstring) > 0:
        nib_ni_tni_txen, resulting_bitstring = resulting_bitstring[:bits_per_int], resulting_bitstring[bits_per_int:]
        next_int_in_bin: str = nib_ni_tni_txen[::-1]
        if len(next_int_in_bin) == bits_per_int:
            resulting_ints += [int(next_int_in_bin, 2)]
    return resulting_ints


def _decode_m_many(x: bytes, bits_per_int: int, k1: int = K, k2: int = 1, n: int = N) -> list[list[list[int]]]:
    result: list[list[list[int]]] = []
    remaining_bytes: bytes = deepcopy(x)
    next_row: bytes
    next_col: bytes
    bytes_per_row: int = len(x) // k1
    bytes_per_col: int = bytes_per_row // k2
    while len(remaining_bytes) > 0:
        result += [[]]
        next_row, remaining_bytes = remaining_bytes[:bytes_per_row], remaining_bytes[bytes_per_row:]
        while len(next_row) > 0:
            next_col, next_row = next_row[:bytes_per_col], next_row[bytes_per_col:]
            next_entry: list[int] = _decode_m_list_of_ints_from_bytes(x=next_col, bits_per_int=bits_per_int)
            result[-1] += [next_entry]
    return result


def _decode_m_matrix(x: bytes, bits_per_int: int, ntt_matrix_flag: bool, q: int = Q, k1: int = K, k2: int = 1, n: int = N) -> PolyCoefs | PolyNTT:
    if ntt_matrix_flag:
        return PolyNTT(vals=_decode_m_many(x=x, bits_per_int=bits_per_int, k1=k1, k2=k2, n=n), k1=k1, k2=k2, n=n, q=q)
    return PolyCoefs(vals=_decode_m_many(x=x, bits_per_int=bits_per_int, k1=k1, k2=k2, n=n), k1=k1, k2=k2, n=n, q=q)


def decode_m(x: bytes, bits_per_int: int, coef_matrix_flag: bool = False, ntt_matrix_flag: bool = False, k1: int = K, k2: int = 1, n: int = N) -> list[int] | list[list[list[int]]] | PolyCoefs | PolyNTT:
    """
    We rename decode_l to decode_m to use m instead of l because l is an ambiguous character. Decodes a bytes object to
    an integer (or a list of integers, or a list of lists of lists of integers) or a PolyCoefs or PolyNTT object. Works
    with the usual serialization... interpret each chunk of m bytes as determined by the binary expansion of an integer,
    written as bytes.

    :param x: Input data
    :type x: bytes
    :param bits_per_int: Number of bytes per integer
    :type bits_per_int: int
    :param coef_matrix_flag: Flag indicating whether input x should be decoded to a PolyCoefs object.
    :type coef_matrix_flag: bool
    :param ntt_matrix_flag: Flag indicating whether input x should be decoded to a PolyNTT object.
    :type ntt_matrix_flag: bool
    :param k1: Integer indicating number of rows of the decoded matrix.
    :type k1: int
    :param k2: Integer indicating number of columsn of the decoded matrix.
    :type k2: int

    :return: Encoded data
    :rtype: bytes
    """
    if isinstance(x, bytes) and len(x) == ceil(n*bits_per_int/8) and isinstance(coef_matrix_flag, bool) and not coef_matrix_flag and isinstance(ntt_matrix_flag, bool) and not ntt_matrix_flag:
        # In this case, we do not care about the coef_matrix_flag, the ntt_matrix_flag, k1, or k2.
        return _decode_m_list_of_ints_from_bytes(x=x, bits_per_int=bits_per_int)
    elif isinstance(x, bytes) and len(x) >= ceil(n*bits_per_int*k1*k2/8) and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and not coef_matrix_flag and not ntt_matrix_flag:
        # In this case, we want a list[list[list[int]]] object returned, not a PolyNTT or a PolyCoefs
        return _decode_m_many(x=x, bits_per_int=bits_per_int)
    elif isinstance(x, bytes) and len(x) >= ceil(n*bits_per_int*k1*k2/8) and len(x) % (k1*k2) == 0 and isinstance(coef_matrix_flag, bool) and isinstance(ntt_matrix_flag, bool) and (coef_matrix_flag ^ ntt_matrix_flag):
        # In this case, coef_matrix_flag XOR ntt_matrix_flag is True, so exactly one of these is true, and we want that sort of object in return.
        return _decode_m_matrix(x=x, bits_per_int=bits_per_int, ntt_matrix_flag=ntt_matrix_flag, k1=k1, k2=k2, n=n)
    elif not isinstance(x, bytes):
        raise TypeError(f'Cannot compute decode_m unless input x is a bytes object, but had type(x)={type(x)}.')
    elif len(x) != ceil(n*bits_per_int/8) and len(x) < ceil(n*bits_per_int*k1*k2/8):
        raise ValueError(f'Cannot compute decode_m unless input x is a bytes object with length exactly {ceil(n*bits_per_int/8)} or at least {ceil(n*bits_per_int*k1*k2/8)} but had len(x) = {len(x)}.')
    elif len(x) != ceil(n*bits_per_int/8) and len(x) % (k1*k2) != 0:
        raise ValueError(f'Cannot compute decode_m with input x with len(x) >= {ceil(n*k1*k2*bits_per_int/8)} unless len(x) % {k1*k2} == 0, but had len(x) % {k1*k2} = {len(x) % (k1*k2)}.')
    elif not isinstance(coef_matrix_flag, bool) or not isinstance(ntt_matrix_flag, bool):
        raise TypeError(f'Cannot compute decode_m unless coef_matrix_flag and ntt_matrix_flag are both bool, but had type(coef_matrix_flag), type(ntt_matrix_flag) = {(type(coef_matrix_flag), type(ntt_matrix_flag))}.')
    raise ValueError(f'Cannot compute decode_m if both coef_matrix_flag and ntt_matrix_flag are True.')


def _ntt_one(x: list[int], inv_flag: bool, const_time: bool = True, q: int = Q, n: int = N, log_n: int = LOG_N, half_q: int = HALF_Q, zetas: list[int] = ZETAS, zeta_inverses: list[int] = ZETA_INVERSES) -> list[int]:
    bit_rev_x: list[int] = bit_rev_cp(x=x, length_in_bits=ceil(log2(len(x))))
    m: int = 1
    for s in range(1, log_n + 1):
        m *= 2
        if inv_flag:
            this_zeta: int = zeta_inverses[s - 1]
        else:
            this_zeta: int = zetas[s - 1]
        for k in range(0, n, m):
            w: int = 1
            for j in range(m // 2):
                t: int = w * bit_rev_x[k + j + m // 2]
                u: int = bit_rev_x[k + j]
                if const_time:
                    bit_rev_x[k + j]: int = reduce(x=u + t, q=q)
                    bit_rev_x[k + j + m // 2]: int = reduce(x=u - t, q=q)
                else:
                    bit_rev_x[k + j]: int = (u + t) % q
                    if bit_rev_x[k + j] > half_q:
                        bit_rev_x[k + j] = bit_rev_x[k + j] - q
                    bit_rev_x[k + j + m // 2]: int = (u - t) % q
                    if bit_rev_x[k + j + m // 2] > half_q:
                        bit_rev_x[k + j + m // 2] = bit_rev_x[k + j + m // 2] - q
                w *= this_zeta
    if inv_flag:
        n_inv: int = 1
        while (n_inv * n) % q != 1:
            n_inv += 1
        if const_time:
            bit_rev_x: list[int] = [reduce(x=n_inv * i, q=q) for i in bit_rev_x]
        else:
            bit_rev_x: list[int] = [(n_inv * i) % q for i in bit_rev_x]
            bit_rev_x = [i if i <= half_q else i - q for i in bit_rev_x]
    elif const_time:
        bit_rev_x: list[int] = [reduce(x=i, q=q) for i in bit_rev_x]
    else:
        bit_rev_x: list[int] = [i % q for i in bit_rev_x]
        bit_rev_x = [i if i <= half_q else i - q for i in bit_rev_x]
    return bit_rev_x


def _ntt_many(x: list[list[list[int]]], inv_flag: bool, const_time: bool = True, q: int = Q, n: int = N, log_n: int = LOG_N, half_q: int = HALF_Q, zetas: list[int] = ZETAS, zeta_inverses: list[int] = ZETA_INVERSES) -> list[list[list[int]]]:
    return [[_ntt_one(x=z, inv_flag=inv_flag, const_time=const_time, q=q, n=n, log_n=log_n, half_q=half_q, zetas=zetas, zeta_inverses=zeta_inverses) for z in y] for y in x]


def is_prime(x: int) -> bool:
    return all([x % i != 0 for i in range(2, ceil(sqrt(x)) + 1)])


def has_prim_rou(q: int, n: int) -> bool:
    return q % n == 1


def is_ntt_friendly_prime(q: int, n: int) -> bool:
    return is_prime(q) and is_pow_two(n) and has_prim_rou(q=q, n=n)


def get_prim_rou_and_rou_inv(q: int, n: int) -> None | tuple[int, int]:
    if not (is_ntt_friendly_prime(q, n)):
        raise ValueError('Input q and d are not ntt-friendly prime and degree.')
    # If we do not raise a ValueError, then there exists a primitive root of unity 2 <= x < q.
    x: int = 2
    while x < q:
        if all(x ** k % q != 1 for k in range(1, n)) and x ** n % q == 1:
            break
        x += 1
    return x, ((x ** (n - 1)) % q)


def make_zetas_and_invs(q: int, n: int, lgn: int) -> tuple[list[int], list[int]]:
    powers: list[int] = [n // (2 ** (s + 1)) for s in range(lgn)]
    zeta, zeta_inv = get_prim_rou_and_rou_inv(q=q, n=n)
    left: list[int] = [int(zeta ** i) % q for i in powers]
    left = [i if i <= q // 2 else i - q for i in left]
    right: list[int] = [int(zeta_inv ** i) % q for i in powers]
    right = [i if i <= q // 2 else i - q for i in right]
    return left, right


def ntt(x: PolyCoefs | PolyNTT, q: int = Q, n: int = N, log_n: int = LOG_N, half_q: int = HALF_Q) -> PolyCoefs | PolyNTT:
    """
    Performs the NTT and the inverse NTT, depending on input. If a PolyCoefs object is input, then the NTT of the input
    data is output. If a PolyNTT object is input, then the inverse NTT of the input data is output.

    See []
    :params x: Input polynomial matrix representation.
    :type x: PolyCoefs | PolyNTT

    :return: Transformed data
    :rtype: PolyCoefs | PolyNTT
    """
    if isinstance(x, PolyCoefs):
        vals=_ntt_many(x.vals, inv_flag=False, const_time=x.const_time_flag, q=q, n=n, log_n=log_n, half_q=half_q, zetas=x.zetas, zeta_inverses=x.zeta_inverses)
        return PolyNTT(vals=vals, q=x.q, n=x.n, k1=x.k1, k2=x.k2, const_time_flag=x.const_time_flag)
    elif isinstance(x, PolyNTT):
        return PolyCoefs(vals=_ntt_many(x.vals, inv_flag=True, const_time=x.const_time_flag, q=q, n=n, log_n=log_n, half_q=half_q, zetas=x.zetas, zeta_inverses=x.zeta_inverses), q=x.q, n=x.n, k1=x.k1, k2=x.k2, const_time_flag=x.const_time_flag)
    raise ValueError(f'Cannot compute NTT (or inverse NTT) for x unless x is a PolyCoefs or a PolynomialNTTMatrix, but had type(x)={x}.')


def transpose(x: PolyCoefs | PolyNTT) -> PolyCoefs | PolyNTT:
    """
    Transposes the input data.

    :params x: Input polynomial matrix representation.
    :type x: PolyCoefs | PolyNTT

    :return: Transposed data
    :rtype: PolyCoefs | PolyNTT
    """
    if isinstance(x, PolyCoefs) or isinstance(x, PolyNTT):
        new_vals: list[list[list[int]]] = [[[] for i in range(x.k1)] for j in range(x.k2)]
        for i in range(x.k1):
            for j in range(x.k2):
                new_vals[j][i] = x.vals[i][j]
        if isinstance(x, PolyCoefs):
            return PolyCoefs(vals=new_vals, k1=x.k2, k2=x.k1)
        return PolyNTT(vals=new_vals, k1=x.k2, k2=x.k1)
    raise ValueError(f'Cannot transpose with x unless x is a PolyCoefs or a PolynomialNTTMatrix, but had type(x)={type(x)}.')


def xof(x: bytes) -> bytes:
    """
    Pass input to whatever XOF your implementation requires.

    :params x: Input data
    :type x: bytes

    :return: XOF of input data
    :type x: bytes
    """
    if isinstance(x, bytes):
        return bytes(0)
    raise ValueError(f'Cannot compute XOF with x unless x is bytes, but had type(x)={type(x)}.')


def prf(x: bytes) -> bytes:
    """
    Pass input to whatever PRF your implementation requires.

    :params x: Input data
    :type x: bytes

    :return: PRF of input data
    :type x: bytes
    """
    if isinstance(x, bytes):
        return bytes(0)
    raise ValueError(f'Cannot compute PRF with x unless x is bytes, but had type(x)={type(x)}.')


def kdf(x: bytes) -> bytes:
    """
    Pass input to whatever KDF your implementation requires.

    :params x: Input data
    :type x: bytes

    :return: KDF of input data
    :type x: bytes
    """
    if isinstance(x, bytes):
        return bytes(0)
    raise ValueError(f'Cannot compute KDF with x unless x is bytes, but had type(x)={type(x)}.')


def hash_h(x: bytes) -> bytes:
    """
    Pass input to whatever hash function H your implementation requires.

    :params x: Input data
    :type x: bytes

    :return: HashFunctionH of input data
    :type x: bytes
    """
    if isinstance(x, bytes):
        return bytes(0)
    raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')


def hash_g(x: bytes) -> bytes:
    """
    Pass input to whatever hash function G your implementation requires.

    :params x: Input data
    :type x: bytes

    :return: HashFunctionG of input data
    :type x: bytes
    """
    if isinstance(x, bytes):
        return bytes(0)
    raise ValueError(f'Cannot compute HashFunctionH with x unless x is bytes, but had type(x)={type(x)}.')


def cpa_pke_keygen() -> bytes:
    """
    Key generation for the Kyber CPA PKE scheme.

    Works by:
        1) Sample a random seed d, then hashing d to two new seeds, rho and sigma.
        2) Feed rho into the XOF to expand out to the matrix A_hat.
        3) Sample s_hat and e_hat with the CBD_eta function
        4) Set t_hat = A_hat * s_hat + e_hat
        5) The secret key is s_hat (encoded) and the public key is rho and t_hat (encoded).
    """
    d: bytes = (randbits(SEED_LEN_IN_BYTES * 8)).to_bytes(length=SEED_LEN_IN_BYTES, byteorder='big')
    rho: bytes  # for typing
    sigma: bytes  # for typing
    rho_and_sigma: bytes = hash_g(d)
    # assert len(rho_and_sigma) == 2 * SEED_LEN
    rho: bytes = rho_and_sigma[:SEED_LEN_IN_BYTES]
    sigma: bytes = rho_and_sigma[SEED_LEN_IN_BYTES:]
    index: int = 0

    # Make A_hat
    a_hat_vals: list[list[list[int]]] = []
    for i in range(K):
        a_hat_vals += [[]]
        for j in range(K):
            next_xof_input: bytes = rho + j.to_bytes(length=1, byteorder='big') + i.to_bytes(length=1, byteorder='big')
            a_hat_vals[-1] += [parse(xof(next_xof_input))]
    a_hat: PolyNTT = PolyNTT(vals=a_hat_vals, k1=K, k2=K)

    # Make s_hat
    s_vals: list[list[list[int]]] = []
    for _ in range(K):
        s_vals += [[]]
        next_x: bytes = sigma + index.to_bytes(length=1, byteorder='big')
        s_vals[-1] += [cbd_eta(x=prf(x=next_x))]
        index += 1
    s_hat: PolyNTT = PolyNTT(vals=s_vals)

    # Make e_hat
    e_vals: list[list[list[int]]] = []
    for _ in range(K):
        e_vals += [[]]
        next_x: bytes = sigma + index.to_bytes(length=1, byteorder='big')
        e_vals[-1] += [cbd_eta(x=next_x)]
        index += 1
    e_hat: PolyNTT = PolyNTT(vals=e_vals)

    # Make t_hat
    t_hat: PolyNTT = (a_hat * s_hat + e_hat) % Q

    encoded_t_hat: tuple[bytes, int, int, int, int, bool] = encode_m(x=t_hat, bits_per_int=LOG_Q)
    pk: bytes = encoded_t_hat[0] + rho

    encoded_s_hat: tuple[bytes, int, int, int, int, bool] = encode_m(x=s_hat, bits_per_int=LOG_Q)
    sk: bytes = encoded_s_hat[0]
    return pk + sk


def _cpa_pke_enc(pk: bytes, plaintext: bytes, randomness: bytes) -> bytes:
    index: int = 0  # this is N in the specs, but we already used capital N for degree n
    encoded_t_hat: bytes = pk[:-SEED_LEN_IN_BYTES]  # split encoded_t_hat from the pk
    rho: bytes = pk[-SEED_LEN_IN_BYTES:]  # split seed from the pk
    t_hat = decode_m(x=encoded_t_hat, bits_per_int=LOG_Q, ntt_matrix_flag=True)
    t_hat_transpose: PolyNTT = transpose(x=t_hat)

    a_hat_transpose_vals: list[list[list[int]]] = []
    for i in range(K):
        a_hat_transpose_vals += [[]]
        for j in range(K):
            next_xof_input: bytes = rho + i.to_bytes(length=1, byteorder='big') + j.to_bytes(length=1, byteorder='big')
            a_hat_transpose_vals[-1] += [parse(xof(next_xof_input))]
    a_hat_transpose: PolyNTT = PolyNTT(vals=a_hat_transpose_vals, k1=K, k2=K)

    r_vals: list[list[list[int]]] = []
    for _ in range(K):
        r_vals += [[]]
        next_prf_input: bytes = randomness + index.to_bytes(length=1, byteorder='big')
        r_vals[-1] += [cbd_eta(x=prf(next_prf_input))]
        index += 1
    r: PolyCoefs = PolyCoefs(vals=r_vals)

    e_one_vals: list[list[list[int]]] = []
    for _ in range(K):
        e_one_vals += [[]]
        next_prf_input: bytes = randomness + index.to_bytes(length=1, byteorder='big')
        e_one_vals[-1] += [cbd_eta(x=prf(next_prf_input))]
        index += 1
    e_one: PolyCoefs = PolyCoefs(vals=e_one_vals)

    e_two_vals: list[list[list[int]]] = []
    next_prf_input: bytes = randomness + index.to_bytes(length=1, byteorder='big')
    e_two_vals += [[cbd_eta(x=prf(next_prf_input))]]
    e_two: PolyCoefs = PolyCoefs(vals=e_two_vals, k1=1, k2=1)

    r_hat: PolyNTT = ntt(x=r)
    partial_u_hat: PolyNTT = a_hat_transpose * r_hat
    u: PolyCoefs = ntt(x=partial_u_hat) + e_one
    partial_v_hat: PolyNTT = t_hat_transpose * r_hat
    partial_v: PolyCoefs = ntt(x=partial_v_hat)
    partial_v += e_two
    decoded_msg: PolyCoefs = decode_m(x=plaintext, bits_per_int=1, k1=1, k2=1, n=N, coef_matrix_flag=True)
    decompressed_decoded_msg: PolyCoefs = decompress(x=decoded_msg, d=1)
    v: PolyCoefs = partial_v + decompressed_decoded_msg

    compressed_u: PolyCoefs = compress(x=u, d=D_U)
    encoded_compressed_u_full: tuple[bytes, int, int, int, int, bool] = encode_m(x=compressed_u, bits_per_int=D_U)
    encoded_compressed_u: bytes = encoded_compressed_u_full[0]

    compressed_v: PolyCoefs = compress(x=v, d=D_V)
    encoded_compressed_v_full: tuple[bytes, int, int, int, int, bool] = encode_m(x=compressed_v, bits_per_int=D_V)
    encoded_compressed_v: bytes = encoded_compressed_v_full[0]

    return encoded_compressed_u + encoded_compressed_v


def cpa_pke_encrypt(pk: bytes, plaintext: bytes, r: bytes) -> bytes:
    """
    Encryption for the Kyber CPA PKE scheme.

    Works by:
        1) Parsing the input pk to get t_hat and rho, and transposing t_hat
        2) Expanding rho to A_hat and then computing the transpose A_hat_transpose
        3) Sampling a random r and e_one with the CBD_eta function
        4) Sampling a polynomial e_two with the CBD_eta function.
        5) Computing u = A_transpose * r + e_1, then compressing and encoding u
        6) Compute v = t_hat_transpose * r + e_2 + Decompress(Decode(plaintext)), then compressing and encoding v.
        7) The ciphertext is encode(compress(u)) + encode(compress(v))

    :params pk: Input encoded public key.
    :type pk: bytes
    :params plaintext: Input plaintext.
    :type plaintext: bytes
    :params: r: Input randomness
    :type r: bytes

    :return: Ciphertext
    :rtype: bytes
    """
    if isinstance(pk, bytes) and len(pk) >= CPA_PKE_PK_LEN and \
            isinstance(plaintext, bytes) and len(plaintext) >= SEED_LEN and \
            isinstance(r, bytes) and len(r) >= SEED_LEN:
        return _cpa_pke_enc(pk=pk, plaintext=plaintext, randomness=r)
    raise ValueError(
        f'Cannot cpa_pke_encrypt pk, plaintext, r unless pk is bytes with len(pk) >= {CPA_PKE_PK_LEN} ' +
        f'and plaintext is bytes with len(m) >= {SEED_LEN_IN_BYTES} and r is bytes len(r) >= {SEED_LEN_IN_BYTES}, but had ' +
        f'(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))=' +
        f'{(type(pk), len(pk), type(plaintext), len(plaintext), type(r), len(r))}.')


def _cpa_pke_dec(sk: bytes, ciphertext: bytes) -> bytes:
    encoded_c_one: bytes = ciphertext[:CPA_PKE_FIRST_CIPHERTEXT_LEN]  # first CPA_PKE_FIRST_CIPHERTEXT_LEN bytes
    encoded_c_two: bytes = ciphertext[CPA_PKE_FIRST_CIPHERTEXT_LEN:]  # next CPA_PKE_SECOND_CIPHERTEXT_LEN bytes
    c_one: PolyCoefs = decode_m(x=encoded_c_one, bits_per_int=D_U, coef_matrix_flag=True)
    u: PolyCoefs = decompress(x=c_one, d=D_U)
    u_hat: PolyNTT = ntt(x=u)
    c_two: PolyCoefs = decode_m(x=encoded_c_two, bits_per_int=D_V, coef_matrix_flag=True, k1=1, k2=1)
    v: PolyCoefs = decompress(x=c_two, d=D_V)
    s_hat: PolyNTT = decode_m(x=sk, bits_per_int=LOG_Q, ntt_matrix_flag=True)

    decoded_msg: PolyCoefs = compress(x=v - ntt(x=transpose(x=s_hat) * u_hat), d=1)
    msg_full: tuple[bytes, int, int, int, int, bool] = encode_m(x=decoded_msg, bits_per_int=1)
    msg: bytes = msg_full[0]
    return msg


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
