"""
A python implementation of a generalized NTT that can handle the NTT approach from Kyber and Dilithium equally well, and handles polynomial arithmetic from the ring $R = \mathbb{Z}[X]/(coef_mod, X^deg_mod + const)$ in constant time. More general NTT approaches that use roots of unity of even smaller degree than that of Kyber should work, too, but at the cost of slower arithmetic operations, as expected from non-standard NTT approaches.

TODO: This code will be moved to lattice-algebra.

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
from math import ceil, sqrt, log2
from copy import deepcopy


def _is_odd_prime(val: int) -> bool:
    """
    Sieve of Eratosthenes
    :param val: Input integer being tested for odd primality.
    :type val: int
    :return: Boolean indicating whether val is an odd prime.
    :rtype: bool
    """
    return all([val % i != 0 for i in range(2, ceil(sqrt(val)) + 1)])


HALF_INTEGERS: dict[int] = dict()
LG_INTEGERS: dict[int] = dict()


def _reduce(val: int, mod: int) -> int:
    """
    Montgomery reduction for constant time centralized modular arithmetic.
    :param val: Input value being reduced
    :type val: int
    :param mod: Input modulus
    :type mod: int
    :return: Reduced val between -(mod//2) and mod//2, inclusive.
    :rtype: int
    """
    if mod not in HALF_INTEGERS:
        HALF_INTEGERS[mod] = mod // 2
    half_mod = HALF_INTEGERS[mod]
    if mod not in LG_INTEGERS:
        LG_INTEGERS[mod] = ceil(log2(mod))
    lg_mod = LG_INTEGERS[mod]
    return (val % mod) - (1 + (((val % mod) - half_mod - 1) >> lg_mod)) * mod


def _is_int_and_gt_one_and_pow_two(val: int) -> bool:
    """
    Test whether val is an integer of the form 2**i for i >= 1.
    :param val: Input value being tested
    :type val: int
    :return: Boolean indicating whether val = 2**i for some i >= 1.
    :rtype: int
    """
    return isinstance(val, int) and \
        val > 1 and \
        not (val & (val - 1))


class Coefs(object):
    """
    Class for coefficient representations of polynomials in the ring $R = Z[X]/(coef_mod, X^deg_mod + const)$ where $Z$ denotes the integers, coef_mod is a prime, deg_mod is a power-of-two, and (coef_mod - 1) % (2**i * deg_mod) == 0 for some integer i. In particular, coefficients are integers modulo q, and if arithmetic causes degree to meet or exceed $deg_mod$, a wrap-around occurs, sending $X^deg_mod$ to $-const$. Arithmetic takes place in constant time.

    Multiplication of Coefs objects is not advisable, because it is slow. Instead, compute the NTT of a Coefs object to get a PtVals object, then do arithmetic with the output of that, and then invert the NTT to get back to a coefficient representation.

    Attributes
    ----------
        coef_mod: int
            Integer modulus for coefficients (0 <= coefficient < modulus or -(coef_mod//2) <= coefficient <= coef_mod)
        deg_mod: int
            Integer degree bound
        const: int
            Integer constant for wrap-around so that X^deg_mod + const = 0.

    Methods
    -------
        __init__(self, coef_mod, deg_mod, const, vals)
            Sets self.coef_mod to coef_mod, self.deg_mod to deg_mod, self.const to const, and then applies Montgomery reduction to entries in val before setting self.vals to the result.
        __repr__()
            Outputs "Coefs(coef_mod, deg_mod, const, vals)"
        __eq__(self, other)
            Compares coef_mod, deg_mod, and const for self and other, and checks if diff'd entries in val are == 0 mod q
        __add__(self, other)
            Compares coef_mod, deg_mod, and const, and returns the component-wise sum of vals modulo q
        __radd__(self, other)
        __neg__(self)
            Deepcopies self, then negates each value in val before returning the result.
        __sub__(self, other)
            Deepcopies other, negates it, and then returns the sum
        __mul__(self, other)
            Multiplication by foiling.
        __rmul__
    """
    coef_mod: int
    deg_mod: int
    const: int
    ord_of_prim_rou: int
    vals: list[int]

    def __init__(self, coef_mod: int, deg_mod: int, const: int, ord_of_prim_rou: int, vals: list[int]):
        if _is_odd_prime(val=coef_mod) and \
                _is_int_and_gt_one_and_pow_two(val=deg_mod) and \
                _mod_has_nth_prim_rou(mod=coef_mod, n=ord_of_prim_rou):
            self.coef_mod = coef_mod
            self.deg_mod: int = deg_mod
            self.const = const
            self.ord_of_prim_rou = ord_of_prim_rou

            # To set vals, we first mod out by X^deg_mod + const then mod out by coef_mod.
            # First mod out by X^deg_mod + const
            these_vals: list[int] = deepcopy(vals)
            right: list[int]
            left: list[int]
            max_len_right_left: int
            while len(these_vals) > deg_mod:
                right = [-const*x for x in vals[deg_mod:]]
                left = vals[:deg_mod]
                max_len_right_left = max(len(right), len(left))
                if len(right) < max_len_right_left:
                    right += [0 for _ in range(max_len_right_left - len(right))]
                if len(left) < max_len_right_left:
                    left += [0 for _ in range(max_len_right_left - len(left))]
                these_vals = [x + y for x, y in zip(left, right)]
            if len(these_vals) < deg_mod:
                these_vals += [0 for _ in range(deg_mod - len(these_vals))]
            # Now mod out by coef_mod
            self.vals = [_reduce(val=x, mod=coef_mod) for x in these_vals]
        elif not _is_odd_prime(val=coef_mod):
            raise ValueError(f'Cannot instantiate a Coefs object without a coef_mod that is a positive odd prime.')
        elif not _is_int_and_gt_one_and_pow_two(val=deg_mod):
            raise ValueError(f'Cannot instantiate Coefs object without a deg_mod that is a power-of-two.')
        else:
            raise ValueError(f'Cannot instantiate a Coefs object unless coef_mod has a primitive n-th root of unity for some power-of-two n <= 2*deg_mod.')

    def __repr__(self) -> str:
        return f'Coefs{(self.coef_mod, self.deg_mod, self.const, self.vals)}'

    def __eq__(self, other) -> bool:
        return isinstance(other, Coefs) and \
            self.coef_mod == other.coef_mod and \
            self.deg_mod == other.deg_mod and \
            self.const == other.const and \
            all((x-y) % self.coef_mod == 0 for x, y in zip(self.vals, other.vals))

    def __add__(self, other):
        if isinstance(other, Coefs) and self.coef_mod == other.coef_mod and self.deg_mod == other.deg_mod and self.const == other.const:
            result = deepcopy(self)
            result.vals = [_reduce(val=x + y, mod=self.coef_mod) for x, y in zip(self.vals, other.vals)]
            return result
        raise ValueError

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        result = deepcopy(self)
        result.vals = [-_ for _ in result.vals]
        return result

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other = -negative_other
        return self + negative_other

    def __mul__(self, other):
        if isinstance(other, Coefs) and self.coef_mod == other.coef_mod and self.deg_mod == other.deg_mod and self.const == other.const:
            result = deepcopy(self)
            result.vals = [0 for _ in range(2 * len(self.vals))]
            for i, a_coef in enumerate(self.vals):
                for j, b_coef in enumerate(other.vals):
                    result.vals[i+j] += a_coef * b_coef
            lower_vals = result.vals[:self.deg_mod]
            upper_vals = result.vals[self.deg_mod:]
            result.vals = [_reduce(val=x - y * self.const, mod=result.coef_mod) for x, y in zip(lower_vals, upper_vals)]
            return result
        raise ValueError

    def __rmul__(self, other):
        return self.__mul__(other=other)


def _mod_has_nth_prim_rou(mod: int, n: int) -> bool:
    """
    Returns a boolean indicating whether a primitive n-th root of unity modulo mod exists.
    :param mod: Input modulus
    :type mod: int
    :param n: Input order of primitive root of unity
    :type n: int
    :return: Boolean indicating whether a primitive n-th root of unity modulo mod exists
    :rtype: bool
    """
    return isinstance(mod, int) and mod > 1 and isinstance(n, int) and n > 1 and (mod - 1) % n == 0


PRIM_NTH_ROUS_AND_ROU_INVS: dict[tuple[int, int]] = dict()


def _get_nth_prim_rou_and_inv(mod: int, n: int) -> tuple[int, int]:
    """
    If a primitive n-th root of unity (rou) exists modulo mod, return rou and the multiplicative inverse (rou_inv). Otherwise raise a runtime error... should this be a value error or something?
    TODO: nth primitive roots of unity ain't unique and we should prefer powers of 2 rous for fast arithmetic
    :param mod: Input modulus
    :type mod: int
    :param n: Input degree of primitive root of unity
    :type n: int
    :return: Tuple (rou, rou_inv) where rou is a primitive n-th root of unity modulo mod, and (rou_inv * rou) % mod == 1
    :rtype: tuple[int, int]
    """
    if _mod_has_nth_prim_rou(mod=mod, n=n) and (mod, n) not in PRIM_NTH_ROUS_AND_ROU_INVS:
        rou: int = 2
        while rou < mod:
            if _reduce(val=rou ** n, mod=mod) == 1 and all(_reduce(val=rou ** k, mod=mod) != 1 for k in range(1, n)):
                break
            rou += 1
        rou_inv = _reduce(val=rou ** (n - 1), mod=mod)
        PRIM_NTH_ROUS_AND_ROU_INVS[(mod, n)] = rou, rou_inv
    elif not _mod_has_nth_prim_rou(mod=mod, n=n):
        raise RuntimeError(f'Cannot compute {n}-th primitive root of unity modulo {mod} unless (mod - 1) % n == 0, but had {(mod - 1) % n}.')
    return PRIM_NTH_ROUS_AND_ROU_INVS[(mod, n)]


class SplitCoefs(object):
    """
    Class for representing lists of coefficient representations of polynomials in a split ring, so that arithmetic takes place coordinate-wise in the list.

    In the case that this list is all polynomials with degree-bound 1 (i.e. constants), this can be interpreted as a point-value representation of a polynomial.

    Attributes
    ----------
        vals: list[Coefs]
            A list of Coefs objects.

    Methods
    -------
        __init__(self, coef_mod, deg_mod, const, vals)
            Sets self.coef_mod to coef_mod, self.deg_mod to deg_mod, self.const to const, and then applies Montgomery reduction to entries in val before setting self.vals to the result.
        __repr__()
            Outputs "SplitCoefs(entry for entry in vals)"
        __eq__(self, other)
            Checks if other is also a SplitCoefs object with the same length vals attribute and entries in val match.
        __add__(self, other)
            Checks if other is also a SplitCoefs object with the same length and each entry in vals share the same coef_mod, deg_mod, and const before returning the component-wise sum.
        __radd__(self, other)
        __neg__(self)
            Deepcopies self, then negates each value in val before returning the result.
        __sub__(self, other)
            Checks if other is also a SplitCoefs object with the same length and each entry in vals share the same coef_mod, deg_mod, and const before negating other and then returning the sum
        __mul__(self, other)
            Checks if other is also a SplitCoefs object with the same length and each entry in vals share the same coef_mod, deg_mod, and const before returning the component-wise product.
        __rmul__
    """
    vals: list[Coefs]

    def __init__(self, vals: list[Coefs]):
        self.vals = vals

    def __repr__(self) -> str:
        return f'SplitCoefs{tuple(self.vals)}'

    def __eq__(self, other) -> bool:
        return isinstance(other, SplitCoefs) and len(self.vals) == len(other.vals) and all(x == y for x, y in zip(self.vals, other.vals))

    def __add__(self, other):
        if isinstance(other, SplitCoefs) and len(self.vals) == len(other.vals) and all(x.coef_mod == y.coef_mod and x.deg_mod == y.deg_mod and x.const == y.const for x, y in zip(self.vals, other.vals)):
            result = deepcopy(self)
            result.vals = [x + y for x, y in zip(self.vals, other.vals)]
            return result
        raise ValueError

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        result = deepcopy(self)
        for i, val in enumerate(result.vals):
            result[i] = -val
        return result

    def __sub__(self, other):
        if isinstance(other, SplitCoefs) and len(self.vals) == len(other.vals) and all(x.coef_mod == y.coef_mod and x.deg_mod == y.deg_mod and x.const == y.const for x, y in zip(self.vals, other.vals)):
            negative_other = deepcopy(other)
            negative_other = -negative_other
            return self + negative_other
        raise ValueError

    def __mul__(self, other):
        if isinstance(other, SplitCoefs) and len(self.vals) == len(other.vals) and all(x.coef_mod == y.coef_mod and x.deg_mod == y.deg_mod and x.const == y.const for x, y in zip(self.vals, other.vals)):
            result = deepcopy(self)
            result.vals = [x * y for x, y in zip(self.vals, other.vals)]
            return result
        raise ValueError

    def __rmul__(self, other):
        return self.__mul__(other=other)


ALREADY_COMPUTED_ZETAS_AND_INVS: dict[tuple[int, int]] = dict()


def _make_zetas_and_invs(mod: int, n: int) -> tuple[list[int], list[int]]:
    """
    Compute the primitive n-th root of unity modulo mod and its multiplicative inverse, and their powers modulo mod in the order used in the NTT algorithm.
    :param mod: Input modulus
    :type mod: int
    :param n: Input degree of primitive root of unity
    :type n: int
    :return: A tuple containing lists of integers. The first list consists of the powers of the primitive n-th root of unity, in the order they are used in the NTT, and the second list consists of the powers of the multiplicative inverse of the root of unity, in the order they are used in the INTT.
    :rtype: tuple[list[int], list[int]]
    """
    if (mod, n) not in ALREADY_COMPUTED_ZETAS_AND_INVS:
        if n not in LG_INTEGERS:
            LG_INTEGERS[n] = ceil(log2(n))
        lg_rou_deg = LG_INTEGERS[n]
        powers: list[int] = [n // (2 ** (s + 1)) for s in range(lg_rou_deg)]
        zeta, zeta_inv = _get_nth_prim_rou_and_inv(mod=mod, n=n)
        rou_powers: list[int] = [int(zeta ** i) % mod for i in powers]
        rou_powers = [i if i <= mod // 2 else i - mod for i in rou_powers]
        rou_inverse_powers: list[int] = [int(zeta_inv ** i) % mod for i in powers]
        rou_inverse_powers = [i if i <= mod // 2 else i - mod for i in rou_inverse_powers]
        ALREADY_COMPUTED_ZETAS_AND_INVS[(mod, n)] = rou_powers, rou_inverse_powers
    return ALREADY_COMPUTED_ZETAS_AND_INVS[(mod, n)]


ALREADY_COMPUTED_BIT_REVS: dict[tuple[int, int]] = dict()


def _bit_rev(val: int, bitstr_len: int) -> int:
    """
    Reverses the bits of val in the binary representation of length bitstr_len.
    :param val: Input value
    :type val: int
    :param bitstr_len: Input bitstring length.
    :type bitstr_len: int
    :return: An integer whose binary representation is the bitstring of the input val, just reversed.
    :rtype: int
    """
    if 0 <= val < 2**bitstr_len and (val, bitstr_len) not in ALREADY_COMPUTED_BIT_REVS:
        val_in_bin: str = bin(val)[2:].zfill(bitstr_len)
        nib_ni_lav: str = val_in_bin[::-1]
        ALREADY_COMPUTED_BIT_REVS[(val, bitstr_len)] = int(nib_ni_lav, 2)
    elif 0 > val or val >= 2**bitstr_len:
        raise ValueError
    return ALREADY_COMPUTED_BIT_REVS[(val, bitstr_len)]


def _bit_rev_cp(val: list[int]) -> list[int]:
    """
    Inputs a list of integers, permutes this list by reversing the bits of the indices, and then outputs the result.
    :param val: Input list of integers
    :type val: list[int]
    :return:
    :rtype: list[int]
    """
    if _is_int_and_gt_one_and_pow_two(val=len(val)):
        bits_per_index: int = ceil(log2(len(val)))
        return [val[_bit_rev(val=i, bitstr_len=bits_per_index)] if _bit_rev(val=i, bitstr_len=bits_per_index) < len(val) else 0 for i in range(len(val))]
    raise ValueError


def _ntt_base(mod: int, n: int, vals: list[int], inverse: bool) -> list[int]:
    """
    The NTT algorithm, which is just the iterative FFT algorithm (but modified to work over the finite field of interest) from Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2022). Introduction to algorithms. MIT press.
    :param mod: Input modulus. The algorithm takes place over the field with this many elements.
    :type mod: int
    :param n: Order of the primitive root of unity in this field.
    :type n: int
    :param vals: Input vector to be NTT'd.
    :return:
    :rtype: list[int]
    """
    if n not in LG_INTEGERS:
        LG_INTEGERS[n] = ceil(log2(n))
    lg_ord_of_prim_rou = LG_INTEGERS[n]
    zetas, zeta_invs = _make_zetas_and_invs(mod=mod, n=n)
    bit_rev_cp_vals: list[int] = _bit_rev_cp(val=vals)
    m: int = 1
    for s in range(1, lg_ord_of_prim_rou + 1):
        m *= 2
        if inverse:
            this_zeta: int = zeta_invs[s - 1]
        else:
            this_zeta: int = zetas[s - 1]
        for k in range(0, n, m):
            w: int = 1
            for j in range(m // 2):
                t: int = w * bit_rev_cp_vals[k + j + m // 2]
                u: int = bit_rev_cp_vals[k + j]
                bit_rev_cp_vals[k + j]: int = _reduce(val=u + t, mod=mod)
                bit_rev_cp_vals[k + j + m // 2]: int = _reduce(val=u - t, mod=mod)
                w *= this_zeta
    if inverse:
        n_inv: int = 1
        while _reduce(val=n_inv * n, mod=mod) != 1 and n_inv < mod:
            n_inv += 1
        bit_rev_cp_vals: list[int] = [_reduce(val=n_inv * i, mod=mod) for i in bit_rev_cp_vals]
    return bit_rev_cp_vals


# def _split_into_sublists_of_len_k(val: list[int], k: int) -> list[list[int]]:
#     """
#     A reshaping function for the NTT. Mainly only matters when (coef_mod - 1) % 2**i != 0 for some i such that 2**i < 2*deg_mod. Works by inputting a list of integers, and outputting a list of lists of integers, such that input elements with index % k == 0 occur in the 0-th list, the elements with index % k == 1 occur in the 1-st list, and so on. This function can be used as a transpose function.
#
#     :param val: Input list of integers to be reshaped.
#     :type val: int
#     :param k: Number of sublists.
#     :type k: int
#     :return:
#     :rtype: list[list[int]]
#     """
#     if k >= 1 and len(val) % k == 0:
#         return [[val[j] for j in range(len(val)) if j % (len(val) // k) == i] for i in range(len(val) // k)]
#     raise ValueError
#
#
# def _unsplit(val: list[list[int]]) -> list[int]:
#     """
#     A reshaping function for the NTT. Inverts the _split_into_sublists_of_len_k function. Note that the _split_into_sublists_of_len_k function is essentially the transpose function, but inputs a flat list. So this basically just flattens the input before re-transposing and then re-flattening..
#
#     TODO: split the flattening into its own function for testing.
#     :param val:
#     :type val: list[list[int]]
#     :return:
#     :rtype: list[int]
#     """
#     k: int = len(val[0])
#     flattened_val: list[int] = []
#     for entry in val:
#         flattened_val += entry
#     transposed_sublists: list[list[int]] = _split_into_sublists_of_len_k(val=flattened_val, k=k)
#     result: list[int] = []
#     for entry in transposed_sublists:
#         result += entry
#     return result
#
#
# def _split_into_sublists(val: Coefs) -> list[list[int]]:
#     twice_deg_mod: int = 2*val.deg_mod
#     k: int = 2*val.deg_mod//val.ord_of_prim_rou
#     padded_coefs = val.vals + [0 for _ in range(twice_deg_mod - len(val.vals))]
#     return [[padded_coefs[i] for i in range(twice_deg_mod) if i % k == j] for j in range(k)]
#
#
# def _split_sublists_to_split_coefs(coef_mod: int, deg_mod_of_ntt: int, n: int, val: list[list[int]], consts: list[int]) -> SplitCoefs:
#     ntts: list[list[int]] = [_ntt_base(mod=coef_mod, n=n, vals=v, inverse=False) for v in val]  # apply the NTTs
#     zntts: list[list[int]] = [list(v) for v in zip(*ntts)]  # zip the results
#     bitrevd_consts: list[int] = _bit_rev_cp(val=consts)  # bit-reverse the constant moduli
#     the_ntts: list[Coefs] = [  # assemble the ntts
#         Coefs(coef_mod=coef_mod, deg_mod=deg_mod_of_ntt, const=x, ord_of_prim_rou=n, vals=y) for x, y in zip(bitrevd_consts, zntts)
#     ]
#     return SplitCoefs(vals=the_ntts)
#
#
# def _merge_sublists(val: list[list[int]]) -> list[int]:
#     zipped_vals: list[tuple[int]] = list(zip(*val))
#     result = []
#     for next_zipped_tuple in zipped_vals:
#         result += list(next_zipped_tuple)
#     return result


# def _ntt_coefs_to_splitcoefs(val: Coefs) -> SplitCoefs:
#     """
#     Input a Coefs object, output a SplitCoefs object representing a (possibly partial, possibly full) NTT of the input Coefs object.
#     :param val: Input Coefs object.
#     :type val: Coefs
#     :param n: Input degree of primitive root of unity for computing the NTT
#     :type n: int
#     :return: NTT of input
#     :rtype: SplitCoefs
#     """
#     coef_mod: int = val.coef_mod
#     deg_mod: int = val.deg_mod
#     n: int = val.ord_of_prim_rou
#     deg_mod_of_ntt: int = 2 * deg_mod // n
#     rou: int
#     rou, _ = _get_nth_prim_rou_and_inv(mod=coef_mod, n=n)
#     exponents: list[int] = list(range(n))
#     if deg_mod_of_ntt > 1:
#         exponents: list[int] = [i for i in range(n) if i % deg_mod_of_ntt == 1]
#     consts: list[int] = [_reduce(val=rou ** exponent, mod=coef_mod) for exponent in exponents]  # compute the "const" attributes for each coordinate
#     split: list[list[int]] = _split_into_sublists(val=val)  # split the list into sublists for applying the NTT
#     return _split_sublists_to_split_coefs(coef_mod=coef_mod, deg_mod_of_ntt=deg_mod_of_ntt, n=n, val=split, consts=consts)
#
#
# def _intt_splitcoefs_to_coefs(val: SplitCoefs, const: int) -> Coefs:
#     """
#     Input a SplitCoefs object, output a Coefs object representing a (possibly partial, possibly full) inverse NTT of the input SplitCoefs object.
#     TODO: Change so that n is computed from val and placed into an ALREADY_COMPUTED dict.
#     TODO: Modify this and the _split function to be less spaghettified.
#     :param val: Input SplitCoefs object.
#     :type val: SplitCoefs
#     :param n: Input degree of primitive root of unity for computing the NTT
#     :type n: int
#     :return: INTT of input
#     :rtype: Coefs
#     """
#     coef_mod: int = val.vals[0].coef_mod
#     deg_mod: int = len(val.vals)*val.vals[0].deg_mod // 2
#     n: int = val.vals[0].ord_of_prim_rou
#     split_ntts: list[list[int]] = [list(v) for v in zip(*[v.vals for v in val.vals])]
#     split_intts: list[list[int]] = [_ntt_base(mod=coef_mod, n=n, vals=val, inverse=True) for val in split_ntts]
#     reduced_split_intts: list[list[int]] = [[_reduce(val=x[i] - const*x[i + n], mod=coef_mod) for i in range(deg_mod)] for x in split_intts]
#     merged_sublists: list[int] = _merge_sublists(val=reduced_split_intts)
#     return Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=1, vals=merged_sublists)
#
#
# def ntt(val: Coefs | SplitCoefs) -> SplitCoefs | Coefs:
#     """
#     General wrapper for _ntt_coefs_to_splitcoefs and _intt_splitcoefs_to_coefs.
#     """
#     if isinstance(val, Coefs):
#         return _ntt_coefs_to_splitcoefs(val=val)
#     return _intt_splitcoefs_to_coefs(val=val)
