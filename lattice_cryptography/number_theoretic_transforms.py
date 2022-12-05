# Handles generalized NTT
from math import ceil, sqrt, log2
from copy import deepcopy
from warnings import warn


def _is_prime(val: int) -> bool:
    # Sieve of Eratosthenes
    return all([val % i != 0 for i in range(2, ceil(sqrt(val)) + 1)])


HALF_INTEGERS: dict[int] = dict()
LG_INTEGERS: dict[int] = dict()


def _reduce(val: int, mod: int) -> int:
    # Montgomery reduction
    if mod not in HALF_INTEGERS:
        HALF_INTEGERS[mod] = mod // 2
    if mod not in LG_INTEGERS:
        LG_INTEGERS[mod] = ceil(log2(mod))
    half_mod = HALF_INTEGERS[mod]
    lg_mod = LG_INTEGERS[mod]
    return (val % mod) - (1 + (((val % mod) - half_mod - 1) >> lg_mod)) * mod


def _is_pos_int_and_pow_two(val: int) -> bool:
    return isinstance(val, int) and val > 0 and not (val & (val - 1))


class Coefs(object):
    # Coefficient representation of a polynomial mod X**max_deg + 1
    coef_mod: int  # coefficient modulus
    half_coef_mod: int  # greatest integer such that 2*half_coefficient_modulus <= coefficient_modulus
    lg_coef_mod: int
    max_deg: int  # maximum degree
    ord_of_prim_rou: int  # order of primitive root-of-unity used to compute the NTT representation
    lg_ord_of_prim_rou: int
    vals: list[int]  # coefficient representation of the polynomial

    def __init__(self, coef_mod: int, max_deg: int, ord_of_prim_rou: int, vals: list[int]):
        self.coef_mod = coef_mod
        if coef_mod not in HALF_INTEGERS:
            HALF_INTEGERS[coef_mod] = coef_mod // 2
        self.half_coef_mod = HALF_INTEGERS[coef_mod]
        if coef_mod not in LG_INTEGERS:
            LG_INTEGERS[coef_mod] = ceil(log2(coef_mod))
        self.lg_coef_mod: int = LG_INTEGERS[coef_mod]
        self.max_deg: int = max_deg
        self.ord_of_prim_rou = ord_of_prim_rou
        if ord_of_prim_rou not in LG_INTEGERS:
            LG_INTEGERS[ord_of_prim_rou] = ceil(log2(ord_of_prim_rou))
        self.lg_ord_of_prim_rou = LG_INTEGERS[ord_of_prim_rou]
        self.vals = [_reduce(val=x, mod=coef_mod) for x in vals]

    def __repr__(self) -> str:
        return f'Coefs{(self.coef_mod, self.max_deg, self.ord_of_prim_rou, self.vals[:self.max_deg+1])}'

    def __eq__(self, other) -> bool:
        return self.__repr__() == other.__repr__()

    def __add__(self, other):
        # check that both have the same max_deg and modulus
        result = deepcopy(self)
        result.vals = [_reduce(val=x + y, mod=self.coef_mod) for x, y in zip(self.vals, other.vals)]
        return result

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [-x for x in negative_other.vals]
        return self + negative_other

    def __mul__(self, other):
        # check that both have the same max_deg and modulus
        result = deepcopy(self)
        result.vals = []
        for i in range(self.max_deg + 1):
            result.vals += [0]
            for j in range(i+1):
                result.vals[-1] += self.vals[j]*other.vals[i-j]
            result.vals[-1] = _reduce(val=result.vals[i], mod=self.coef_mod)
        lower_vals = result.vals[:self.max_deg+1]
        upper_vals = result.vals[self.max_deg+1:]
        result.vals = [x - y for x, y in zip(lower_vals, upper_vals)]
        return result

    def __rmul__(self, other):
        return self.__mul__(other=other)


class PtVals(object):
    # Coefficient representation of a polynomial mod X**max_deg + 1
    coef_mod: int  # coefficient modulus
    half_coef_mod: int  # greatest integer such that 2*half_coefficient_modulus <= coefficient_modulus
    lg_coef_mod: int
    max_deg: int  # maximum degree
    ord_of_prim_rou: int  # order of primitive root-of-unity used to compute the NTT representation
    lg_ord_of_prim_rou: int
    vals: list[Coefs]  # "point-value" representation of the polynomial

    def __init__(self, coef_mod: int, max_deg: int, ord_of_prim_rou: int, vals: list[Coefs]):
        self.coef_mod = coef_mod
        if coef_mod not in HALF_INTEGERS:
            HALF_INTEGERS[coef_mod] = coef_mod // 2
        self.half_coef_mod = HALF_INTEGERS[coef_mod]
        if coef_mod not in LG_INTEGERS:
            LG_INTEGERS[coef_mod] = ceil(log2(coef_mod))
        self.lg_coef_mod: int = LG_INTEGERS[coef_mod]
        self.max_deg: int = max_deg
        self.ord_of_prim_rou = ord_of_prim_rou
        if ord_of_prim_rou not in LG_INTEGERS:
            LG_INTEGERS[ord_of_prim_rou] = ceil(log2(ord_of_prim_rou))
        self.lg_ord_of_prim_rou = LG_INTEGERS[ord_of_prim_rou]
        self.vals = vals  # check that len(vals) * ord_of_prim_rou == 2*(max_deg + 1)?

    def __repr__(self) -> str:
        return f'Coefs{(self.coef_mod, self.max_deg, self.ord_of_prim_rou, self.vals)}'

    def __eq__(self, other) -> bool:
        return self.__repr__() == other.__repr__()

    def __add__(self, other):
        result = deepcopy(self)
        result.vals = [x + y for x, y in zip(self.vals, other.vals)]
        return result

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [-x for x in negative_other.vals]
        return self + negative_other

    def __mul__(self, other):
        result = deepcopy(self)
        result.vals = [x * y for x, y in zip(self.vals, other.vals)]
        return result

    def __rmul__(self, other):
        return self.__mul__(other=other)


def _mod_has_nth_prim_rou(mod: int, n: int) -> bool:
    return isinstance(mod, int) and mod > 1 and isinstance(n, int) and n > 1 and (mod - 1) % n == 0


PRIM_NTH_ROUS_AND_ROU_INVS: dict[tuple[int, int]] = dict()


def _get_nth_prim_rou_and_inv(mod: int, n: int) -> tuple[int, int]:
    if (mod, n) not in PRIM_NTH_ROUS_AND_ROU_INVS:
        rou: int = 2
        while rou < mod:
            if _reduce(val=rou ** n, mod=mod) == 1 and all(_reduce(val=rou ** k, mod=mod) != 1 for k in range(1, n)):
                break
            rou += 1
        rou_inv = _reduce(val=rou ** (n - 1), mod=mod)
        PRIM_NTH_ROUS_AND_ROU_INVS[(mod, n)] = rou, rou_inv
    return PRIM_NTH_ROUS_AND_ROU_INVS[(mod, n)]


ALREADY_COMPUTED_ZETAS_AND_INVS: dict[tuple[int, int]] = dict()


def _make_zetas_and_invs(mod: int, ord_of_prim_rou: int) -> tuple[list[int], list[int]]:
    if (mod, ord_of_prim_rou) not in ALREADY_COMPUTED_ZETAS_AND_INVS:
        lg_rou_deg: int = ceil(log2(ord_of_prim_rou))
        powers: list[int] = [ord_of_prim_rou // (2 ** (s + 1)) for s in range(lg_rou_deg)]
        zeta, zeta_inv = _get_nth_prim_rou_and_inv(mod=mod, n=ord_of_prim_rou)
        rou_powers: list[int] = [int(zeta ** i) % mod for i in powers]
        rou_powers = [i if i <= mod // 2 else i - mod for i in rou_powers]
        rou_inverse_powers: list[int] = [int(zeta_inv ** i) % mod for i in powers]
        rou_inverse_powers = [i if i <= mod // 2 else i - mod for i in rou_inverse_powers]
        ALREADY_COMPUTED_ZETAS_AND_INVS[(mod, ord_of_prim_rou)] = rou_powers, rou_inverse_powers
    return ALREADY_COMPUTED_ZETAS_AND_INVS[(mod, ord_of_prim_rou)]


ALREADY_COMPUTED_BIT_REVS: dict[tuple[int, int]] = dict()


def _bit_rev(val: int, bitstr_len: int) -> int:
    # assumes 0 <= val < 2**len_of_val_as_bitstr
    if (val, bitstr_len) not in ALREADY_COMPUTED_BIT_REVS:
        val_in_bin: str = bin(val)[2:].zfill(bitstr_len)
        nib_ni_lav: str = val_in_bin[::-1]
        ALREADY_COMPUTED_BIT_REVS[(val, bitstr_len)] = int(nib_ni_lav, 2)
    return ALREADY_COMPUTED_BIT_REVS[(val, bitstr_len)]


def _bit_rev_cp(val: list[int]) -> list[int]:
    # assumes len(val) is a power-of-two
    bits_per_index: int = ceil(log2(len(val)))
    return [val[_bit_rev(val=i, bitstr_len=bits_per_index)] for i in range(len(val))]


def _ntt_base(mod: int, ord_of_prim_rou: int, vals: list[int], inverse: bool) -> list[int]:
    if ord_of_prim_rou not in LG_INTEGERS:
        LG_INTEGERS[ord_of_prim_rou] = ceil(log2(ord_of_prim_rou))
    lg_ord_of_prim_rou = LG_INTEGERS[ord_of_prim_rou]
    zetas, zeta_invs = _make_zetas_and_invs(mod=mod, ord_of_prim_rou=ord_of_prim_rou)
    bit_rev_cp_vals: list[int] = _bit_rev_cp(val=vals)
    m: int = 1
    for s in range(1, lg_ord_of_prim_rou + 1):
        m *= 2
        if inverse:
            this_zeta: int = zeta_invs[s - 1]
        else:
            this_zeta: int = zetas[s - 1]
        for k in range(0, ord_of_prim_rou, m):
            w: int = 1
            for j in range(m // 2):
                t: int = w * bit_rev_cp_vals[k + j + m // 2]
                u: int = bit_rev_cp_vals[k + j]
                bit_rev_cp_vals[k + j]: int = _reduce(val=u + t, mod=mod)
                bit_rev_cp_vals[k + j + m // 2]: int = _reduce(val=u - t, mod=mod)
                w *= this_zeta
    if inverse:
        n_inv: int = 1
        while _reduce(val=n_inv * ord_of_prim_rou, mod=mod) != 1:
            n_inv += 1
        bit_rev_cp_vals: list[int] = [_reduce(val=n_inv * i, mod=mod) for i in bit_rev_cp_vals]
    return bit_rev_cp_vals


def _ntt(val: Coefs | PtVals, ord_of_prim_rou: int) -> PtVals | Coefs:
    q: int = val.coef_mod  # rename for shorter lines of code
    n: int = ord_of_prim_rou
    if isinstance(val, Coefs):
        k: int = 2 * len(val.vals) // ord_of_prim_rou
        pad_len: int = 2*val.max_deg + 1
        pad_vals: list[int] = val.vals + [0 for _ in range(pad_len - len(val.vals))]
        split_vals: list[list[int]] = [[pad_vals[j] for j in range(pad_len) if j % k == i] for i in range(k)]
        split_ntts: list[list[int]] = [_ntt_base(mod=q, ord_of_prim_rou=n, vals=x, inverse=False) for x in split_vals]
        ntt_vals: list[Coefs] = [Coefs(coef_mod=q, max_deg=n - 1, vals=[x[i] for x in split_ntts]) for i in range(n)]
        return PtVals(coef_mod=q, max_deg=val.max_deg, ord_of_prim_rou=n, vals=ntt_vals)
    k: int = len(val.vals) // ord_of_prim_rou
    split_ntts: list[list[int]] = [[val.vals[j].vals[i] for j in range(k)] for i in range(n)]
    split_vals: list[list[int]] = [_ntt_base(mod=q, ord_of_prim_rou=n, vals=x, inverse=True) for x in split_ntts]
    pad_vals: list[int] = [split_vals[i][j] for j in range(len(split_vals[0])) for i in range(len(split_vals))]
    lower_vals = pad_vals[:val.max_deg + 1]
    upper_vals = pad_vals[val.max_deg + 1:]
    coefs_vals = [_reduce(val=x - y, mod=q) for x, y in zip(lower_vals, upper_vals)]
    return Coefs(coef_mod=q, max_deg=val.max_deg, ord_of_prim_rou=n, vals=coefs_vals)
