# Handles generalized NTT
from math import ceil, sqrt, log2
from copy import deepcopy


def _is_odd_prime(val: int) -> bool:
    # Sieve of Eratosthenes
    return all([val % i != 0 for i in range(2, ceil(sqrt(val)) + 1)])


HALF_INTEGERS: dict[int] = dict()
LG_INTEGERS: dict[int] = dict()


def _reduce(val: int, mod: int) -> int:
    # Montgomery reduction
    if mod not in HALF_INTEGERS:
        HALF_INTEGERS[mod] = mod // 2
    half_mod = HALF_INTEGERS[mod]
    if mod not in LG_INTEGERS:
        LG_INTEGERS[mod] = ceil(log2(mod))
    lg_mod = LG_INTEGERS[mod]
    return (val % mod) - (1 + (((val % mod) - half_mod - 1) >> lg_mod)) * mod


def _is_int_gt_one_and_is_pow_two(val: int) -> bool:
    return isinstance(val, int) and \
        val > 1 and \
        not (val & (val - 1))


class Coefs(object):
    # Coefficient representation of a polynomial mod (coef_mod, X**deg_mod + const) where
    # deg_mod is power of 2, coef_mod is a prime, and (coef_mod-1) % ord_of_prim_rou == 0
    coef_mod: int  # coefficient modulus
    deg_mod: int  # maximum degree
    const: int
    vals: list[int]  # coefficient representation of the polynomial

    def __init__(self, coef_mod: int, deg_mod: int, const: int, vals: list[int]):
        self.coef_mod = coef_mod
        self.deg_mod: int = deg_mod
        self.const = const
        self.vals = [_reduce(val=x, mod=coef_mod) for x in vals]

    def __repr__(self) -> str:
        return f'Coefs(coef_mod={self.coef_mod}, deg_mod={self.deg_mod}, const={self.const}, vals={self.vals})'

    def __eq__(self, other) -> bool:
        return isinstance(other, Coefs) and \
            self.coef_mod == other.coef_mod and \
            self.deg_mod == other.deg_mod and \
            self.const == other.const and \
            all((x-y) % self.coef_mod == 0 for x, y in zip(self.vals, other.vals))

    def __add__(self, other):
        # check that both have the same max_deg and modulus
        result = deepcopy(self)
        result.vals = [_reduce(val=x + y, mod=self.coef_mod) for x, y in zip(self.vals, other.vals)]
        return result

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
        negative_other.vals = [-x for x in negative_other.vals]
        return self + negative_other

    def __mul__(self, other):
        # check that both have the same max_deg and modulus
        result = deepcopy(self)
        result.vals = [0 for _ in range(2 * self.deg_mod)]
        for i, a_coef in enumerate(self.vals):
            for j, b_coef in enumerate(other.vals):
                result.vals[i+j] += a_coef * b_coef
        lower_vals = result.vals[:self.deg_mod]
        upper_vals = result.vals[self.deg_mod:]
        result.vals = [_reduce(val=x - y * self.const, mod=result.coef_mod) for x, y in zip(lower_vals, upper_vals)]
        return result

    def __rmul__(self, other):
        return self.__mul__(other=other)


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


class PtVals(object):
    # "Point-Value" representation of a polynomial mod (coef_mod, X**deg_mod + const) where
    # deg_mod is power of 2, coef_mod is a prime, and (coef_mod-1) % ord_of_prim_rou == 0
    # Not really a point-value representation.
    vals: list[Coefs]  # "point-value" representation of the polynomial

    def __init__(self, vals: list[Coefs]):
        self.vals = vals  # check that len(vals) * ord_of_prim_rou == 2*(max_deg + 1)?

    def __repr__(self) -> str:
        return f'PtVals(vals={[entry for entry in self.vals]})'

    def __eq__(self, other) -> bool:
        return isinstance(other, PtVals) and all(x == y for x, y in zip(self.vals, other.vals))

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


def _mod_has_nth_prim_rou(mod: int, ord_of_prim_rou: int) -> bool:
    return isinstance(mod, int) and \
        mod > 1 and \
        isinstance(ord_of_prim_rou, int) and \
        ord_of_prim_rou > 1 and \
        (mod - 1) % ord_of_prim_rou == 0


ALREADY_COMPUTED_ZETAS_AND_INVS: dict[tuple[int, int]] = dict()


def _make_zetas_and_invs(mod: int, ord_of_prim_rou: int) -> tuple[list[int], list[int]]:
    if (mod, ord_of_prim_rou) not in ALREADY_COMPUTED_ZETAS_AND_INVS:
        if ord_of_prim_rou not in LG_INTEGERS:
            LG_INTEGERS[ord_of_prim_rou] = ceil(log2(ord_of_prim_rou))
        lg_rou_deg = LG_INTEGERS[ord_of_prim_rou]
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
        while _reduce(val=n_inv * ord_of_prim_rou, mod=mod) != 1 and n_inv < mod:
            n_inv += 1
        bit_rev_cp_vals: list[int] = [_reduce(val=n_inv * i, mod=mod) for i in bit_rev_cp_vals]
    return bit_rev_cp_vals


def _split_into_sublists_of_len_k(val: list[int], k: int) -> list[list[int]]:
    return [[val[j] for j in range(len(val)) if j % (len(val) // k) == i] for i in range(len(val) // k)]


def _unsplit(val: list[list[int]]) -> list[int]:
    k: int = len(val[0])
    flattened_val: list[int] = []
    for entry in val:
        flattened_val += entry
    transposed_sublists: list[list[int]] = _split_into_sublists_of_len_k(val=flattened_val, k=k)
    result: list[int] = []
    for entry in transposed_sublists:
        result += entry
    return result


def _ntt_coefs_to_ptvals(val: Coefs, ord_of_prim_rou: int) -> PtVals:
    coef_mod: int = val.coef_mod
    deg_mod: int = val.deg_mod
    deg_mod_of_ntt: int = 2 * val.deg_mod // ord_of_prim_rou
    rou: int
    rou_inv: int
    rou, rou_inv = _get_nth_prim_rou_and_inv(mod=coef_mod, n=ord_of_prim_rou)
    padvals: list[int] = val.vals
    while len(padvals) < 2 * deg_mod:
        padvals += [0]
    split: list[list[int]] = _split_into_sublists_of_len_k(val=padvals, k=ord_of_prim_rou)
    ntts: list[list[int]] = [
        _ntt_base(mod=coef_mod, ord_of_prim_rou=ord_of_prim_rou, vals=val, inverse=False) for val in split
    ]
    zntts: list[list[int]] = [list(v) for v in zip(*ntts)]
    consts: list[int] = [_reduce(val=rou ** (2 * i + 1), mod=coef_mod) for i in range(len(zntts))]
    bitrevd_consts: list[int] = _bit_rev_cp(val=consts)
    the_ntts: list[Coefs] = [
        Coefs(coef_mod=coef_mod, deg_mod=deg_mod_of_ntt, const=x, vals=y) for x, y in zip(bitrevd_consts, zntts)
    ]
    return PtVals(vals=the_ntts)


def _ntt_ptvals_to_coefs(val: PtVals, ord_of_prim_rou: int) -> Coefs:
    coef_mod: int = val.vals[0].coef_mod
    deg_mod: int = len(val.vals)*val.vals[0].deg_mod // 2
    zippedntts: list[list[int]] = [v.vals for v in val.vals]
    splitntts: list[list[int]] = [list(v) for v in zip(*zippedntts)]
    unsplitvals: list[int] = _unsplit(val=splitntts)
    lower_vals: list[int] = unsplitvals[:deg_mod]
    upper_vals: list[int] = unsplitvals[deg_mod:]
    merged_vals: list[int] = [_reduce(val=x - y, mod=coef_mod) for x, y in zip(lower_vals, upper_vals)]
    return Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=1, vals=merged_vals)


def ntt(val: Coefs | PtVals, ord_of_prim_rou: int) -> PtVals | Coefs:
    if isinstance(val, Coefs):
        return _ntt_coefs_to_ptvals(val=val, ord_of_prim_rou=ord_of_prim_rou)
    return _ntt_ptvals_to_coefs(val=val, ord_of_prim_rou=ord_of_prim_rou)
