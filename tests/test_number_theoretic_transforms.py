import pytest
from math import ceil, log2
from random import randrange
from lattice_cryptography.number_theoretic_transforms import _is_odd_prime, \
    _reduce, \
    _is_int_gt_one_and_is_pow_two, \
    Coefs, \
    SplitCoefs, \
    _mod_has_nth_prim_rou, \
    _get_nth_prim_rou_and_inv, \
    _make_zetas_and_invs, \
    _bit_rev, \
    _bit_rev_cp, \
    _ntt_base, \
    _split_into_sublists_of_len_k, \
    ntt

SAMPLE_SIZE: int = 32  # for use later

IS_ODD_PRIME_CASES: list[tuple[int, bool]] = [
    (3, True),
    (5, True),
    (7, True),
    (11, True),
    (13, True),
    (17, True),
    (3329, True),  # kyber modulus
    (8380417, True),  # dilithium modulus
    (8675309, True),  # Jennys number
]
IS_NOT_ODD_PRIME_CASES: list[tuple[int, bool]] = [
    (4, False),
    (6, False),
    (8, False),
    (9, False),
    (10, False),
    (12, False),
    (3328, False),
    (3330, False),
    (8380416, False),
    (8380418, False),
    (8675308, False),
    (8675310, False)
]


# @pytest.mark.skip()
@pytest.mark.parametrize("val,expected_output", IS_ODD_PRIME_CASES + IS_NOT_ODD_PRIME_CASES)
def test_is_prime(val, expected_output):
    assert _is_odd_prime(val=val) == expected_output


REDUCE_CASES: list[tuple[int, int, int]] = [
    (i, q[0], i if i <= q[0] // 2 else i - q[0]) for q in IS_ODD_PRIME_CASES[:6] for i in range(q[0])
]


# @pytest.mark.skip()
@pytest.mark.parametrize("val,mod,expected_output", REDUCE_CASES)
def test_reduce(val, mod, expected_output):
    assert _reduce(val=val, mod=mod) == expected_output


IS_POS_INT_AND_POW_TWO_CASES: list[tuple[int, bool]] = [
    (2 ** i, True) for i in range(3, SAMPLE_SIZE + 3)
]
IS_NOT_POS_INT_OR_NOT_POW_TWO_CASES: list[tuple[int, bool]] = [
    (2 ** i + 2 ** j, False) for i in range(3, SAMPLE_SIZE + 3) for j in range(1, i)
]


# @pytest.mark.skip()
@pytest.mark.parametrize("val,expected_output", IS_POS_INT_AND_POW_TWO_CASES + IS_NOT_POS_INT_OR_NOT_POW_TWO_CASES)
def test_is_pos_int_and_pow_two(val, expected_output):
    assert _is_int_gt_one_and_is_pow_two(val=val) == expected_output


COEFS_PRECASES: list[tuple[int, int, int, int, list[int]]] = [
    (3329, 128, 1, 256, list(range(i, 128+i))) for i in range(SAMPLE_SIZE)
] + [
    (3329, 256, 1, 256, list(range(i, 256+i))) for i in range(SAMPLE_SIZE)
]
COEFS_CASES: list[tuple[int, int, int, int, list[int], Coefs]] = [
    i + tuple([
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[4])
    ]) for i in COEFS_PRECASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("coef_mod,deg_mod,const,n,vals,initialized_object", COEFS_CASES)
def test_coefs_init_and_repr_and_eq(coef_mod, deg_mod, const, n, vals, initialized_object):
    assert isinstance(initialized_object, Coefs)
    assert initialized_object.coef_mod == coef_mod
    assert initialized_object.deg_mod == deg_mod
    assert initialized_object.const == const
    assert len(initialized_object.vals) <= deg_mod
    assert all((x - y) % coef_mod == 0 for x, y in zip(initialized_object.vals, vals))
    assert str(initialized_object) == f'Coefs{(coef_mod, deg_mod, const, vals)}'
    assert Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=[i + coef_mod for i in vals]) == initialized_object


COEFS_ARITHMETIC_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
    (3329, 128, 1, 256, list(range(i, 128+i, 1)), list(range(512 + i, 640 + i, 1)), list(range(512+2*i, 768 + 2*i, 2))) for i in range(SAMPLE_SIZE)
] + [
    (3329, 256, 1, 256, list(range(i, 128+i, 1)), list(range(512 + i, 640 + i, 1)), list(range(512+2*i, 768 + 2*i, 2))) for i in range(SAMPLE_SIZE)
]
COEFS_ARITHMETIC_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs]] = [
    i + tuple([
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[4]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[5]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[6])
    ]) for i in COEFS_ARITHMETIC_PRECASES
]


@pytest.mark.parametrize("coef_mod,deg_mod,const,ord_of_prim_rou,vals_a,vals_b,vals_c,coefs_a,coefs_b,coefs_c", COEFS_ARITHMETIC_CASES)
def test_coefs_arithmetic(coef_mod, deg_mod, const, ord_of_prim_rou, vals_a, vals_b, vals_c, coefs_a, coefs_b, coefs_c):
    assert isinstance(coefs_a, Coefs)
    assert isinstance(coefs_b, Coefs)
    assert all(_.coef_mod == coef_mod for _ in [coefs_a, coefs_b, coefs_c])
    assert all(_.deg_mod == deg_mod for _ in [coefs_a, coefs_b, coefs_c])
    assert all(_.const == const for _ in [coefs_a, coefs_b, coefs_c])

    assert coefs_a.vals == vals_a
    assert coefs_b.vals == vals_b
    assert coefs_c.vals == vals_c

    observed_coefs_a_plus_b: Coefs = coefs_a + coefs_b
    assert observed_coefs_a_plus_b == coefs_c
    observed_coefs_b_plus_a: Coefs = coefs_b + coefs_a
    assert observed_coefs_b_plus_a == observed_coefs_a_plus_b

    neg_coefs_a: Coefs = -coefs_a
    neg_coefs_b: Coefs = -coefs_b
    neg_coefs_c: Coefs = -coefs_c

    assert neg_coefs_a.vals == [-_ for _ in coefs_a.vals]
    assert neg_coefs_b.vals == [-_ for _ in coefs_b.vals]
    assert neg_coefs_c.vals == [-_ for _ in coefs_c.vals]

    zero: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=[0 for _ in range(deg_mod)])
    zero_a: Coefs = coefs_a + neg_coefs_a
    zero_b: Coefs = coefs_b + neg_coefs_b
    zero_c: Coefs = coefs_c + neg_coefs_c

    assert zero == zero_a == zero_b == zero_c

    coefs_a_minus_b: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=[a - b for a, b in zip(coefs_a.vals, coefs_b.vals)])
    coefs_a_minus_coefs_b: Coefs = coefs_a - coefs_b
    assert coefs_a_minus_b == coefs_a_minus_coefs_b

    coefs_b_minus_a: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=[b - a for a, b in zip(coefs_a.vals, coefs_b.vals)])
    coefs_b_minus_coefs_a: Coefs = coefs_b - coefs_a
    assert coefs_b_minus_a == coefs_b_minus_coefs_a

    coefs_ab_vals: list[int] = [0 for _ in range(2*deg_mod)]
    for i, a_coef in enumerate(coefs_a.vals):
        for j, b_coef in enumerate(coefs_b.vals):
            coefs_ab_vals[i+j] += a_coef*b_coef
    lower: list[int] = coefs_ab_vals[:deg_mod]
    upper: list[int] = coefs_ab_vals[deg_mod:]
    coefs_ab_vals = [x - y for x, y in zip(lower, upper)]
    coefs_ab: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=coefs_ab_vals)

    assert coefs_ab == coefs_a*coefs_b == coefs_b*coefs_a


COEFS_MONOMIALS_ARITHMETIC_PREPREPRECASES: list[tuple[int, int, int, int, int, int, int, int]] = [
    (3329, 256, 1, 256, randrange(256), randrange(256), randrange(3329), randrange(3329)) for i in range(SAMPLE_SIZE)
]
COEFS_MONOMIALS_ARITHMETIC_PREPRECASES: list[tuple[int, int, int, int, int, int, int, int, list[int], list[int], list[int]]] = [
    i + tuple([
        [i[5] if j == i[3] else 0 for j in range(256)],
        [i[6] if j == i[4] else 0 for j in range(256)],
        [(2*int(i[3]+i[4]<256)-1)*((i[5]*i[6]) % 3329) if j == (i[3]+i[4]) % 256 else 0 for j in range(256)]
    ]) for i in COEFS_MONOMIALS_ARITHMETIC_PREPREPRECASES
]
COEFS_MONOMIALS_ARITHMETIC_PRECASES: list[tuple[int, int, int, int, int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs]] = [
    i + tuple([
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[-3]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[-2]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[2], ord_of_prim_rou=i[3], vals=i[-1])
    ]) for i in COEFS_MONOMIALS_ARITHMETIC_PREPRECASES
]
COEFS_MONOMIALS_ARITHMETIC_CASES: list[tuple[int, int, int, int, int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
    i + tuple([
        i[-3]*i[-2]
    ]) for i in COEFS_MONOMIALS_ARITHMETIC_PRECASES
]


@pytest.mark.parametrize("coef_mod,deg_mod,const,n,idx_a,idx_b,int_a,int_b,vals_a,vals_b,vals_ab,coefs_a,coefs_b,coefs_ab,coefs_a_times_coefs_b", COEFS_MONOMIALS_ARITHMETIC_CASES)
def test_monomial_products(coef_mod,deg_mod,const,n,idx_a,idx_b,int_a,int_b,vals_a,vals_b,vals_ab,coefs_a,coefs_b,coefs_ab,coefs_a_times_coefs_b):
    assert coefs_ab == coefs_a_times_coefs_b


MOD_HAS_NTH_PRIM_ROU_CASES: list[tuple[int, int, bool]] = [
    (17, 2, True),
    (17, 3, False),
    (17, 4, True),
    (17, 5, False),
    (17, 6, False),
    (17, 7, False),
    (17, 8, True),
    (17, 9, False),
    (17, 10, False),
    (17, 11, False),
    (17, 12, False),
    (17, 13, False),
    (17, 14, False),
    (17, 15, False),
    (17, 16, True)
]


@pytest.mark.parametrize("mod,n,expected_output", MOD_HAS_NTH_PRIM_ROU_CASES)
def test_mod_has_nth_prim_rou(mod, n, expected_output):
    assert _mod_has_nth_prim_rou(mod=mod, n=n) == expected_output


GET_NTH_PRIM_ROU_AND_INV_PRECASES: list[tuple[int, int, int, int]] = [
    (17, 2, 16, -1),
    (17, 4, 4, -4),
    (17, 8, 2, -8),
    (17, 16, 3, 6)
]
GET_NTH_PRIM_ROU_AND_INV_CASES: list[tuple[int, int, int, int, int, int]] = [
    i + _get_nth_prim_rou_and_inv(mod=i[0], n=i[1]) for i in GET_NTH_PRIM_ROU_AND_INV_PRECASES
]


@pytest.mark.parametrize("mod,n,expected_rou,expected_rou_inv,observed_rou,observed_rou_inv", GET_NTH_PRIM_ROU_AND_INV_CASES)
def test_get_nth_prim_rou_and_inv(mod,n,expected_rou,expected_rou_inv,observed_rou,observed_rou_inv):
    assert _get_nth_prim_rou_and_inv(mod=mod, n=n) == (observed_rou, observed_rou_inv)
    assert expected_rou == observed_rou
    assert expected_rou_inv == observed_rou_inv
    assert (observed_rou * observed_rou_inv) % mod == 1
    assert all((observed_rou**i % mod) != 1 for i in range(1, n))
    assert (observed_rou**n) % mod == 1


SPLIT_COEFS_CASES: list[tuple[Coefs, Coefs, SplitCoefs]] = [
    (i[-1], j[-1], SplitCoefs(vals=[i[-1], j[-1]])) for i, j in zip(COEFS_CASES, COEFS_CASES[::-1])
]


@pytest.mark.parametrize("coefs_a,coefs_b,split_coefs", SPLIT_COEFS_CASES)
def test_split_coefs_init_and_repr_and_eq(coefs_a, coefs_b, split_coefs):
    assert split_coefs == SplitCoefs(vals=[coefs_a, coefs_b])
    assert isinstance(split_coefs.vals[0], Coefs)
    assert isinstance(split_coefs.vals[1], Coefs)
    assert str(split_coefs) == f'SplitCoefs{tuple([coefs_a, coefs_b])}'


SPLIT_COEFS_ARITHMETIC_PRECASES: list[tuple[Coefs, Coefs, Coefs, Coefs, SplitCoefs, SplitCoefs, SplitCoefs, SplitCoefs]] = [
    (i[-3], i[-2], j[-3], j[-2],
     SplitCoefs(vals=[i[-3], i[-2]]),
     SplitCoefs(vals=[j[-3], j[-2]]),
     SplitCoefs(vals=[i[-3]+j[-3], i[-2]+j[-2]]),
     SplitCoefs(vals=[i[-3]*j[-3], i[-2]*j[-2]]))
    for i, j in zip(COEFS_ARITHMETIC_CASES[:len(COEFS_ARITHMETIC_CASES)//2], (COEFS_ARITHMETIC_CASES[:len(COEFS_ARITHMETIC_CASES)//2])[::-1])
]
SPLIT_COEFS_ARITHMETIC_CASES = [
    i + tuple([
        i[4] + i[5],
        i[4] * i[5]
    ]) for i in SPLIT_COEFS_ARITHMETIC_PRECASES
]


@pytest.mark.parametrize("coefs_a,coefs_b,coefs_c,coefs_d,split_coefs_a_and_b,split_coefs_c_and_d,expected_sum,expected_product,observed_sum,observed_product", SPLIT_COEFS_ARITHMETIC_CASES)
def test_split_coefs_arithmetic(coefs_a,coefs_b,coefs_c,coefs_d,split_coefs_a_and_b,split_coefs_c_and_d,expected_sum,expected_product,observed_sum,observed_product):
    assert split_coefs_a_and_b == SplitCoefs(vals=[coefs_a, coefs_b])
    assert split_coefs_c_and_d == SplitCoefs(vals=[coefs_c, coefs_d])
    assert expected_sum == split_coefs_a_and_b + split_coefs_c_and_d
    assert expected_product == split_coefs_a_and_b * split_coefs_c_and_d


# coef_mod,deg_mod,const,idx_a,idx_b,int_a,int_b,vals_a,vals_b,vals_ab,coefs_a,coefs_b,coefs_ab,coefs_a_times_coefs_b

@pytest.mark.parametrize("coef_mod,deg_mod,const,n,idx_a,idx_b,int_a,int_b,vals_a,vals_b,vals_ab,coefs_a,coefs_b,coefs_ab,coefs_a_times_coefs_b", COEFS_MONOMIALS_ARITHMETIC_CASES)
def test_monomial_products_via_ntt(coef_mod,deg_mod,const,n,idx_a,idx_b,int_a,int_b,vals_a,vals_b,vals_ab,coefs_a,coefs_b,coefs_ab,coefs_a_times_coefs_b):
    ntt_a = ntt(coefs_a)
    ntt_b = ntt(coefs_b)
    ntt_ab = ntt(coefs_ab)
    ntt_a_times_ntt_b = ntt_a*ntt_b
    assert ntt_ab == ntt_a_times_ntt_b
    intt_ntt_a_times_ntt_b = ntt(ntt_a_times_ntt_b)
    assert intt_ntt_a_times_ntt_b == coefs_a * coefs_b














# # @pytest.mark.skip()
# @pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_plus_b", COEFS_ADD_CASES)
# def test_coefs_add_and_radd(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_plus_b):
#     assert a + b == c == b + a
#     assert a_plus_b == c
#     assert all((x - y) % coef_mod == 0 for x, y in zip(a_plus_b.vals, c.vals))
#
#
# COEFS_SUB_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
#     (i[0], i[1], i[2], i[3], i[4], i[5], [_reduce(val=x - y, mod=i[0]) for x, y in zip(i[4], i[5])]) for i in
#     COEFS_ADD_PREPRECASES
# ]
# COEFS_SUB_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
#     (i[0], i[1], i[2], i[3], i[4], i[5], i[6],
#      Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4]),
#      Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[5]),
#      Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[6]),
#      Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=[x - y for x, y in zip(i[4], i[5])])) for
#     i in COEFS_SUB_PRECASES
# ]
#
#
# # @pytest.mark.skip()
# @pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_minus_b", COEFS_SUB_CASES)
# def test_coefs_sub_and_neg(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_minus_b):
#     assert a - b == c == -(b - a)
#     assert a_minus_b == c
#     assert all((x - y) % coef_mod == 0 for x, y in zip(a_minus_b.vals, c.vals))
#
#
# COEFS_MUL_PREPRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
#     (17, 4, 8, 1, [1, 2, 3, 4], [7, 8, -8, -7], [-4, 7, 6, -5]),
#     (17, 4, 8, 1, [1, 0, 0, 0], [7, 8, -8, -7], [7, 8, -8, -7]),
#     (17, 4, 8, 1, [0, 1, 0, 0], [7, 8, -8, -7], [7, 7, 8, -8]),
#     (17, 4, 8, 1, [0, 0, 1, 0], [7, 8, -8, -7], [8, 7, 7, 8]),
#     (17, 4, 8, 1, [0, 0, 0, 1], [7, 8, -8, -7], [-8, 8, 7, 7]),
#     (17, 4, 8, 2, [1, 2, 3, 4], [7, 8, -8, -7], [2, 9, 0, -5]),
#     (17, 4, 8, 2, [1, 0, 0, 0], [7, 8, -8, -7], [7, 8, -8, -7]),
#     (17, 4, 8, 2, [0, 1, 0, 0], [7, 8, -8, -7], [-3, 7, 8, -8]),
#     (17, 4, 8, 2, [0, 0, 1, 0], [7, 8, -8, -7], [-1, 14, 7, 8]),
#     (17, 4, 8, 2, [0, 0, 0, 1], [7, 8, -8, -7], [1, -1, -3, 7]),
# ]
# COEFS_MUL_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
#     i + tuple([
#         Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4]),
#         Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[5]),
#         Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[6])
#     ]) for i in COEFS_MUL_PREPRECASES
# ]
# COEFS_MUL_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
#     i + tuple([
#         i[-3] * i[-2]
#     ]) for i in COEFS_MUL_PRECASES
# ]
#
#
# @pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_times_b", COEFS_MUL_CASES)
# def test_coefs_mul_and_rmul(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_times_b):
#     assert a * b == c
#     assert b * a == c
#     assert a_times_b == c
#     assert all((x - y) % coef_mod == 0 for x, y in zip(a_times_b.vals, c.vals))
#
#
# PTVALS_PREPRECASES: list[tuple[int, int, int, int, list[list[int]]]] = [
#     (q[0], d[0], n, _get_nth_prim_rou_and_inv(mod=q[0], n=n)[0],
#      [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)]) for q in
#     IS_ODD_PRIME_CASES for d in IS_POS_INT_AND_POW_TWO_CASES for n in
#     [2 ** i for i in range(ceil(log2(d[0])) + 1, 0, -1)]
# ]
# PTVALS_PRECASES: list[tuple[int, int, int, int, list[list[int]], list[Coefs]]] = [
#     i + tuple([
#         [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
#          enumerate(i[-1])]
#     ]) for i in PTVALS_PREPRECASES
# ]
# PTVALS_CASES: list[tuple[int, int, int, int, list[list[int]], list[Coefs], SplitCoefs]] = [
#     i + tuple([SplitCoefs(vals=i[-1])]) for i in PTVALS_PRECASES
# ]
#
#
# @pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const, coefs,vals,initialized_object", PTVALS_CASES)
# def test_ptvals_init_and_repr_and_eq(coef_mod, max_deg, ord_of_prim_rou, const, coefs, vals, initialized_object):
#     assert initialized_object == SplitCoefs(vals=vals)
#
#
# PTVALS_ADD_PREPREPRECASES: list[tuple[int, int, int, int, list[list[int]], list[list[int]]]] = [
#     (q[0], d[0], n, _get_nth_prim_rou_and_inv(mod=q[0], n=n)[0],
#      [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)],
#      [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)]) for q in
#     IS_ODD_PRIME_CASES for d in IS_POS_INT_AND_POW_TWO_CASES for n in
#     [2 ** i for i in range(ceil(log2(d[0])) + 1, 0, -1)]
# ]
# PTVALS_ADD_PREPRECASES: list[
#     tuple[int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs]]] = [
#     i + tuple([
#         [[_reduce(val=x + y, mod=i[0]) for x, y in zip(a, b)] for a, b in zip(i[-1], i[-2])],
#         [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
#          enumerate(i[-2])],
#         [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
#          enumerate(i[-1])]
#     ]) for i in PTVALS_ADD_PREPREPRECASES
# ]
# PTVALS_ADD_PRECASES: list[tuple[
#     int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
#         Coefs], SplitCoefs, SplitCoefs]] = [
#     i + tuple([
#         [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
#          enumerate(i[-3])],
#         SplitCoefs(vals=i[-2]),
#         SplitCoefs(vals=i[-1])
#     ]) for i in PTVALS_ADD_PREPRECASES
# ]
# PTVALS_ADD_CASES: list[tuple[
#     int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
#         Coefs], SplitCoefs, SplitCoefs, SplitCoefs, SplitCoefs]] = [
#     i + tuple([
#         SplitCoefs(vals=i[-3]),
#         i[-1] + i[-2]
#     ]) for i in PTVALS_ADD_PRECASES
# ]
#
#
# @pytest.mark.parametrize(
#     "coef_mod,max_deg,ord_of_prim_rou,const,coefs_a,coefs_b,coefs_a_plus_b,vals_a,vals_b,vals_a_plus_b,a,b,c,a_plus_b",
#     PTVALS_ADD_CASES)
# def test_ptvals_add(coef_mod, max_deg, ord_of_prim_rou, const, coefs_a, coefs_b, coefs_a_plus_b, vals_a, vals_b,
#                     vals_a_plus_b, a, b, c, a_plus_b):
#     assert a + b == a_plus_b == b + a == c
#
#
# # (17, 4, 8, _get_nth_prim_rou_and_inv(mod=17, n=8)[0], [[1], [2], [3], [4], [5], [6], [7], [8]], [[7], [8], [-8], [-7],
# # [-6], [-5], [-4], [-3]], [[7], [-1], [-7], [6], [4], [4], [6], [-7]]),
# PTVALS_MUL_PREPREPRECASES: list[tuple[int, int, int, int, list[list[int]], list[list[int]], list[list[int]]]] = [
#     (17, 4, 4, _get_nth_prim_rou_and_inv(mod=17, n=4)[0], [[1, 2], [3, 4], [5, 6], [7, 8]],
#      [[7, 8], [-8, -7], [-6, -5], [-4, -3]], [[-6, 5], [-1, -2], [3, 7], [-4, -2]]),
# ]
# PTVALS_MUL_PREPRECASES: list[tuple[
#     int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[Coefs]]] = [
#     i + tuple([
#         [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2], const=i[3] ** (k + 1), vals=vals) for
#          k, vals in enumerate(i[4])],
#         [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2], const=i[3] ** (k + 1), vals=vals) for
#          k, vals in enumerate(i[5])],
#         [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2],  const=i[3] ** (k + 1), vals=vals) for
#          k, vals in enumerate(i[6])]
#     ]) for i in PTVALS_MUL_PREPREPRECASES
# ]
# PTVALS_MUL_PRECASES: list[tuple[
#     int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
#         Coefs], SplitCoefs, SplitCoefs, SplitCoefs]] = [
#     i + tuple([
#         SplitCoefs(vals=i[7]),
#         SplitCoefs(vals=i[8]),
#         SplitCoefs(vals=i[9]),
#     ]) for i in PTVALS_MUL_PREPRECASES
# ]
# PTVALS_MUL_CASES: list[tuple[
#     int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
#         Coefs], SplitCoefs, SplitCoefs, SplitCoefs, SplitCoefs]] = [
#     i + tuple([
#         i[-2] * i[-3]
#     ]) for i in PTVALS_MUL_PRECASES
# ]
# PTVALS_MUL_CASES = [
#     i + tuple([
#         i[-1].vals
#     ]) for i in PTVALS_MUL_CASES
# ]
#
#
# @pytest.mark.parametrize(
#     "coef_mod,max_deg,ord_of_prim_rou,rou,coefs_a,coefs_b,coefs_c,vals_a,vals_b,vals_c,a,b,c,a_times_b,a_times_b_vals",
#     PTVALS_MUL_CASES)
# def test_ptvals_mul(coef_mod, max_deg, ord_of_prim_rou, rou, coefs_a, coefs_b, coefs_c, vals_a, vals_b, vals_c, a, b, c,
#                     a_times_b, a_times_b_vals):
#     assert a * b == c
#     assert b * a == c
#     assert a_times_b == c
#
#
# def test_split_into_sublists_of_len_k():
#     val: list[int] = list(range(16))
#     k: int = 1
#     result: list[list[int]] = _split_into_sublists_of_len_k(val=val, k=k)
#     assert isinstance(result, list)
#     assert all(isinstance(v, list) for v in result)
#     assert all(len(v) == k for v in result)
#     assert all(isinstance(w, int) for v in result for w in v)
#     assert result == [[i] for i in range(16)]
#
#     k: int = 2
#     result: list[list[int]] = _split_into_sublists_of_len_k(val=val, k=k)
#     assert isinstance(result, list)
#     assert all(isinstance(v, list) for v in result)
#     assert all(len(v) == k for v in result)
#     assert all(isinstance(w, int) for v in result for w in v)
#     assert result == [[i, i + len(val) // k] for i in range(len(val) // k)]
#
#
# def test_coefs_mul_matches_ntt_mul():
#     coef_mod: int = 17
#     deg_mod: int = 8
#     ord_of_prim_rou: int = 8  # this is a half-NTT
#     const: int = 1
#
#     one_vals: list[int] = [1] + [0 for _ in range(deg_mod - 1)]
#     one_poly: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=one_vals)
#
#     some_vals: list[int] = [randrange(coef_mod) for _ in range(deg_mod)]
#     a_poly: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=some_vals)
#
#     assert one_poly * a_poly == a_poly
#
#     one_ntt: SplitCoefs = ntt(val=one_poly, ord_of_prim_rou=ord_of_prim_rou)
#     assert all(isinstance(v, Coefs) for v in one_ntt.vals)
#     expected_degree_mod_of_vals: int = 2 * deg_mod // ord_of_prim_rou
#     assert all(len(v.vals) == expected_degree_mod_of_vals for v in one_ntt.vals)
#     assert all(v.coef_mod == coef_mod for v in one_ntt.vals)
#     assert all(v.deg_mod <= expected_degree_mod_of_vals for v in one_ntt.vals)
#     # assert all(v.ord_of_prim_rou == ord_of_prim_rou for v in one_ntt.vals)
#     assert all(isinstance(v, Coefs) for v in one_ntt.vals)
#     assert all(isinstance(v.vals, list) for v in one_ntt.vals)
#     assert all(len(v.vals) == 2 for v in one_ntt.vals)
#     assert all(v.vals == [1, 0] for v in one_ntt.vals)
#
#     some_ntt: SplitCoefs = ntt(val=a_poly, ord_of_prim_rou=ord_of_prim_rou)
#     assert all(isinstance(v, Coefs) for v in some_ntt.vals)
#     expected_degree_mod_of_vals: int = 2 * deg_mod // ord_of_prim_rou
#     assert all(len(v.vals) == expected_degree_mod_of_vals for v in some_ntt.vals)
#     assert all(v.coef_mod == coef_mod for v in some_ntt.vals)
#     assert all(v.deg_mod == expected_degree_mod_of_vals for v in some_ntt.vals)
#     # assert all(v.ord_of_prim_rou == ord_of_prim_rou for v in some_ntt.vals)
#     assert all(isinstance(v, Coefs) for v in some_ntt.vals)
#     assert all(isinstance(v.vals, list) for v in some_ntt.vals)
#     assert all(len(v.vals) == expected_degree_mod_of_vals for v in some_ntt.vals)
#
#     product: SplitCoefs = one_ntt * some_ntt
#     intt_product: Coefs = ntt(val=product, ord_of_prim_rou=ord_of_prim_rou)
#
#     assert intt_product == a_poly
#     assert product == some_ntt
