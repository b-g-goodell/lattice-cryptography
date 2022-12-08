import pytest
from math import ceil, log2
from random import randrange
from lattice_cryptography.number_theoretic_transforms import _is_odd_prime, \
    _reduce, \
    _is_int_gt_one_and_is_pow_two, \
    Coefs, \
    PtVals, \
    _mod_has_nth_prim_rou, \
    _get_nth_prim_rou_and_inv, \
    _make_zetas_and_invs, \
    _bit_rev, \
    _bit_rev_cp, \
    _ntt_base, \
    _split_into_sublists_of_len_k, \
    ntt

SAMPLE_SIZE: int = 3  # for use later

IS_ODD_PRIME_CASES: list[tuple[int, bool]] = [
    (3, True),
    (5, True),
    (7, True),
    (11, True),
    (13, True),
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


COEFS_PARAMS_CASES: list[tuple[int, int, int, int, list[int]]] = [
                                                                     (q[0], d[0], n, 1,
                                                                      [randrange(q[0]) for i in range(d[0])]) for q in
                                                                     IS_ODD_PRIME_CASES for d in
                                                                     IS_POS_INT_AND_POW_TWO_CASES for n in
                                                                     [d[0] >> i for i in range(ceil(log2(d[0]))) for j
                                                                      in range(SAMPLE_SIZE)]
                                                                 ] + [
                                                                     (q[0], d[0], n, 2,
                                                                      [randrange(q[0]) for i in range(d[0])]) for q in
                                                                     IS_ODD_PRIME_CASES for d in
                                                                     IS_POS_INT_AND_POW_TWO_CASES for n in
                                                                     [d[0] >> i for i in range(ceil(log2(d[0]))) for j
                                                                      in range(SAMPLE_SIZE)]
                                                                 ]
COEFS_CASES: list[tuple[int, int, int, int, list[int], Coefs]] = [
    (i[0], i[1], i[2], i[3], i[4], Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4])) for
    i in COEFS_PARAMS_CASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals,initialized_object", COEFS_CASES)
def test_coefs_init_and_repr_and_eq(coef_mod, max_deg, ord_of_prim_rou, const, vals, initialized_object):
    assert isinstance(initialized_object, Coefs)
    assert initialized_object.coef_mod == coef_mod
    assert initialized_object.deg_mod == max_deg
    assert initialized_object.const == const
    assert len(initialized_object.vals) <= max_deg
    assert all((x - y) % coef_mod == 0 for x, y in zip(initialized_object.vals, vals[:max_deg + 1]))
    assert Coefs(coef_mod=coef_mod, deg_mod=max_deg, const=const, vals=vals) == initialized_object


COEFS_ADD_PREPRECASES: list[tuple[int, int, int, int, list[int], list[int]]] = [
                                                                                   (q[0], d[0], n, 1,
                                                                                    [randrange(q[0]) for i in
                                                                                     range(d[0])],
                                                                                    [randrange(q[0]) for i in
                                                                                     range(d[0])]) for q in
                                                                                   IS_ODD_PRIME_CASES for d in
                                                                                   IS_POS_INT_AND_POW_TWO_CASES for n in
                                                                                   [d[0] >> i for i in
                                                                                    range(ceil(log2(d[0]))) for j in
                                                                                    range(SAMPLE_SIZE)]
                                                                               ] + [
                                                                                   (q[0], d[0], n, 2,
                                                                                    [randrange(q[0]) for i in
                                                                                     range(d[0])],
                                                                                    [randrange(q[0]) for i in
                                                                                     range(d[0])]) for q in
                                                                                   IS_ODD_PRIME_CASES for d in
                                                                                   IS_POS_INT_AND_POW_TWO_CASES for n in
                                                                                   [d[0] >> i for i in
                                                                                    range(ceil(log2(d[0]))) for j in
                                                                                    range(SAMPLE_SIZE)]
                                                                               ]
COEFS_ADD_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
    (i[0], i[1], i[2], i[3], i[4], i[5], [_reduce(val=x + y, mod=i[0]) for x, y in zip(i[4], i[5])]) for i in
    COEFS_ADD_PREPRECASES
]
COEFS_ADD_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
    (i[0], i[1], i[2], i[3], i[4], i[5], i[6],
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[5]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[6]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=[x + y for x, y in zip(i[4], i[5])])) for
    i in COEFS_ADD_PRECASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_plus_b", COEFS_ADD_CASES)
def test_coefs_add_and_radd(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_plus_b):
    assert a + b == c == b + a
    assert a_plus_b == c
    assert all((x - y) % coef_mod == 0 for x, y in zip(a_plus_b.vals, c.vals))


COEFS_SUB_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
    (i[0], i[1], i[2], i[3], i[4], i[5], [_reduce(val=x - y, mod=i[0]) for x, y in zip(i[4], i[5])]) for i in
    COEFS_ADD_PREPRECASES
]
COEFS_SUB_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
    (i[0], i[1], i[2], i[3], i[4], i[5], i[6],
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[5]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[6]),
     Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=[x - y for x, y in zip(i[4], i[5])])) for
    i in COEFS_SUB_PRECASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_minus_b", COEFS_SUB_CASES)
def test_coefs_sub_and_neg(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_minus_b):
    assert a - b == c == -(b - a)
    assert a_minus_b == c
    assert all((x - y) % coef_mod == 0 for x, y in zip(a_minus_b.vals, c.vals))


COEFS_MUL_PREPRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int]]] = [
    (17, 4, 8, 1, [1, 2, 3, 4], [7, 8, -8, -7], [-4, 7, 6, -5]),
    (17, 4, 8, 1, [1, 0, 0, 0], [7, 8, -8, -7], [7, 8, -8, -7]),
    (17, 4, 8, 1, [0, 1, 0, 0], [7, 8, -8, -7], [7, 7, 8, -8]),
    (17, 4, 8, 1, [0, 0, 1, 0], [7, 8, -8, -7], [8, 7, 7, 8]),
    (17, 4, 8, 1, [0, 0, 0, 1], [7, 8, -8, -7], [-8, 8, 7, 7]),
    (17, 4, 8, 2, [1, 2, 3, 4], [7, 8, -8, -7], [2, 9, 0, -5]),
    (17, 4, 8, 2, [1, 0, 0, 0], [7, 8, -8, -7], [7, 8, -8, -7]),
    (17, 4, 8, 2, [0, 1, 0, 0], [7, 8, -8, -7], [-3, 7, 8, -8]),
    (17, 4, 8, 2, [0, 0, 1, 0], [7, 8, -8, -7], [-1, 14, 7, 8]),
    (17, 4, 8, 2, [0, 0, 0, 1], [7, 8, -8, -7], [1, -1, -3, 7]),
]
COEFS_MUL_PRECASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
    i + tuple([
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[4]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[5]),
        Coefs(coef_mod=i[0], deg_mod=i[1], const=i[3], vals=i[6])
    ]) for i in COEFS_MUL_PREPRECASES
]
COEFS_MUL_CASES: list[tuple[int, int, int, int, list[int], list[int], list[int], Coefs, Coefs, Coefs, Coefs]] = [
    i + tuple([
        i[-3] * i[-2]
    ]) for i in COEFS_MUL_PRECASES
]


@pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const,vals_a,vals_b,vals_c,a,b,c,a_times_b", COEFS_MUL_CASES)
def test_coefs_mul_and_rmul(coef_mod, max_deg, ord_of_prim_rou, const, vals_a, vals_b, vals_c, a, b, c, a_times_b):
    assert a * b == c
    assert b * a == c
    assert a_times_b == c
    assert all((x - y) % coef_mod == 0 for x, y in zip(a_times_b.vals, c.vals))


PTVALS_PREPRECASES: list[tuple[int, int, int, int, list[list[int]]]] = [
    (q[0], d[0], n, _get_nth_prim_rou_and_inv(mod=q[0], n=n)[0],
     [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)]) for q in
    IS_ODD_PRIME_CASES for d in IS_POS_INT_AND_POW_TWO_CASES for n in
    [2 ** i for i in range(ceil(log2(d[0])) + 1, 0, -1)]
]
PTVALS_PRECASES: list[tuple[int, int, int, int, list[list[int]], list[Coefs]]] = [
    i + tuple([
        [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
         enumerate(i[-1])]
    ]) for i in PTVALS_PREPRECASES
]
PTVALS_CASES: list[tuple[int, int, int, int, list[list[int]], list[Coefs], PtVals]] = [
    i + tuple([PtVals(vals=i[-1])]) for i in PTVALS_PRECASES
]


@pytest.mark.parametrize("coef_mod,max_deg,ord_of_prim_rou,const, coefs,vals,initialized_object", PTVALS_CASES)
def test_ptvals_init_and_repr_and_eq(coef_mod, max_deg, ord_of_prim_rou, const, coefs, vals, initialized_object):
    assert initialized_object == PtVals(vals=vals)


PTVALS_ADD_PREPREPRECASES: list[tuple[int, int, int, int, list[list[int]], list[list[int]]]] = [
    (q[0], d[0], n, _get_nth_prim_rou_and_inv(mod=q[0], n=n)[0],
     [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)],
     [[_reduce(val=randrange(q[0]), mod=q[0]) for k in range((2 * d[0]) // n)] for j in range(n)]) for q in
    IS_ODD_PRIME_CASES for d in IS_POS_INT_AND_POW_TWO_CASES for n in
    [2 ** i for i in range(ceil(log2(d[0])) + 1, 0, -1)]
]
PTVALS_ADD_PREPRECASES: list[
    tuple[int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs]]] = [
    i + tuple([
        [[_reduce(val=x + y, mod=i[0]) for x, y in zip(a, b)] for a, b in zip(i[-1], i[-2])],
        [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
         enumerate(i[-2])],
        [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
         enumerate(i[-1])]
    ]) for i in PTVALS_ADD_PREPREPRECASES
]
PTVALS_ADD_PRECASES: list[tuple[
    int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
        Coefs], PtVals, PtVals]] = [
    i + tuple([
        [Coefs(coef_mod=i[0], deg_mod=len(j) - 1, const=i[3] ** (k + 1), vals=j) for k, j in
         enumerate(i[-3])],
        PtVals(vals=i[-2]),
        PtVals(vals=i[-1])
    ]) for i in PTVALS_ADD_PREPRECASES
]
PTVALS_ADD_CASES: list[tuple[
    int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
        Coefs], PtVals, PtVals, PtVals, PtVals]] = [
    i + tuple([
        PtVals(vals=i[-3]),
        i[-1] + i[-2]
    ]) for i in PTVALS_ADD_PRECASES
]


@pytest.mark.parametrize(
    "coef_mod,max_deg,ord_of_prim_rou,const,coefs_a,coefs_b,coefs_a_plus_b,vals_a,vals_b,vals_a_plus_b,a,b,c,a_plus_b",
    PTVALS_ADD_CASES)
def test_ptvals_add(coef_mod, max_deg, ord_of_prim_rou, const, coefs_a, coefs_b, coefs_a_plus_b, vals_a, vals_b,
                    vals_a_plus_b, a, b, c, a_plus_b):
    assert a + b == a_plus_b == b + a == c


# (17, 4, 8, _get_nth_prim_rou_and_inv(mod=17, n=8)[0], [[1], [2], [3], [4], [5], [6], [7], [8]], [[7], [8], [-8], [-7],
# [-6], [-5], [-4], [-3]], [[7], [-1], [-7], [6], [4], [4], [6], [-7]]),
PTVALS_MUL_PREPREPRECASES: list[tuple[int, int, int, int, list[list[int]], list[list[int]], list[list[int]]]] = [
    (17, 4, 4, _get_nth_prim_rou_and_inv(mod=17, n=4)[0], [[1, 2], [3, 4], [5, 6], [7, 8]],
     [[7, 8], [-8, -7], [-6, -5], [-4, -3]], [[-6, 5], [-1, -2], [3, 7], [-4, -2]]),
]
PTVALS_MUL_PREPRECASES: list[tuple[
    int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[Coefs]]] = [
    i + tuple([
        [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2], const=i[3] ** (k + 1), vals=vals) for
         k, vals in enumerate(i[4])],
        [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2], const=i[3] ** (k + 1), vals=vals) for
         k, vals in enumerate(i[5])],
        [Coefs(coef_mod=i[0], deg_mod=(2 * i[1]) // i[2],  const=i[3] ** (k + 1), vals=vals) for
         k, vals in enumerate(i[6])]
    ]) for i in PTVALS_MUL_PREPREPRECASES
]
PTVALS_MUL_PRECASES: list[tuple[
    int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
        Coefs], PtVals, PtVals, PtVals]] = [
    i + tuple([
        PtVals(vals=i[7]),
        PtVals(vals=i[8]),
        PtVals(vals=i[9]),
    ]) for i in PTVALS_MUL_PREPRECASES
]
PTVALS_MUL_CASES: list[tuple[
    int, int, int, int, list[list[int]], list[list[int]], list[list[int]], list[Coefs], list[Coefs], list[
        Coefs], PtVals, PtVals, PtVals, PtVals]] = [
    i + tuple([
        i[-2] * i[-3]
    ]) for i in PTVALS_MUL_PRECASES
]
PTVALS_MUL_CASES = [
    i + tuple([
        i[-1].vals
    ]) for i in PTVALS_MUL_CASES
]


@pytest.mark.parametrize(
    "coef_mod,max_deg,ord_of_prim_rou,rou,coefs_a,coefs_b,coefs_c,vals_a,vals_b,vals_c,a,b,c,a_times_b,a_times_b_vals",
    PTVALS_MUL_CASES)
def test_ptvals_mul(coef_mod, max_deg, ord_of_prim_rou, rou, coefs_a, coefs_b, coefs_c, vals_a, vals_b, vals_c, a, b, c,
                    a_times_b, a_times_b_vals):
    assert a * b == c
    assert b * a == c
    assert a_times_b == c


def test_split_into_sublists_of_len_k():
    val: list[int] = list(range(16))
    k: int = 1
    result: list[list[int]] = _split_into_sublists_of_len_k(val=val, k=k)
    assert isinstance(result, list)
    assert all(isinstance(v, list) for v in result)
    assert all(len(v) == k for v in result)
    assert all(isinstance(w, int) for v in result for w in v)
    assert result == [[i] for i in range(16)]

    k: int = 2
    result: list[list[int]] = _split_into_sublists_of_len_k(val=val, k=k)
    assert isinstance(result, list)
    assert all(isinstance(v, list) for v in result)
    assert all(len(v) == k for v in result)
    assert all(isinstance(w, int) for v in result for w in v)
    assert result == [[i, i + len(val) // k] for i in range(len(val) // k)]


def test_coefs_mul_matches_ntt_mul():
    coef_mod: int = 17
    deg_mod: int = 8
    ord_of_prim_rou: int = 8  # this is a half-NTT
    const: int = 1

    one_vals: list[int] = [1] + [0 for _ in range(deg_mod - 1)]
    one_poly: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=one_vals)

    some_vals: list[int] = [randrange(coef_mod) for _ in range(deg_mod)]
    a_poly: Coefs = Coefs(coef_mod=coef_mod, deg_mod=deg_mod, const=const, vals=some_vals)

    assert one_poly * a_poly == a_poly

    one_ntt: PtVals = ntt(val=one_poly, ord_of_prim_rou=ord_of_prim_rou)
    assert all(isinstance(v, Coefs) for v in one_ntt.vals)
    expected_degree_mod_of_vals: int = 2 * deg_mod // ord_of_prim_rou
    assert all(len(v.vals) == expected_degree_mod_of_vals for v in one_ntt.vals)
    assert all(v.coef_mod == coef_mod for v in one_ntt.vals)
    assert all(v.deg_mod <= expected_degree_mod_of_vals for v in one_ntt.vals)
    # assert all(v.ord_of_prim_rou == ord_of_prim_rou for v in one_ntt.vals)
    assert all(isinstance(v, Coefs) for v in one_ntt.vals)
    assert all(isinstance(v.vals, list) for v in one_ntt.vals)
    assert all(len(v.vals) == 2 for v in one_ntt.vals)
    assert all(v.vals == [1, 0] for v in one_ntt.vals)

    some_ntt: PtVals = ntt(val=a_poly, ord_of_prim_rou=ord_of_prim_rou)
    assert all(isinstance(v, Coefs) for v in some_ntt.vals)
    expected_degree_mod_of_vals: int = 2 * deg_mod // ord_of_prim_rou
    assert all(len(v.vals) == expected_degree_mod_of_vals for v in some_ntt.vals)
    assert all(v.coef_mod == coef_mod for v in some_ntt.vals)
    assert all(v.deg_mod == expected_degree_mod_of_vals for v in some_ntt.vals)
    # assert all(v.ord_of_prim_rou == ord_of_prim_rou for v in some_ntt.vals)
    assert all(isinstance(v, Coefs) for v in some_ntt.vals)
    assert all(isinstance(v.vals, list) for v in some_ntt.vals)
    assert all(len(v.vals) == expected_degree_mod_of_vals for v in some_ntt.vals)

    product: PtVals = one_ntt * some_ntt
    intt_product: Coefs = ntt(val=product, ord_of_prim_rou=ord_of_prim_rou)

    assert intt_product == a_poly
    assert product == some_ntt
