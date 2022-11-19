from crystals.kyber import SEED_LEN_IN_BYTES, Q, bit_rev, is_pow_two, _bit_rev_cp, bit_rev_cp, _reduce, reduce, _round_up, round_up, N, LOG_Q, _our_parse_one, _our_parse_many, K, parse, is_arithmetic_legal, PolyCoefs, PolyNTT, _cbd_eta, cbd_eta, _cbd_polycoefs, cbd_polycoefs, _compress_one_int, _decompress_one_int, _should_compress_many, compress, decompress, _encode_m_list_of_ints_to_bytes, _encode_m_many, encode_m, _decode_m_list_of_ints_from_bytes, decode_m, cpa_pke_keygen, CPA_PKE_SK_LEN, CPA_PKE_PK_LEN, CPA_PKE_CIPHERTEXT_LEN, _cpa_pke_enc, _encode_m_matrix, _cpa_pke_dec, _ntt_one, ntt
from random import getrandbits, randrange
import pytest
from math import ceil, log2


SAMPLE_SIZE: int = 2**10
LOG_SAMPLE_SIZE: int = 10


def test_bit_rev():
    with pytest.raises(TypeError):
        bit_rev(x='hello world', length_in_bits=LOG_SAMPLE_SIZE)

    with pytest.raises(TypeError):
        bit_rev(x='hello world', length_in_bits=0.001)

    with pytest.raises(TypeError):
        bit_rev(x='hello world', length_in_bits='hello world')

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length_in_bits=LOG_SAMPLE_SIZE)

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length_in_bits=0.001)

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length_in_bits='hello world')

    with pytest.raises(ValueError):
        bit_rev(x=1, length_in_bits=0)

    with pytest.raises(ValueError):
        bit_rev(x=7, length_in_bits=2)

    length_in_bytes = max(1, ceil(LOG_SAMPLE_SIZE/8))
    for i in range(SAMPLE_SIZE):
        reversed_i = bit_rev(x=i, length_in_bits=LOG_SAMPLE_SIZE)
        i_as_bits = bin(i)[2:].zfill(LOG_SAMPLE_SIZE)
        reversed_i_as_bits = bin(reversed_i)[2:].zfill(LOG_SAMPLE_SIZE)
        for j in range(len(i_as_bits)):
            assert i_as_bits[j] == reversed_i_as_bits[LOG_SAMPLE_SIZE - 1 - j]


IS_POW_TWO_CASES = [
    (2, True),
    (4, True),
    (8, True),
    (16, True),
    (3, False),
    (5, False),
    (9, False),
    (17, False),
    (0.01, False),
    ('Hello world', False)
]


@pytest.mark.parametrize("x,expected_output", IS_POW_TWO_CASES)
def test_is_pow_two(x, expected_output):
    assert is_pow_two(x=x) == expected_output


BIT_REV_CP_CASES = [
    ([0, 1], [0, 1]),
    ([0, 1, 2, 3], [0, 2, 1, 3]),
    ([0, 1, 2, 3, 4, 5, 6, 7], [0, 4, 2, 6, 1, 5, 3, 7]),
    ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15])
]


@pytest.mark.parametrize("x,expected_output", BIT_REV_CP_CASES)
def test_bit_rev_cp(x, expected_output):
    assert _bit_rev_cp(x=x, length_in_bits=ceil(log2(len(x)))) == expected_output


def test_bit_rev_cp_failures():
    with pytest.raises(TypeError):
        bit_rev_cp(x=3, length_in_bits=2)

    with pytest.raises(TypeError):
        bit_rev_cp(x=['hello world', 0.001], length_in_bits=2)

    with pytest.raises(TypeError):
        bit_rev_cp(x=['hello world', 0.001], length_in_bits=0.001)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), length_in_bits=0)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), length_in_bits=2)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), length_in_bits=4)

@pytest.mark.parametrize("x,expected_output", BIT_REV_CP_CASES)
def test_bit_rev_cp_full(x, expected_output):
    assert bit_rev_cp(x=x, length_in_bits=ceil(log2(len(x)))) == expected_output


REDUCE_CASES = [(_, Q, _) if _ <= Q//2 else (_, Q, _ - Q) for _ in list(range(Q))]

# @pytest.mark.parametrize("x, q, expected_output", REDUCE_CASES)
# def test_reduce(x, q, expected_output):
#     assert _reduce(x=x, q=q) == expected_output
#     assert reduce(x=x, q=q) == expected_output


def test_reduce_fail():
    with pytest.raises(TypeError):
        reduce(x='Hello world')

    with pytest.raises(TypeError):
        reduce(x=0.001)


def test_round_up():
    assert _round_up(x=0.5) == 1
    assert _round_up(x=0.499) == 0
    for x in range(10):
        assert _round_up(x=x) == x

    assert round_up(x=0.5) == 1
    assert round_up(x=0.499) == 0
    for x in range(10):
        assert round_up(x=x) == x

    with pytest.raises(TypeError):
        assert round_up(x='Hello world')

    with pytest.raises(TypeError):
        assert round_up(x=None)


def test_parse_one():
    test_integers = [randrange(2**LOG_Q) for _ in range(10*N)]
    test_integers_to_bytes = []
    for i in range(0, len(test_integers)-1, 2):
        first_byte = test_integers[i] % 256

        tmp = test_integers[i+1] % 16
        third_byte = (test_integers[i+1] - tmp)//16

        second_byte_mod_16 = (test_integers[i] - first_byte)//256
        second_byte_div_by_16_rd_down = test_integers[i+1] - 16*third_byte
        second_byte = second_byte_mod_16 + 16*second_byte_div_by_16_rd_down

        test_integers_to_bytes += [first_byte, second_byte, third_byte]

    resulting_ints: list[int]
    index: int
    resulting_ints, index = _our_parse_one(x=bytes(test_integers_to_bytes))
    expected_ints = [_ for _ in test_integers if _ < Q]
    expected_ints = expected_ints[:N]
    assert len(resulting_ints) == N
    assert all(isinstance(x, int) for x in resulting_ints)
    assert all(0 <= x < Q for x in resulting_ints)
    for i, (next_result_int, next_test_int) in enumerate(zip(resulting_ints, expected_ints)):
        assert next_result_int == next_test_int

    with pytest.raises(RuntimeError):
        _our_parse_one(x=b'hello world')


def test_parse_many():
    test_integers = [randrange(2**LOG_Q) for _ in range(10*K*K*N)]
    test_integers_to_bytes = []
    for i in range(0, len(test_integers)-1, 2):
        first_byte = test_integers[i] % 256

        tmp = test_integers[i+1] % 16
        third_byte = (test_integers[i+1] - tmp)//16

        second_byte_mod_16 = (test_integers[i] - first_byte)//256
        second_byte_div_by_16_rd_down = test_integers[i+1] - 16*third_byte
        second_byte = second_byte_mod_16 + 16*second_byte_div_by_16_rd_down

        test_integers_to_bytes += [first_byte, second_byte, third_byte]

    resulting_ints: list[list[list[int]]] = _our_parse_many(x=bytes(test_integers_to_bytes))
    # We only check shape, size, and modulus because we call _parse_one multiple times to make _parse_many
    assert isinstance(resulting_ints, list)
    assert len(resulting_ints) == K
    assert all(isinstance(x, list) for x in resulting_ints)
    assert all(len(x) == K for x in resulting_ints)
    assert all(isinstance(y, list) for x in resulting_ints for y in x)
    assert all(len(y) == N for x in resulting_ints for y in x)
    assert all(isinstance(z, int) for x in resulting_ints for y in x for z in y)
    assert all(0 <= z < Q for x in resulting_ints for y in x for z in y)

    with pytest.raises(RuntimeError):
        _our_parse_many(x=b'hello world')


def test_parse():
    with pytest.raises(TypeError):
        parse(x='hello world')
    with pytest.raises(ValueError):
        parse(x=b'hello world')


def test_is_arithmetic_legal():
    assert is_arithmetic_legal(a_vals = [[[]]], b_vals=[[[]]], q_a=17, q_b=17)
    assert not is_arithmetic_legal(a_vals=[[[1]]], b_vals=[[[]]], q_a=17, q_b=17)
    assert not is_arithmetic_legal(a_vals=[[[1]]], b_vals=[[[1]]], q_a=17, q_b=19)


def test_add():
    # too simple to bother testing
    pass


def test_polynomial_mul():
    identity_matrix_coefs: list[list[list[int]]] = [[[0] * N] * 1] * 1
    identity_matrix_coefs[0][0][0] = 1
    identity_matrix: PolyCoefs = PolyCoefs(vals=identity_matrix_coefs, q=Q, n=N, k1=1, k2=1)
    identity_matrix_ntt: PolyNTT = ntt(identity_matrix)
    assert len(identity_matrix_ntt.vals) == 1
    assert len(identity_matrix_ntt.vals[0]) == 1
    assert len(identity_matrix_ntt.vals[0][0]) == N
    assert all(z == 1 for x in identity_matrix_ntt.vals for y in x for z in y)

    new_random_matrix_coefs: list[list[list[int]]] = [[[randrange(Q) for _ in range(N)]]]
    new_random_matrix: PolyCoefs = PolyCoefs(vals=new_random_matrix_coefs, q=Q, n=N, k1=1, k2=1)
    new_random_matrix_ntt: PolyNTT = ntt(new_random_matrix)
    assert len(new_random_matrix_ntt.vals) == 1
    assert len(new_random_matrix_ntt.vals[0]) == 1
    assert len(new_random_matrix_ntt.vals[0][0]) == N

    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix
    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix_ntt
    with pytest.raises(NotImplementedError):
        new_random_matrix_ntt * identity_matrix

    id_times_new_ntt: PolyNTT = identity_matrix_ntt * new_random_matrix_ntt
    new_times_id_ntt: PolyNTT = new_random_matrix_ntt * identity_matrix_ntt
    id_times_new: PolyCoefs = ntt(id_times_new_ntt)
    new_times_id: PolyCoefs = ntt(new_times_id_ntt)

    assert id_times_new.q == new_times_id.q == new_random_matrix.q == Q
    assert id_times_new.k1 == new_times_id.k1 == new_random_matrix.k1 == 1
    assert id_times_new.k2 == new_times_id.k2 == new_random_matrix.k2 == 1
    for i in range(1):
        for j in range(1):
            for l in range(N):
                assert (id_times_new.vals[i][j][l] - new_random_matrix.vals[i][j][l]) % Q == 0
    assert id_times_new == new_random_matrix
    assert new_times_id == new_random_matrix



    identity_matrix_coefs: list[list[list[int]]] = [[[0] * N] * K] * K
    for i in range(K):
        for j in range(K):
            if i == j:
                identity_matrix_coefs[i][j][0] = 1
    identity_matrix: PolyCoefs = PolyCoefs(vals=identity_matrix_coefs, q=Q, n=N, k1=K, k2=K)
    identity_matrix_ntt: PolyNTT = ntt(identity_matrix)
    assert len(identity_matrix_ntt.vals) == K
    assert all(len(x) == K for x in identity_matrix_ntt.vals)
    assert all(len(y) == N for x in identity_matrix_ntt.vals for y in x)
    assert all(z == 1 for x in identity_matrix_ntt.vals for y in x for z in y)

    new_random_matrix_coefs: list[list[list[int]]] = [[[randrange(Q) for l in range(N)] for j in range(K)] for i in range(K)]
    new_random_matrix: PolyCoefs = PolyCoefs(vals=new_random_matrix_coefs, q=Q, n=N, k1=K, k2=K)
    new_random_matrix_ntt: PolyNTT = ntt(new_random_matrix)
    assert len(new_random_matrix_ntt.vals) == K
    assert len(new_random_matrix_ntt.vals[0]) == K
    assert len(new_random_matrix_ntt.vals[0][0]) == N

    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix
    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix_ntt
    with pytest.raises(NotImplementedError):
        new_random_matrix_ntt * identity_matrix

    id_times_new_ntt: PolyNTT = identity_matrix_ntt * new_random_matrix_ntt
    new_times_id_ntt: PolyNTT = new_random_matrix_ntt * identity_matrix_ntt
    id_times_new: PolyCoefs = ntt(id_times_new_ntt)
    new_times_id: PolyCoefs = ntt(new_times_id_ntt)

    assert id_times_new.q == new_times_id.q == new_random_matrix.q == Q
    assert id_times_new.k1 == new_times_id.k1 == new_random_matrix.k1 == K
    assert id_times_new.k2 == new_times_id.k2 == new_random_matrix.k2 == K
    for i in range(K):
        for j in range(K):
            for l in range(N):
                assert (id_times_new.vals[i][j][l] - new_random_matrix.vals[i][j][l]) % Q == 0
    assert id_times_new == new_random_matrix
    assert new_times_id == new_random_matrix


#
# # @pytest.mark.skip
# def test_polynomial_mul_small_without_const_time(one_without_const_time,
#                                                  pairs_of_random_polys_and_their_products_without_const_time):
#     lp = lp_for_testing
#     # First, let's make an identity polynomials and add it to itself
#     assert one_without_const_time * one_without_const_time == one_without_const_time
#
#     # Now let's do some addition with some random linear polynomials (AND the unity)
#     for next_item in pairs_of_random_polys_and_their_products_without_const_time:
#         f_dat, g_dat, expected_h_dat, observed_h_dat = next_item
#         a_f, b_f, f = f_dat
#         a_g, b_g, g = g_dat
#         a_h, b_h, c_h, exp_h, exp_h_norm, exp_h_wt = expected_h_dat
#         obs_h_coefs, obs_h, obs_h_norm, obs_h_wt = observed_h_dat
#         assert one_without_const_time * f == f
#         assert f * one_without_const_time == f
#         assert one_without_const_time * g == g
#         assert g * one_without_const_time == g
#         assert f * g == exp_h == obs_h
#         assert len(obs_h_coefs) == 3
#         assert 0 in obs_h_coefs
#         assert 1 in obs_h_coefs
#         assert 2 in obs_h_coefs
#         assert (obs_h_coefs[0] - a_h) % lp.modulus == 0
#         assert (obs_h_coefs[1] - b_h) % lp.modulus == 0
#         assert (obs_h_coefs[2] - c_h) % lp.modulus == 0
#         assert (a_f * a_g - a_h) % lp.modulus == 0
#         assert (a_f * b_g + b_f * a_g - b_h) % lp.modulus == 0
#         assert (b_f * b_g - c_h) % lp.modulus == 0
#

def test_polycoefs():
    x: PolyCoefs = PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    assert x.q == 17
    assert x.n == 2
    assert x.k1 == 3
    assert x.k2 == 2

    y: PolyCoefs = PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[5, 6], [7, 8]], [[9, 10], [11, 12]], [[13, 14], [15, 16]]])
    assert y.q == 17
    assert y.n == 2
    assert y.k1 == 3
    assert y.k2 == 2

    z: PolyCoefs = x + y
    assert z.q == 17
    assert z.n == 2
    assert z.k1 == 3
    assert z.k2 == 2
    assert z.vals == [[[5, 7], [-8, -6]], [[-4, -2], [0, 2]], [[4, 6], [8, -7]]]

    # Some failure tests. not complete, but sufficient for now.
    with pytest.raises(TypeError):
        PolyCoefs(q='hello world', n=2, k1=3, k2=2, vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n='hello world', k1=3, k2=2, vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1='hello world', k2=2, vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2='hello world', vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[['hello world', 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 'hello world'], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 1], ['hello world', 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 1], [2, 'hello world']], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[['hello world', 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 'hello world'], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 1], [2, 3]], [['hello world', 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(TypeError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[0, 1], [2, 3]], [[4, 'hello world'], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(ValueError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[17, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])
    with pytest.raises(ValueError):
        PolyCoefs(q=17, n=2, k1=3, k2=2, vals=[[[-9, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]])


def test_polyntt():
    x: PolyNTT = PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    assert x.q == 17
    assert x.n == 2
    assert x.k1 == 3
    assert x.k2 == 2

    y: PolyNTT = PolyNTT(vals=[[[5, 6], [7, 8]], [[9, 10], [11, 12]], [[13, 14], [15, 16]]], q=17, n=2, k1=3, k2=2)
    assert y.q == 17
    assert y.n == 2
    assert y.k1 == 3
    assert y.k2 == 2

    z: PolyNTT = x + y
    assert z.q == 17
    assert z.n == 2
    assert z.k1 == 3
    assert z.k2 == 2
    assert z.vals == [[[5, 7], [-8, -6]], [[-4, -2], [0, 2]], [[4, 6], [8, -7]]]

    x: PolyNTT = PolyNTT(vals=[[[3, 4]]], q=17, n=2, k1=1, k2=1)
    assert x.q == 17
    assert x.n == 2
    assert x.k1 == 1
    assert x.k2 == 1

    z: PolyNTT = x * y
    assert z.q == 17
    assert z.n == 2
    assert z.k1 == 3
    assert z.k2 == 2
    assert z.vals == [[[-2, 7], [4, -2]], [[-7, 6], [-1, -3]], [[5, 5], [-6, -4]]]

    x: PolyNTT = PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    assert x.q == 17
    assert x.n == 2
    assert x.k1 == 3
    assert x.k2 == 2
    y: PolyNTT = PolyNTT(vals=[[[1, 2]], [[3, 4]]], q=17, n=2, k1=2, k2=1)
    assert y.q == 17
    assert y.n == 2
    assert y.k1 == 2
    assert y.k2 == 1
    z: PolyNTT = x * y
    assert z.q == 17
    assert z.n == 2
    assert z.k1 == 3
    assert z.k2 == 1
    assert z.vals == [[[6, -3]], [[5, 4]], [[4, -6]]]

    # Some failure tests. not complete, but sufficient for now.
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q='hello world', n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n='hello world', k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1='hello world', k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2='hello world')
    with pytest.raises(TypeError):
        PolyNTT(vals=[[['hello world', 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 'hello world'], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], ['hello world', 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 'hello world']], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[['hello world', 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 'hello world'], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [['hello world', 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(TypeError):
        PolyNTT(vals=[[[0, 1], [2, 3]], [[4, 'hello world'], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(ValueError):
        PolyNTT(vals=[[[17, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)
    with pytest.raises(ValueError):
        PolyNTT(vals=[[[-9, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]], q=17, n=2, k1=3, k2=2)


def test_cbd_eta():
    test_eta: int = 1
    give_me_a_one: str = '10' * N
    give_me_a_one_as_bytes: bytes = bytes(int(give_me_a_one[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_one)//8))
    expected_result: list[int] = [1] * N
    assert _cbd_eta(x=give_me_a_one_as_bytes, eta=test_eta) == expected_result
    assert cbd_eta(x=give_me_a_one_as_bytes, eta=test_eta) == expected_result

    test_eta: int = 2
    give_me_a_two: str = '1100' * N
    give_me_a_two_as_bytes: bytes = bytes(int(give_me_a_two[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_two)//8))
    expected_result: list[int] = [2] * N
    assert _cbd_eta(x=give_me_a_two_as_bytes, eta=test_eta) == expected_result
    assert cbd_eta(x=give_me_a_two_as_bytes, eta=test_eta) == expected_result

    give_me_a_one: str = '1000' * N
    give_me_a_one_as_bytes: bytes = bytes(int(give_me_a_one[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_one)//8))
    expected_result: list[int] = [1] * N
    assert _cbd_eta(x=give_me_a_one_as_bytes, eta=test_eta) == expected_result
    assert cbd_eta(x=give_me_a_one_as_bytes, eta=test_eta) == expected_result

    test_eta: int = 3
    give_me_a_negative_two: str = '000011' * N
    give_me_a_negative_two_as_bytes: bytes = bytes(int(give_me_a_negative_two[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_negative_two)//8))
    expected_result: list[int] = [-2] * N
    assert _cbd_eta(x=give_me_a_negative_two_as_bytes, eta=test_eta) == expected_result
    assert cbd_eta(x=give_me_a_negative_two_as_bytes, eta=test_eta) == expected_result


def test_cbd_polycoefs():
    test_num_rows: int = 2
    test_num_cols: int = 3
    test_eta: int = 3
    give_me_a_negative_two: str = ''
    for _ in range(N*test_num_cols*test_num_rows):
        give_me_a_negative_two += '000011'
    length_of_give_me_a_negative_two = len(give_me_a_negative_two)
    give_me_a_negative_two_as_ints: list[int] = [int(give_me_a_negative_two[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_negative_two)//8)]
    give_me_a_negative_two_as_bytes: bytes = bytes(give_me_a_negative_two_as_ints)
    expected_result: list[int] = [-2] * N

    result: PolyCoefs = _cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)
    assert isinstance(result, PolyCoefs)
    assert result.k1 == test_num_rows
    assert result.k2 == test_num_cols
    assert all(z == -2 for x in result.vals for y in x for z in y)

    result: PolyCoefs = cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)
    assert isinstance(result, PolyCoefs)

    with pytest.raises(TypeError):
        cbd_polycoefs(x='hello world', eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(ValueError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes[:-1], eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta='hello world', num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows='test_num_rows', num_cols=test_num_cols)

    with pytest.raises(ValueError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows=0, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows=test_num_rows, num_cols='test_num_cols')

    with pytest.raises(ValueError):
        cbd_polycoefs(x=give_me_a_negative_two_as_bytes, eta=test_eta, num_rows=test_num_rows, num_cols=0)


COMPRESS_ONE_INT_CASES = [
    (0, 1, 17, 0),
    (1, 1, 17, 0),
    (16, 1, 17, 0),
    (2, 1, 17, 0),
    (15, 1, 17, 0),
    (3, 1, 17, 0),
    (14, 1, 17, 0),
    (4, 1, 17, 0),
    (13, 1, 17, 0),
    (5, 1, 17, 1),
    (12, 1, 17, 1),
    (6, 1, 17, 1),
    (11, 1, 17, 1),
    (7, 1, 17, 1),
    (10, 1, 17, 1),
    (8, 1, 17, 1),
    (9, 1, 17, 1),
]


@pytest.mark.parametrize("x,d,p,expected_result", COMPRESS_ONE_INT_CASES)
def test_compress_one_int(x, d, p, expected_result):
    assert _compress_one_int(x=x, d=d, p=p) == expected_result


def test_compress_list_of_ints():
    # applies compress_one_int many times
    pass


def test_compress_many_ints():
    # applies compress_list_of_ints many times
    pass


def test_compress_polycoefs():
    # returns a PolyCoefs object with compress_many_ints applied to the vals.
    pass


def test_should_compress_many():
    x: list[list[list[int]]] = [[[7]]]
    d: int = 1
    p: int = 17
    assert _should_compress_many(x=x, d=d, p=p)

    x: int = 7
    d: int = 1
    p: int = 17
    assert not _should_compress_many(x=x, d=d, p=p)

    x: list[int] = [7]
    d: int = 1
    p: int = 17
    assert not _should_compress_many(x=x, d=d, p=p)

    d: int = 0
    assert not _should_compress_many(x=x, d=d, p=p)

    d: str = 'hello world'
    assert not _should_compress_many(x=x, d=d, p=p)

    d: int = 1
    p: int = 1
    assert not _should_compress_many(x=x, d=d, p=p)

    p: str = 'hello world'
    assert not _should_compress_many(x=x, d=d, p=p)


COMPRESS_CASES = [
    (0, 1, 17, 0),
    (1, 1, 17, 0),
    (2, 1, 17, 0),
    (3, 1, 17, 0),
    (4, 1, 17, 0),
    (5, 1, 17, 1),
    (6, 1, 17, 1),
    (7, 1, 17, 1),
    (8, 1, 17, 1),
    (9, 1, 17, 1),
    (10, 1, 17, 1),
    (11, 1, 17, 1),
    (12, 1, 17, 1),
    (13, 1, 17, 0),
    (14, 1, 17, 0),
    (15, 1, 17, 0),
    (16, 1, 17, 0),
    ([0], 1, 17, [0]),
    ([1], 1, 17, [0]),
    ([2], 1, 17, [0]),
    ([3], 1, 17, [0]),
    ([4], 1, 17, [0]),
    ([5], 1, 17, [1]),
    ([6], 1, 17, [1]),
    ([7], 1, 17, [1]),
    ([8], 1, 17, [1]),
    ([9], 1, 17, [1]),
    ([10], 1, 17, [1]),
    ([11], 1, 17, [1]),
    ([12], 1, 17, [1]),
    ([13], 1, 17, [0]),
    ([14], 1, 17, [0]),
    ([15], 1, 17, [0]),
    ([16], 1, 17, [0]),
    ([0, 5], 1, 17, [0, 1]),
    ([1, 6], 1, 17, [0, 1]),
    ([2, 7], 1, 17, [0, 1]),
    ([3, 8], 1, 17, [0, 1]),
    ([4, 9], 1, 17, [0, 1]),
    ([5, 10], 1, 17, [1, 1]),
    ([6, 11], 1, 17, [1, 1]),
    ([7, 12], 1, 17, [1, 1]),
    ([8, 13], 1, 17, [1, 0]),
    ([9, 14], 1, 17, [1, 0]),
    ([10, 15], 1, 17, [1, 0]),
    ([11, 16], 1, 17, [1, 0]),
    ([12, 0], 1, 17, [1, 0]),
    ([13, 1], 1, 17, [0, 0]),
    ([14, 2], 1, 17, [0, 0]),
    ([15, 3], 1, 17, [0, 0]),
    ([16, 4], 1, 17, [0, 0]),
    ([[[0]]], 1, 17, [[[0]]]),
    ([[[1]]], 1, 17, [[[0]]]),
    ([[[2]]], 1, 17, [[[0]]]),
    ([[[3]]], 1, 17, [[[0]]]),
    ([[[4]]], 1, 17, [[[0]]]),
    ([[[5]]], 1, 17, [[[1]]]),
    ([[[6]]], 1, 17, [[[1]]]),
    ([[[7]]], 1, 17, [[[1]]]),
    ([[[8]]], 1, 17, [[[1]]]),
    ([[[9]]], 1, 17, [[[1]]]),
    ([[[10]]], 1, 17, [[[1]]]),
    ([[[11]]], 1, 17, [[[1]]]),
    ([[[12]]], 1, 17, [[[1]]]),
    ([[[13]]], 1, 17, [[[0]]]),
    ([[[14]]], 1, 17, [[[0]]]),
    ([[[15]]], 1, 17, [[[0]]]),
    ([[[16]]], 1, 17, [[[0]]]),
    ([[[0, 5]]], 1, 17, [[[0, 1]]]),
    ([[[1, 6]]], 1, 17, [[[0, 1]]]),
    ([[[2, 7]]], 1, 17, [[[0, 1]]]),
    ([[[3, 8]]], 1, 17, [[[0, 1]]]),
    ([[[4, 9]]], 1, 17, [[[0, 1]]]),
    ([[[5, 10]]], 1, 17, [[[1, 1]]]),
    ([[[6, 11]]], 1, 17, [[[1, 1]]]),
    ([[[7, 12]]], 1, 17, [[[1, 1]]]),
    ([[[8, 13]]], 1, 17, [[[1, 0]]]),
    ([[[9, 14]]], 1, 17, [[[1, 0]]]),
    ([[[10, 15]]], 1, 17, [[[1, 0]]]),
    ([[[11, 16]]], 1, 17, [[[1, 0]]]),
    ([[[12, 0]]], 1, 17, [[[1, 0]]]),
    ([[[13, 1]]], 1, 17, [[[0, 0]]]),
    ([[[14, 2]]], 1, 17, [[[0, 0]]]),
    ([[[15, 3]]], 1, 17, [[[0, 0]]]),
    ([[[16, 4]]], 1, 17, [[[0, 0]]]),
]


@pytest.mark.parametrize("x,d,p,expected_result", COMPRESS_CASES)
def test_compress(x, d, p, expected_result):
    assert compress(x=x, d=d, p=p) == expected_result


DECOMPRESS_ONE_INT_CASES = [
    (0, 1, 17, 0),
    (1, 1, 17, 9),
]


@pytest.mark.parametrize("x,d,p,expected_result", DECOMPRESS_ONE_INT_CASES)
def test_decompress_one_int(x, d, p, expected_result):
    assert _decompress_one_int(x=x, d=d, p=p) == expected_result


def test_decompress_list_of_ints():
    pass


def test_decompress_many_ints():
    pass


def test_decompress_polycoefs():
    pass


DECOMPRESS_CASES = [
    (0, 1, 17, 0),
    (1, 1, 17, 9),
    ([0], 1, 17, [0]),
    ([1], 1, 17, [9]),
    ([[[0]]], 1, 17, [[[0]]]),
    ([[[1]]], 1, 17, [[[9]]]),
]


@pytest.mark.parametrize("x,d,p,expected_result", DECOMPRESS_CASES)
def test_decompress(x, d, p, expected_result):
    decompress(x=x, d=d, p=p) == expected_result


ENCODE_M_LIST_OF_INTS_CASES = [
    (list(range(4)), 2, bytes([39])),
    (list(range(10)), 4, bytes([8, 76, 42, 110, 25]))
]


@pytest.mark.parametrize("x,m,expected_result", ENCODE_M_LIST_OF_INTS_CASES)
def test_encode_m_list_of_ints(x, m, expected_result):
    assert _encode_m_list_of_ints_to_bytes(x=x, bits_per_int=m) == expected_result
    assert len(expected_result) == ceil(m*len(x)/8)


ENCODE_M_MANY_CASES = [
    ([[list(range(16))], [list(range(16))]], 4, bytes([8, 76, 42, 110, 25, 93, 59, 127, 8, 76, 42, 110, 25, 93, 59, 127]))
]


@pytest.mark.parametrize("x,m,expected_result", ENCODE_M_MANY_CASES)
def test_encode_m_many(x,m,expected_result):
    assert _encode_m_many(x=x, bits_per_int=m) == expected_result


ENCODE_M_MATRIX_CASES = [
    (PolyNTT(vals=[[list(range(16))], [list(range(16))]], q=17, n=16, k1=2, k2=1), 4, (bytes([8, 76, 42, 110, 25, 93, 59, 127, 8, 76, 42, 110, 25, 93, 59, 127]), 17, 16, 2, 1, False)),
    (PolyCoefs(vals=[[list(range(16))], [list(range(16))]], q=17, n=16, k1=2, k2=1), 4, (bytes([8, 76, 42, 110, 25, 93, 59, 127, 8, 76, 42, 110, 25, 93, 59, 127]), 17, 16, 2, 1, True))
]


@pytest.mark.parametrize("x,m,expected_result", ENCODE_M_MATRIX_CASES)
def test_encode_m_matrix(x, m, expected_result):
    assert _encode_m_matrix(x=x, bits_per_int=m) == expected_result


def test_encode_m():
    pass


def test_decode_m_list_of_ints_from_bytes():
    bits_per_int: int = 10
    some_ints: list[int] = [getrandbits(bits_per_int) for _ in range(2)]
    encoded_ints: bytes = encode_m(x=some_ints, bits_per_int=bits_per_int)
    decoded_encoded_ints: list[int] = _decode_m_list_of_ints_from_bytes(x=encoded_ints, bits_per_int=bits_per_int)
    assert decoded_encoded_ints == some_ints


def test_decode_m_many():
    pass


def test_decode_m_matrix():
    pass


def test_decode_m():
    pass


# def test_ntt_one():
#     pass
#
#
# def test_ntt_many():
#     pass
#
#
# def test_ntt_raw():
#     pass

# We only test the NTT and inverse of constant polynomials here; more thorough tests are advisable
NTT_CASES = [
    (False, [1] + [0] * (N-1), [1] * N),
    (False, [2] + [0] * (N-1), [2] * N),
    (False, [3] + [0] * (N-1), [3] * N),
    (False, [4] + [0] * (N-1), [4] * N),
    (False, [5] + [0] * (N-1), [5] * N),
    (False, [6] + [0] * (N-1), [6] * N),
    (False, [7] + [0] * (N-1), [7] * N),
    (False, [8] + [0] * (N-1), [8] * N),
    (False, [-8] + [0] * (N-1), [-8] * N),
    (False, [-7] + [0] * (N-1), [-7] * N),
    (False, [-6] + [0] * (N-1), [-6] * N),
    (False, [-5] + [0] * (N-1), [-5] * N),
    (False, [-4] + [0] * (N-1), [-4] * N),
    (False, [-3] + [0] * (N-1), [-3] * N),
    (False, [-2] + [0] * (N-1), [-2] * N),
    (False, [-1] + [0] * (N-1), [-1] * N),
    (True, [1] * N, [1] + [0] * (N-1)),
    (True, [2] * N, [2] + [0] * (N-1)),
    (True, [3] * N, [3] + [0] * (N-1)),
    (True, [4] * N, [4] + [0] * (N-1)),
    (True, [5] * N, [5] + [0] * (N-1)),
    (True, [6] * N, [6] + [0] * (N-1)),
    (True, [7] * N, [7] + [0] * (N-1)),
    (True, [8] * N, [8] + [0] * (N-1)),
    (True, [-8] * N, [-8] + [0] * (N-1)),
    (True, [-7] * N, [-7] + [0] * (N-1)),
    (True, [-6] * N, [-6] + [0] * (N-1)),
    (True, [-5] * N, [-5] + [0] * (N-1)),
    (True, [-4] * N, [-4] + [0] * (N-1)),
    (True, [-3] * N, [-3] + [0] * (N-1)),
    (True, [-2] * N, [-2] + [0] * (N-1)),
    (True, [-1] * N, [-1] + [0] * (N-1)),
]


@pytest.mark.parametrize("inv_flag,x,expected_output", NTT_CASES)
def test_ntt_one(inv_flag, x, expected_output):
    assert expected_output == _ntt_one(x=x, inv_flag=inv_flag, const_time=False)
    assert expected_output == _ntt_one(x=x, inv_flag=inv_flag)


def test_transpose():
    pass


def test_xof():
    pass


def test_prf():
    pass


def test_kdf():
    pass


def test_hash_h():
    pass


def test_hash_g():
    pass


def test_cpa_pke_keygen(mocker):
    mocker.patch('crystals.kyber.hash_g', return_value=bytes(range(2 * SEED_LEN_IN_BYTES)))
    mocker.patch('crystals.kyber.parse', return_value=list(range(N)))
    mocker.patch('crystals.kyber.cbd_eta', return_value=list(range(N)))
    some_keys = cpa_pke_keygen()
    assert isinstance(some_keys, bytes)
    assert len(some_keys) == CPA_PKE_PK_LEN + CPA_PKE_SK_LEN


def test_cpa_pke_enc(mocker):
    mocker.patch('crystals.kyber.hash_g', return_value=bytes(list(range(2 * SEED_LEN_IN_BYTES))))
    mocker.patch('crystals.kyber.parse', return_value=list(range(N)))
    mocker.patch('crystals.kyber.cbd_eta', return_value=list(range(N)))

    # Roll a random 32-byte plaintext. Note: all plaintexts must be 32-byte messages.
    plaintext: bytes = getrandbits(SEED_LEN_IN_BYTES * 8).to_bytes(length=SEED_LEN_IN_BYTES, byteorder='big')

    some_keys: bytes = cpa_pke_keygen()
    assert len(some_keys) == CPA_PKE_PK_LEN + CPA_PKE_SK_LEN
    pk: bytes = some_keys[:CPA_PKE_PK_LEN]
    randomness: bytes = getrandbits(SEED_LEN_IN_BYTES * 8).to_bytes(length=SEED_LEN_IN_BYTES, byteorder='big')
    ciphertext: bytes = _cpa_pke_enc(pk=pk, plaintext=plaintext, randomness=randomness)
    assert isinstance(ciphertext, bytes)
    assert len(ciphertext) == CPA_PKE_CIPHERTEXT_LEN


def test_cpa_pka_encrypt():
    pass


def test_cpa_pke_dec(mocker):
    mocker.patch('crystals.kyber.hash_g', return_value=bytes(list(range(2 * SEED_LEN_IN_BYTES))))
    mocker.patch('crystals.kyber.parse', return_value=list(range(N)))
    mocker.patch('crystals.kyber.cbd_eta', return_value=list(range(N)))

    # Roll a random 32-byte plaintext. Note: all plaintexts must be 32-byte messages.
    plaintext: bytes = getrandbits(SEED_LEN_IN_BYTES * 8).to_bytes(length=SEED_LEN_IN_BYTES, byteorder='big')

    some_keys: bytes = cpa_pke_keygen()
    assert len(some_keys) == CPA_PKE_PK_LEN + CPA_PKE_SK_LEN

    pk: bytes = some_keys[:CPA_PKE_PK_LEN]
    randomness: bytes = getrandbits(SEED_LEN_IN_BYTES * 8).to_bytes(length=SEED_LEN_IN_BYTES, byteorder='big')
    ciphertext: bytes = _cpa_pke_enc(pk=pk, plaintext=plaintext, randomness=randomness)
    assert isinstance(ciphertext, bytes)
    assert len(ciphertext) == CPA_PKE_CIPHERTEXT_LEN

    sk: bytes = some_keys[CPA_PKE_PK_LEN:]
    decrypted_ciphertext: bytes = _cpa_pke_dec(sk=sk, ciphertext=ciphertext)
    assert isinstance(decrypted_ciphertext, bytes)
    assert len(decrypted_ciphertext) == len(plaintext)
    assert decrypted_ciphertext == plaintext


def test_cpa_pke_decrypt():
    pass

#
# def test_cca_kem_keygen():
#     pass
#
#
# def test_cca_kem_enc():
#     pass
#
#
# def test_cca_kem_encapsulate():
#     pass
#
#
# def test_cca_kem_dec():
#     pass
#
#
# def test_cca_kem_decapsulate():
#     pass


