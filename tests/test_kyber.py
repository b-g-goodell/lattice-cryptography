from crystals.kyber import SEED_LEN_IN_BYTES, Q, bit_rev, is_pow_two, _bit_rev_cp, bit_rev_cp, reduce, _round_up, round_up, N, LOG_Q, _our_parse_one, _our_parse_many, K, parse, is_arithmetic_legal, PolyCoefs, PolyNTT, _cbd_eta, cbd_eta, _cbd_polycoefs, cbd_polycoefs, _compress_one_int, _decompress_one_int, _should_compress_many, compress, decompress, _encode_m_list_of_ints_to_bytes, _encode_m_many, encode_m, _decode_m_list_of_ints_from_bytes, cpa_pke_keygen, CPA_PKE_SK_LEN, CPA_PKE_PK_LEN, CPA_PKE_CIPHERTEXT_LEN, _cpa_pke_enc, _encode_m_matrix, _cpa_pke_dec, _ntt_one, ntt, make_zetas_and_invs
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
    assert is_arithmetic_legal(a_vals=[[[]]], b_vals=[[[]]], q_a=17, q_b=17)
    assert not is_arithmetic_legal(a_vals=[[[1]]], b_vals=[[[]]], q_a=17, q_b=17)
    assert not is_arithmetic_legal(a_vals=[[[1]]], b_vals=[[[1]]], q_a=17, q_b=19)


def test_add():
    # too simple to bother testing
    pass


def test_against_nayuki():
    # https://www.nayuki.io/page/number-theoretic-transform-integer-dft
    q = 17
    n = 8
    lgn = 3
    rou = 2
    while rou < q:
        if all((rou ** i) % q != 1 for i in range(1, n)) and (rou ** n) % q == 1:
            break
        rou += 1
    rou_inverse = (rou ** (n-1)) % q
    powers = [n // (2**(s+1)) for s in range(3)]
    zetas = [(rou ** i) % q for i in powers]
    zetas = [i if i <= q else i - q for i in zetas]
    zeta_inverses = [(rou_inverse ** i) % q for i in powers]
    zeta_inverses = [i if i <= q else i - q for i in zeta_inverses]

    input_vector = [6, 0, 10, 7, 2, 0, 0, 0]
    observed_ntt = _ntt_one(x=input_vector, inv_flag=False, q=q, n=n, log_n=lgn, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses)
    expected_ntt = [8, 15, 4, 12, 11, 5, 9, 1]
    assert all((i-j) % q == 0 for i, j in zip(observed_ntt, expected_ntt))

    observed_intt_of_ntt = _ntt_one(x=observed_ntt, inv_flag=True, q=q, n=n, log_n=lgn, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses)
    assert all((i-j) % q == 0 for i, j in zip(observed_intt_of_ntt, input_vector))


def test_identities():
    q = 17
    n = 2
    log_n = 1
    k1 = 1
    k2 = 1
    zetas: list[int]
    zeta_inverses: list[int]
    zetas, zeta_inverses = make_zetas_and_invs(q=q, n=n, lgn=log_n)
    identity_matrix_coefs: list[list[list[int]]] = [[[0] * n] * k2] * k1
    identity_matrix_coefs[0][0][0] = 1
    identity_matrix: PolyCoefs = PolyCoefs(vals=identity_matrix_coefs, q=q, n=n, k1=k1, k2=k1, zetas=zetas, zeta_inverses=zeta_inverses)
    identity_matrix_ntt: PolyNTT = ntt(x=identity_matrix, q=q, n=n, log_n=log_n, half_q=q//2)
    assert len(identity_matrix_ntt.vals) == k1
    assert len(identity_matrix_ntt.vals[0]) == k2
    assert len(identity_matrix_ntt.vals[0][0]) == n
    assert all(z == 1 for x in identity_matrix_ntt.vals for y in x for z in y)

    new_random_matrix_coefs: list[list[list[int]]] = [[[randrange(q) for l in range(n)] for j in range(k2)] for i in range(k1)]
    new_random_matrix: PolyCoefs = PolyCoefs(vals=new_random_matrix_coefs, q=q, n=n, k1=k1, k2=k2, zetas=zetas, zeta_inverses=zeta_inverses)
    new_random_matrix_ntt: PolyNTT = ntt(x=new_random_matrix, q=q, n=n, log_n=log_n, half_q=q//2)
    assert len(new_random_matrix_ntt.vals) == k1
    assert len(new_random_matrix_ntt.vals[0]) == k2
    assert len(new_random_matrix_ntt.vals[0][0]) == n

    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix
    with pytest.raises(NotImplementedError):
        new_random_matrix * identity_matrix_ntt
    with pytest.raises(NotImplementedError):
        new_random_matrix_ntt * identity_matrix

    id_times_new_ntt: PolyNTT = identity_matrix_ntt * new_random_matrix_ntt
    new_times_id_ntt: PolyNTT = new_random_matrix_ntt * identity_matrix_ntt
    id_times_new: PolyCoefs = ntt(id_times_new_ntt, q=q, n=n, log_n=log_n, half_q=q//2)
    new_times_id: PolyCoefs = ntt(new_times_id_ntt, q=q, n=n, log_n=log_n, half_q=q//2)

    assert id_times_new.q == new_times_id.q == new_random_matrix.q == q
    assert id_times_new.k1 == new_times_id.k1 == new_random_matrix.k1 == k1
    assert id_times_new.k2 == new_times_id.k2 == new_random_matrix.k2 == k2
    for i in range(k1):
        for j in range(k2):
            for l in range(n):
                assert (id_times_new.vals[i][j][l] - new_random_matrix.vals[i][j][l]) % q == 0
    assert id_times_new == new_random_matrix
    assert new_times_id == new_random_matrix


def test_poly_mul():
    # Polynomials from Z_q[X]/(X**4 + 1)
    q = 17
    half_q = q//2
    n = 4
    log_n = 2
    zetas, zeta_inverses = make_zetas_and_invs(q=q, n=2*n, lgn=log_n+1)
    rou = zetas[-1]

    a = [randrange(q) for _ in range(n)] + [0 for _ in range(n)]  # a[0] + a[1]*X + a[2]*X**2 + a[3]*X**3
    b = [randrange(q) for _ in range(n)] + [0 for _ in range(n)]  # b[0] + b[1]*X + b[2]*X**2 + b[3]*X**3
    ab_full = [a[0]*b[0],
               a[0]*b[1] + a[1]*b[0],
               a[0]*b[2] + a[1]*b[1] + a[2]*b[0],
               a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0],
               a[1]*b[3] + a[2]*b[2] + a[3]*b[1],
               a[2]*b[3] + a[3]*b[2],
               a[3]*b[3]]  # a[0]*b[0] + (a[0]*b[1] + a[1]*b[0])*X + (a[0]*b[2] + a[1]*b[1] + a[2]*b[0])*X**2 + (a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0])*X**3 + (a[1]*b[3] + a[2]*b[2] + a[3]*b[1])*X**4 + (a[2]*b[3] + a[3]*b[2])*X**5 + a[3]*b[3]*X**6
    ab_full += [0 for _ in range(2*n - len(ab_full))]

    a_hat = _ntt_one(x=a, inv_flag=False, const_time=False, q=q, n=n, log_n=log_n, half_q=half_q, zetas=zetas, zeta_inverses=zeta_inverses)
    b_hat = _ntt_one(x=b, inv_flag=False, const_time=False, q=q, n=n, log_n=log_n, half_q=half_q, zetas=zetas, zeta_inverses=zeta_inverses)
    ab_hat = _ntt_one(x=ab_full, inv_flag=False, const_time=False, q=q, n=n, log_n=log_n, half_q=half_q, zetas=zetas, zeta_inverses=zeta_inverses)

    coordinate_wise_a_hat_times_b_hat = [(x * y) % q if (x*y) % q <= half_q else ((x*y) % q) - q for x, y in zip(a_hat, b_hat)]
    intt_of_coordinate_wise = _ntt_one(x=coordinate_wise_a_hat_times_b_hat, inv_flag=True, const_time=False, q=q, n=n, log_n=log_n, half_q=half_q, zetas=zetas, zeta_inverses=zeta_inverses)


    assert all((x-y) % q == 0 for x, y in zip(ab_hat, coordinate_wise_a_hat_times_b_hat))

def test_matrix_mul():
    q = 17
    n = 2
    log_n = 1
    zetas, zeta_inverses = make_zetas_and_invs(q=q, n=n, lgn=log_n)

    left_k1 = 3
    left_k2 = 2
    left_matrix_coefs: list[list[list[int]]] = [[[0 for k in range(n)] for j in range(left_k2)] for i in range(left_k1)]  # 3x2 matrix
    for i in range(left_k1):
        for j in range(left_k2):
            for k in range(n):
                left_matrix_coefs[i][j][k] = randrange(q)
    left_matrix: PolyCoefs = PolyCoefs(vals=left_matrix_coefs, q=q, n=n, k1=left_k1, k2=left_k2)  # 3x2 matrix

    a_poly = left_matrix_coefs[0][0]
    b_poly = left_matrix_coefs[0][1]
    c_poly = left_matrix_coefs[1][0]
    d_poly = left_matrix_coefs[1][1]
    e_poly = left_matrix_coefs[2][0]
    f_poly = left_matrix_coefs[2][1]

    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[0][0], a_poly))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[0][1], b_poly))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[1][0], c_poly))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[1][1], d_poly))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[2][0], e_poly))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix.vals[2][1], f_poly))

    a_poly_ntt = _ntt_one(x=a_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    b_poly_ntt = _ntt_one(x=b_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    c_poly_ntt = _ntt_one(x=c_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    d_poly_ntt = _ntt_one(x=d_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    e_poly_ntt = _ntt_one(x=e_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    f_poly_ntt = _ntt_one(x=f_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)

    a_poly_intt_of_ntt = _ntt_one(x=a_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    b_poly_intt_of_ntt = _ntt_one(x=b_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    c_poly_intt_of_ntt = _ntt_one(x=c_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    d_poly_intt_of_ntt = _ntt_one(x=d_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    e_poly_intt_of_ntt = _ntt_one(x=e_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    f_poly_intt_of_ntt = _ntt_one(x=f_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)

    assert all((x-y) % q == 0 for x, y in zip(a_poly_intt_of_ntt, a_poly))
    assert all((x-y) % q == 0 for x, y in zip(b_poly_intt_of_ntt, b_poly))
    assert all((x-y) % q == 0 for x, y in zip(c_poly_intt_of_ntt, c_poly))
    assert all((x-y) % q == 0 for x, y in zip(d_poly_intt_of_ntt, d_poly))
    assert all((x-y) % q == 0 for x, y in zip(e_poly_intt_of_ntt, e_poly))
    assert all((x-y) % q == 0 for x, y in zip(f_poly_intt_of_ntt, f_poly))

    left_matrix_ntt: PolyNTT = ntt(x=left_matrix, q=q, n=n, log_n=log_n, half_q=q//2)

    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[0][0], a_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[0][1], b_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[1][0], c_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[1][1], d_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[2][0], e_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(left_matrix_ntt.vals[2][1], f_poly_ntt))

    left_matrix_ntt_manual: PolyNTT = PolyNTT(vals=[[a_poly_ntt, b_poly_ntt], [c_poly_ntt, d_poly_ntt], [e_poly_ntt, f_poly_ntt]], q=q, n=n, k1=left_k1, k2=left_k2, zetas=zetas, zeta_inverses=zeta_inverses)
    assert left_matrix_ntt == left_matrix_ntt_manual

    left_matrix_intt_of_ntt: PolyCoefs = ntt(x=left_matrix_ntt, q=q, n=n, log_n=log_n, half_q=q//2)
    assert left_matrix_intt_of_ntt == left_matrix

    right_k1 = 2
    right_k2 = 2
    right_matrix_coefs: list[list[list[int]]] = [[[0 for k in range(n)] for j in range(right_k2)] for i in range(right_k1)]  # 3x2 matrix
    for i in range(right_k1):
        for j in range(right_k2):
            for k in range(n):
                right_matrix_coefs[i][j][k] = randrange(q)
    right_matrix: PolyCoefs = PolyCoefs(vals=right_matrix_coefs, q=q, n=n, k1=right_k1, k2=right_k2)  # 3x2 matrix

    g_poly = right_matrix_coefs[0][0]
    h_poly = right_matrix_coefs[0][1]
    i_poly = right_matrix_coefs[1][0]
    j_poly = right_matrix_coefs[1][1]

    assert all((x-y) % q == 0 for x, y in zip(right_matrix.vals[0][0], g_poly))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix.vals[0][1], h_poly))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix.vals[1][0], i_poly))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix.vals[1][1], j_poly))

    g_poly_ntt = _ntt_one(x=g_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    h_poly_ntt = _ntt_one(x=h_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    i_poly_ntt = _ntt_one(x=i_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)
    j_poly_ntt = _ntt_one(x=j_poly, q=q, n=n, log_n=log_n, half_q=q//2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=False)

    g_poly_intt_of_ntt = _ntt_one(x=g_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    h_poly_intt_of_ntt = _ntt_one(x=h_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    i_poly_intt_of_ntt = _ntt_one(x=i_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)
    j_poly_intt_of_ntt = _ntt_one(x=j_poly_ntt, q=q, n=n, log_n=log_n, half_q=q // 2, zetas=zetas, zeta_inverses=zeta_inverses, inv_flag=True)

    assert all((x-y) % q == 0 for x, y in zip(g_poly_intt_of_ntt, g_poly))
    assert all((x-y) % q == 0 for x, y in zip(h_poly_intt_of_ntt, h_poly))
    assert all((x-y) % q == 0 for x, y in zip(i_poly_intt_of_ntt, i_poly))
    assert all((x-y) % q == 0 for x, y in zip(j_poly_intt_of_ntt, j_poly))

    right_matrix_ntt: PolyNTT = ntt(x=right_matrix, q=q, n=n, log_n=log_n, half_q=q//2)

    assert all((x-y) % q == 0 for x, y in zip(right_matrix_ntt.vals[0][0], g_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix_ntt.vals[0][1], h_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix_ntt.vals[1][0], i_poly_ntt))
    assert all((x-y) % q == 0 for x, y in zip(right_matrix_ntt.vals[1][1], j_poly_ntt))

    right_matrix_ntt_manual: PolyNTT = PolyNTT(vals=[[g_poly_ntt, h_poly_ntt], [i_poly_ntt, j_poly_ntt]], q=q, n=n, k1=right_k1, k2=right_k2, zetas=zetas, zeta_inverses=zeta_inverses)
    assert right_matrix_ntt == right_matrix_ntt_manual

    right_matrix_intt_of_ntt: PolyCoefs = ntt(x=right_matrix_ntt, q=q, n=n, log_n=log_n, half_q=q//2)
    assert right_matrix_intt_of_ntt == right_matrix

    expected_n: int = left_matrix.n
    expected_k1: int = left_matrix.k1
    expected_k2: int = right_matrix.k2
    inner_k: int = left_matrix.k2
    expected_product_coefs: list[list[list[int]]] = [[[0, 0] for j in range(right_matrix.k2)] for i in range(left_matrix.k1)]

    for i in range(expected_k1):
        for j in range(expected_k2):
            expected_product_coefs[i][j] = [
                sum(left_matrix.vals[i][l][0] * right_matrix.vals[l][j][0] - left_matrix.vals[i][l][1] * right_matrix.vals[l][j][1] for l in range(inner_k)) % q,
                sum(left_matrix.vals[i][l][1] * right_matrix.vals[l][j][0] + left_matrix.vals[i][l][0] * right_matrix.vals[l][j][1] for l in range(inner_k)) % q
            ]
    expected_product: PolyCoefs = PolyCoefs(vals=expected_product_coefs, q=q, n=n, k1=expected_k1, k2=expected_k2)
    expected_product_ntt: PolyNTT = ntt(x=expected_product, q=q, n=n, log_n=log_n, half_q=q//2)

    left_matrix_ntt_times_right_matrix_ntt: PolyNTT = left_matrix_ntt * right_matrix_ntt
    intt_of_left_matrix_ntt_times_right_matrix_ntt: PolyCoefs = ntt(x=left_matrix_ntt_times_right_matrix_ntt, q=q, n=n, log_n=log_n, half_q=q//2)
    assert intt_of_left_matrix_ntt_times_right_matrix_ntt.k1 == expected_product.k1 == expected_k1
    assert intt_of_left_matrix_ntt_times_right_matrix_ntt.k2 == expected_product.k2 == expected_k2
    assert intt_of_left_matrix_ntt_times_right_matrix_ntt.n == expected_product.n == n
    assert intt_of_left_matrix_ntt_times_right_matrix_ntt.q == expected_product.q == q
    for row_a, row_b in zip(intt_of_left_matrix_ntt_times_right_matrix_ntt.vals, expected_product.vals):
        for col_a, col_b in zip(row_a, row_b):
            for coef_a, coef_b in zip(col_a, col_b):
                assert (coef_a - coef_b) % expected_product.q == 0
    assert intt_of_left_matrix_ntt_times_right_matrix_ntt == expected_product
    assert expected_product_ntt == left_matrix_ntt_times_right_matrix_ntt
    assert left_matrix_ntt_times_right_matrix_ntt.k1 == expected_k1
    assert left_matrix_ntt_times_right_matrix_ntt.k2 == expected_k2
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[0][0][0] - (a_poly_ntt[0] * g_poly_ntt[0] + b_poly_ntt[0] * i_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[0][0][1] - (a_poly_ntt[1] * g_poly_ntt[1] + b_poly_ntt[1] * i_poly_ntt[1])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[0][1][0] - (a_poly_ntt[0] * h_poly_ntt[0] + b_poly_ntt[0] * j_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[0][1][1] - (a_poly_ntt[1] * h_poly_ntt[1] + b_poly_ntt[1] * j_poly_ntt[1])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[1][0][0] - (c_poly_ntt[0] * g_poly_ntt[0] + d_poly_ntt[0] * i_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[1][0][1] - (c_poly_ntt[1] * g_poly_ntt[1] + d_poly_ntt[1] * i_poly_ntt[1])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[1][1][0] - (c_poly_ntt[0] * h_poly_ntt[0] + d_poly_ntt[0] * j_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[1][1][1] - (c_poly_ntt[1] * h_poly_ntt[1] + d_poly_ntt[1] * j_poly_ntt[1])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[2][0][0] - (e_poly_ntt[0] * g_poly_ntt[0] + f_poly_ntt[0] * i_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[2][0][1] - (e_poly_ntt[1] * g_poly_ntt[1] + f_poly_ntt[1] * i_poly_ntt[1])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[2][1][0] - (e_poly_ntt[0] * h_poly_ntt[0] + f_poly_ntt[0] * j_poly_ntt[0])) % q == 0
    assert (left_matrix_ntt_times_right_matrix_ntt.vals[2][1][1] - (e_poly_ntt[1] * h_poly_ntt[1] + f_poly_ntt[1] * j_poly_ntt[1])) % q == 0

    intt_of_left_matrix_ntt_times_right_matrix_ntt: PolyCoefs = ntt(x=left_matrix_ntt_times_right_matrix_ntt, q=q, n=n, log_n=log_n, half_q=q//2)

    expected_coefs: list[list[list[int]]] = [[[0 for k in range(n)] for j in range(expected_k2)] for i in range(expected_k1)]  # 3x2 matrix
    expected_coefs[0][0][0] = ((a_poly[0] * g_poly[0] - a_poly[1] * g_poly[1]) + (b_poly[0] * i_poly[0] - b_poly[1] * i_poly[1])) % q
    expected_coefs[0][0][1] = ((a_poly[1] * g_poly[0] + a_poly[0] * g_poly[1]) + (b_poly[1] * i_poly[0] + b_poly[0] * i_poly[1])) % q
    expected_coefs[0][1][0] = ((a_poly[0] * h_poly[0] - a_poly[1] * h_poly[1]) + (b_poly[0] * j_poly[0] - b_poly[1] * j_poly[1])) % q
    expected_coefs[0][1][1] = ((a_poly[1] * h_poly[0] + a_poly[0] * h_poly[1]) + (b_poly[1] * j_poly[0] + b_poly[0] * j_poly[1])) % q
    expected_coefs[1][0][0] = ((c_poly[0] * g_poly[0] - c_poly[1] * g_poly[1]) + (d_poly[0] * i_poly[0] - d_poly[1] * i_poly[1])) % q
    expected_coefs[1][0][1] = ((c_poly[1] * g_poly[0] + c_poly[0] * g_poly[1]) + (d_poly[1] * i_poly[0] + d_poly[0] * i_poly[1])) % q
    expected_coefs[1][1][0] = ((c_poly[0] * h_poly[0] - c_poly[1] * h_poly[1]) + (d_poly[0] * j_poly[0] - d_poly[1] * j_poly[1])) % q
    expected_coefs[1][1][1] = ((c_poly[1] * h_poly[0] + c_poly[0] * h_poly[1]) + (d_poly[1] * j_poly[0] + d_poly[0] * j_poly[1])) % q
    expected_coefs[2][0][0] = ((e_poly[0] * g_poly[0] - e_poly[1] * g_poly[1]) + (f_poly[0] * i_poly[0] - f_poly[1] * i_poly[1])) % q
    expected_coefs[2][0][1] = ((e_poly[1] * g_poly[0] + e_poly[0] * g_poly[1]) + (f_poly[1] * i_poly[0] + f_poly[0] * i_poly[1])) % q
    expected_coefs[2][1][0] = ((e_poly[0] * h_poly[0] - e_poly[1] * h_poly[1]) + (f_poly[0] * j_poly[0] - f_poly[1] * j_poly[1])) % q
    expected_coefs[2][1][1] = ((e_poly[1] * h_poly[0] + e_poly[0] * h_poly[1]) + (f_poly[1] * j_poly[0] + f_poly[0] * j_poly[1])) % q

    expected: PolyCoefs = PolyCoefs(vals=expected_coefs, q=q, n=n, k1=expected_k1, k2=expected_k2)

    assert all(all((x - y) % q == 0 for x, y in zip(expected.vals[i][j], intt_of_left_matrix_ntt_times_right_matrix_ntt.vals[i][j])) for i in range(expected_k1) for j in range(expected_k2))

    assert intt_of_left_matrix_ntt_times_right_matrix_ntt == expected


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
    assert z.vals == [[[15, 7], [4, 15]], [[10, 6], [16, 14]], [[5, 5], [11, 13]]]

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
    expected_vals = [[[6, 14]], [[5, 4]], [[4, 11]]]
    for a, b in zip(z.vals, expected_vals):
        for c, d in zip(a, b):
            for e, f in zip(c, d):
                assert (e-f) % x.q == 0
    # assert z.vals == [[[6, -3]], [[5, 4]], [[4, -6]]]

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
    give_me_a_negative_two_as_ints: list[int] = [int(give_me_a_negative_two[_*8: (_+1)*8], 2) for _ in range(len(give_me_a_negative_two)//8)]
    give_me_a_negative_two_as_bytes: bytes = bytes(give_me_a_negative_two_as_ints)

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
    assert decompress(x=x, d=d, p=p) == expected_result


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
def test_encode_m_many(x, m, expected_result):
    assert _encode_m_many(x=x, bits_per_int=m) == expected_result


zetas, zeta_inverses = make_zetas_and_invs(q=17, n=16, lgn=4)
ENCODE_M_MATRIX_CASES = [
    (PolyNTT(vals=[[list(range(16))], [list(range(16))]], q=17, n=16, k1=2, k2=1, zetas=zetas, zeta_inverses=zeta_inverses), 4, (bytes([8, 76, 42, 110, 25, 93, 59, 127, 8, 76, 42, 110, 25, 93, 59, 127]), 17, 16, 2, 1, False)),
    (PolyCoefs(vals=[[list(range(16))], [list(range(16))]], q=17, n=16, k1=2, k2=1, zetas=zetas, zeta_inverses=zeta_inverses), 4, (bytes([8, 76, 42, 110, 25, 93, 59, 127, 8, 76, 42, 110, 25, 93, 59, 127]), 17, 16, 2, 1, True))
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

