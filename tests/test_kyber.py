from crystals.kyber import Q, _int2bytes, int2bytes, _bytes2int, bytes2int, _bit_rev, bit_rev, is_pow_two, _bit_rev_cp, bit_rev_cp, _reduce, reduce, _round_up, round_up, N, LOG_Q, _parse_one, _parse_many, K, parse, is_arithmetic_legal, PolyCoefs, PolyNTT, _cbd_eta, cbd_eta, _cbd_polycoefs, cbd_polycoefs, _compress_one_int, _decompress_one_int, _encode_m_one_int
from random import getrandbits, randrange
import pytest
from math import ceil, log2


SAMPLE_SIZE: int = 2**10
LOG_SAMPLE_SIZE: int = 10


def test_int2bytes():
    assert _int2bytes(x=7, length=1) == bytes([7])

    with pytest.raises(TypeError):
        int2bytes(x='hello world', length=17)
    with pytest.raises(TypeError):
        int2bytes(x='hello world', length=0.001)
    with pytest.raises(TypeError):
        int2bytes(x=7, length=0.001)
    with pytest.raises(ValueError):
        int2bytes(x=-1, length=17)
    with pytest.raises(ValueError):
        int2bytes(x=7, length=0)

    assert int2bytes(x=7, length=3) == bytes([0, 0, 7])
    assert int2bytes(x=7, length=4) == bytes([0, 0, 0, 7])
    assert int2bytes(x=7, length=17) == bytes([0] * 16 + [7])


def test_bytes2int():
    with pytest.raises(TypeError):
        bytes2int(x='hello world')

    with pytest.raises(TypeError):
        bytes2int(x=0.001)

    assert _bytes2int(x=bytes([7])) == 7
    assert _bytes2int(x=bytes([0, 0, 7])) == 7
    assert _bytes2int(x=bytes([0, 0, 0, 7])) == 7
    assert _bytes2int(x=bytes([0] * 16 + [7])) == 7

    assert bytes2int(x=bytes([7])) == 7
    assert bytes2int(x=bytes([0, 0, 7])) == 7
    assert bytes2int(x=bytes([0, 0, 0, 7])) == 7
    assert bytes2int(x=bytes([0] * 16 + [7])) == 7


def test_bytes2int_and_int2bytes_are_inverses():
    length: int = LOG_SAMPLE_SIZE // 8 + 1
    for i in range(SAMPLE_SIZE):
        i2bytes = int2bytes(x=i, length=length)
        i2bytes2i = bytes2int(x=i2bytes)
        assert i2bytes2i == i

    for i in range(SAMPLE_SIZE):
        next_bytes_object = bytes(0)
        for j in range(LOG_SAMPLE_SIZE):
            next_bit = getrandbits(1)
            if next_bit:
                next_bytes_object += b'1'
            else:
                next_bytes_object += b'0'
        assert int2bytes(x=bytes2int(x=next_bytes_object), length=LOG_SAMPLE_SIZE) == next_bytes_object


def test_bit_rev():
    for i in range(SAMPLE_SIZE):
        x = int2bytes(x=i, length=LOG_SAMPLE_SIZE)
        reversed_i = _bit_rev(x=i, length=LOG_SAMPLE_SIZE)
        reversed_x = int2bytes(x=reversed_i, length=LOG_SAMPLE_SIZE)
        assert len(x) == len(reversed_x) == LOG_SAMPLE_SIZE
        for i in range(LOG_SAMPLE_SIZE):
            assert x[i] == reversed_x[LOG_SAMPLE_SIZE - 1 - i]

    with pytest.raises(TypeError):
        bit_rev(x='hello world', length=LOG_SAMPLE_SIZE)

    with pytest.raises(TypeError):
        bit_rev(x='hello world', length=0.001)

    with pytest.raises(TypeError):
        bit_rev(x='hello world', length='hello world')

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length=LOG_SAMPLE_SIZE)

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length=0.001)

    with pytest.raises(TypeError):
        bit_rev(x=0.001, length='hello world')

    with pytest.raises(ValueError):
        bit_rev(x=1, length=0)

    with pytest.raises(ValueError):
        bit_rev(x=7, length=2)

    for i in range(SAMPLE_SIZE):
        x = int2bytes(x=i, length=LOG_SAMPLE_SIZE)
        reversed_i = bit_rev(x=i, length=LOG_SAMPLE_SIZE)
        reversed_x = int2bytes(x=reversed_i, length=LOG_SAMPLE_SIZE)
        assert len(x) == len(reversed_x) == LOG_SAMPLE_SIZE
        for i in range(LOG_SAMPLE_SIZE):
            assert x[i] == reversed_x[LOG_SAMPLE_SIZE - 1 - i]


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
    assert _bit_rev_cp(x=x, num_bits=ceil(log2(len(x)))) == expected_output


def test_bit_rev_cp_failures():
    with pytest.raises(TypeError):
        bit_rev_cp(x=3, num_bits=2)

    with pytest.raises(TypeError):
        bit_rev_cp(x=['hello world', 0.001], num_bits=2)

    with pytest.raises(TypeError):
        bit_rev_cp(x=['hello world', 0.001], num_bits=0.001)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), num_bits=0)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), num_bits=2)

    with pytest.raises(ValueError):
        bit_rev_cp(x=list(range(8)), num_bits=4)

@pytest.mark.parametrize("x,expected_output", BIT_REV_CP_CASES)
def test_bit_rev_cp_full(x, expected_output):
    assert bit_rev_cp(x=x, num_bits=ceil(log2(len(x)))) == expected_output


REDUCE_CASES = [(_, _) if _ <= Q//2 else (_, _ - Q) for _ in list(range(Q))]

@pytest.mark.parametrize("x, expected_output", REDUCE_CASES)
def test_reduce(x, expected_output):
    assert _reduce(x=x) == expected_output
    assert reduce(x=x) == expected_output


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
    resulting_ints, index = _parse_one(x=bytes(test_integers_to_bytes))
    expected_ints = [_ for _ in test_integers if _ < Q]
    expected_ints = expected_ints[:N]
    assert len(resulting_ints) == N
    assert all(isinstance(x, int) for x in resulting_ints)
    assert all(0 <= x < Q for x in resulting_ints)
    for i, (next_result_int, next_test_int) in enumerate(zip(resulting_ints, expected_ints)):
        assert next_result_int == next_test_int

    with pytest.raises(RuntimeError):
        _parse_one(x=b'hello world')


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

    resulting_ints: list[list[list[int]]] = _parse_many(x=bytes(test_integers_to_bytes))
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
        _parse_many(x=b'hello world')


def test_parse():
    with pytest.raises(TypeError):
        parse(x='hello world')
    with pytest.raises(ValueError):
        parse(x=b'hello world')


def test_is_arithmetic_legal():
    assert is_arithmetic_legal(a_vals = [[[]]], b_vals=[[[]]])
    assert not is_arithmetic_legal(a_vals=[[[1]]], b_vals=[[[]]])


def test_add():
    # too simple to bother testing
    pass


def test_mul():
    # too simple to bother testing
    pass


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
    test_eta: int = 3
    give_me_a_three: bytes = b'111000'
    test_val: bytes = bytes(0)
    for _ in range(N):
        test_val += give_me_a_three
    expected_result: int = [3 for _ in range(N)]
    assert _cbd_eta(x=test_val, eta=test_eta) == expected_result
    assert cbd_eta(x=test_val, eta=test_eta) == expected_result

    give_me_a_negative_three: bytes = b'000111'
    test_val: bytes = bytes(0)
    for _ in range(N):
        test_val += give_me_a_negative_three
    expected_result: int = [-3 for _ in range(N)]
    assert _cbd_eta(x=test_val, eta=test_eta) == expected_result
    assert cbd_eta(x=test_val, eta=test_eta) == expected_result


def test_cbd_polycoefs():
    test_eta: int = 3
    test_num_rows: int = 2
    test_num_cols: int = 3
    give_me_a_three: bytes = b'111000'
    test_val: bytes = bytes(0)
    for _ in range(test_num_rows*test_num_cols*2*N*test_eta):
        test_val += give_me_a_three

    result = _cbd_polycoefs(x=test_val, eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)
    assert isinstance(result, PolyCoefs)
    assert result.k1 == test_num_rows
    assert result.k2 == test_num_cols
    assert all(z == 3 for x in result.vals for y in x for z in y)

    result: PolyCoefs = cbd_polycoefs(x=test_val, eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)
    assert isinstance(result, PolyCoefs)

    with pytest.raises(TypeError):
        cbd_polycoefs(x='hello world', eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(ValueError):
        cbd_polycoefs(x=give_me_a_three, eta=test_eta, num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=test_val, eta='hello world', num_rows=test_num_rows, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=test_val, eta=test_eta, num_rows='test_num_rows', num_cols=test_num_cols)

    with pytest.raises(ValueError):
        cbd_polycoefs(x=test_val, eta=test_eta, num_rows=0, num_cols=test_num_cols)

    with pytest.raises(TypeError):
        cbd_polycoefs(x=test_val, eta=test_eta, num_rows=test_num_rows, num_cols='test_num_cols')

    with pytest.raises(ValueError):
        cbd_polycoefs(x=test_val, eta=test_eta, num_rows=test_num_rows, num_cols=0)


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
    assert _compress_one_int(x=x,d=d,p=p) == expected_result


DECOMPRESS_ONE_INT_CASES = [
    (0, 1, 17, 0),
    (1, 1, 17, 9),
]

@pytest.mark.parametrize("x,d,p,expected_result", DECOMPRESS_ONE_INT_CASES)
def test_decompress_one_int(x, d, p, expected_result):
    assert _decompress_one_int(x=x, d=d, p=p) == expected_result


ENCODE_M_ONE_INT_CASES = [
    (i, j, int2bytes(x=i, length=j)) for j in range(1, 3) for i in range(2**j)
]

@pytest.mark.parametrize("x,m,expected_result", ENCODE_M_ONE_INT_CASES)
def test_encode_m_one_int(x, m, expected_result):
    assert _encode_m_one_int(x=x, m=m) == expected_result
