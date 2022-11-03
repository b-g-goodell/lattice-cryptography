from crystals.kyber import _int2bytes, int2bytes, _bytes2int, bytes2int, _bit_rev, bit_rev, is_pow_two, _bit_rev_cp, bit_rev_cp
from random import getrandbits
import pytest
from math import ceil, log2


SAMPLE_SIZE: int = 2**10
LOG_SAMPLE_SIZE: int = 10


def test_int2bytes():
    with pytest.raises(TypeError):
        _int2bytes(x='hello world', length=17)
    with pytest.raises(TypeError):
        _int2bytes(x=0.001, length=17)

    assert _int2bytes(x=7, length=17) == b'00000000000000111'

    with pytest.raises(TypeError):
        int2bytes(x='hello world', length=17)
    with pytest.raises(TypeError):
        int2bytes(x='hello world', length=0.001)
    with pytest.raises(TypeError):
        int2bytes(x=7, length=0.001)
    with pytest.raises(ValueError):
        int2bytes(x=7, length=2)

    assert int2bytes(x=7, length=3) == b'111'
    assert int2bytes(x=7, length=4) == b'0111'
    assert int2bytes(x=7, length=17) == b'00000000000000111'


def test_bytes2int():
    with pytest.raises(ValueError):
        _bytes2int(x='hello world')

    with pytest.raises(TypeError):
        _bytes2int(x=0.001)

    assert _bytes2int(x=b'00000000000000111') == 7

    with pytest.raises(TypeError):
        bytes2int(x='hello world')

    with pytest.raises(TypeError):
        bytes2int(x=0.001)

    assert bytes2int(x=b'00000000000000111') == 7


def test_bytes2int_and_int2bytes_are_inverses():
    for i in range(SAMPLE_SIZE):
        assert bytes2int(int2bytes(x=i, length=LOG_SAMPLE_SIZE)) == i

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
