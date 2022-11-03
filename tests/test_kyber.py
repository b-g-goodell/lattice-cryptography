from crystals.kyber import _int2bytes, int2bytes, _bytes2int, bytes2int, _bit_rev
from random import getrandbits
import pytest
from random import randbytes


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
        print()
