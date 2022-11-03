from crystals.kyber import _int2bytes
import pytest
from random import randbytes


SAMPLE_SIZE: int = 2**10


def test_int2bytes():
    with pytest.raises(TypeError):
        _int2bytes(x='hello world', length=17)
        _int2bytes(x=0.001, length=17)
    assert _int2bytes(x=7, length=17) == b'00000000000000111'


# def int2bytes_inverse_of_bytes2int():
#     for length in range(8, 32):
#         for next_int in range(2**length):
#             assert next_int == bytes2int(x=int2bytes(x=next_int, length=length))
#         some_random_bytes: list[bytes] = [randbytes(length) for _ in range(SAMPLE_SIZE)]
#         for next_bytes in some_random_bytes:
#             assert next_bytes == int2bytes(x=bytes2int(x=next_bytes), length=length)
#
#
# def test_encode_decode_inverses():
#     for m in range(2, 17):
#         a_polynomial: list[int] = [randbits(m) for _ in range(N)]
#         encoded_a_polynomial: bytes = encode_m(x=a_polynomial, m=m)
#         decoded_a_polynomial: list[int] = decode_m(x=encoded_a_polynomial, m=m)
#         assert a_polynomial == decoded_a_polynomial
#
#         a_matrix_of_polynomials: list[list[list[int]]] = [[[randbits(m) for _ in range(N)]] for j in range(K)]
#         encoded_a_matrix_of_polynomials: bytes = encode_m(x=a_matrix_of_polynomials, m=m)
#         decoded_a_matrix_of_polynomials: list[list[list[int]]] = decode_m(x=encoded_a_matrix_of_polynomials, m=m)
#         assert a_matrix_of_polynomials == decoded_a_matrix_of_polynomials
#
#
# def test_cpa_pke_keygen():
#     pass
#
#
# def test_cpa_pke_encrypt():
#     pass
#
#
# def test_cpa_pke_decrypt():
#     pass
#
#
# def test_cca_kem_keygen():
#     pass
#
#
# def test_cca_kem_encapsulate():
#     pass
#
#
# def test_cca_kem_decapsulate():
#     pass
#
