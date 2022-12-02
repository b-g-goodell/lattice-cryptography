# Handles generalized NTT
from math import ceil, sqrt, log2
from copy import deepcopy
from warnings import warn


def is_prime(x: int) -> bool:
    return all([x % i != 0 for i in range(2, ceil(sqrt(x)) + 1)])


def _reduce(x: int, q: int) -> int:
    y: int = x % q
    z: int = y - q//2 - 1
    w: int = y - (1 + (z >> ceil(log2(q)))) * q
    return w


def reduce(x: int, q: int) -> int:
    """
    Compute some integer w such that -Q//2 <= w <= Q//2 and such that (x-w) % Q == 0, in constant
    time.
    :param x: Input integer
    :type x: int
    :param q: Input modulus
    :type q: int
    :return: "Centered" representative of x % Q.
    :rtype: int
    """
    if isinstance(x, int) and isinstance(q, int) and q >= 2:
        return _reduce(x=x, q=q)
    elif not isinstance(x, int):
        raise TypeError(f'Cannot compute reduce for x unless x is an integer, but had type(x)={type(x)}.')
    elif not isinstance(q, int):
        raise TypeError(f'Cannot compute reduce with q unless q is an integer, but had type(q)={type(q)}.')
    elif q < 2:
        raise ValueError(f'Cannot compute reduce with q unless q >= 2, but had q={q}.')


def is_pow_two(x: int) -> bool:
    """
    Check whether input integer is a (positive) power-of-two.
    :param x: Input integer
    :type x: int
    :return: Boolean indicating whether x is a positive power of two.
    :rtype: bool
    """
    if isinstance(x, int) and x > 0:
        return not (x & (x - 1))
    return False


class PolyCoefRep(object):
    # Coefficient representation of a polynomial modulo (q, X^n+1) for a prime q and a power-of-2 n
    q: int  # coefficient modulus
    half_q: int  # half the coefficient modulus
    n: int  # degree of polynomial modulus, must be power-of-2 (and we must track this many coefficients)
    vals: list[int]

    def __init__(self, q: int, n: int, vals: list[int]):
        q_is_prime: bool = is_prime(x=q)
        half_q: int = q//2
        if isinstance(q, int) and q > 1 and q_is_prime and isinstance(n, int) and n >= 1 and is_pow_two(x=n) and isinstance(vals, list) and len(vals) <= n and all(isinstance(x, int) for x in vals) and all(-half_q <= x <= half_q for x in vals):
            self.q = q
            self.half_q = half_q
            self.n = n
            self.vals = vals + [0 for _ in range(n - len(vals))]  # pad with zeros
        elif not isinstance(q, int):
            raise TypeError(f'Cannot instantiate a PolyCoefRep with a non-integer modulus, but had type(coef_modulus)={type(q)}.')
        elif q <= 1:
            raise ValueError(f'Cannot instantiate a PolyCoefRep with a modulus <=1, but had coef_modulus={q}.')
        elif not q_is_prime:
            raise ValueError(f'Cannot instantiate a PolyCoefRep with a non-prime modulus, but had coef_modulus={q}.')
        elif not isinstance(n, int):
            raise TypeError(f'Cannot instantiate a PolyCoefRep with a non-integer maximum degree n, but had n={n}.')
        elif n < 1:
            raise ValueError(f'Cannot instantiate a PolyCoefRep with maximum degree n < 1, but had n={n}.')
        elif not is_pow_two(x=n):
            raise ValueError(f'Cannot instantiate a PolyCoefRep with maximum degree n that is not a power of two, but had n={n}.')
        elif not isinstance(vals, list):
            raise TypeError(f'Cannot instantiate a PolyCoefRep unless input vals is a list, but had type(vals)={type(vals)}.')
        elif len(vals) > n:
            raise ValueError(f'Cannot instantiate a PolyCoefRep with more than n coefficients, but had len(vals)={len(vals)}.')
        elif not all(isinstance(x, int) for x in vals):
            raise TypeError(f'Cannot instantiate a PolyCoefRep unless input vals is a list of integers, but had [type(x) for x in vals]={[type(x) for x in vals]}.')
        raise ValueError(f'Cannot instantiate a PolyCoefRep unless input vals is a list of integers absolutely bounded by coef_modulus//2={half_q}, but had max(abs(x) for x in vals)={max(abs(x) for x in vals)}.')

    def __repr__(self):
        return f"PolyCoefRep{(self.q, self.n, self.vals)}"

    def __eq__(self, other):
        return isinstance(other, PolyCoefRep) and other.q == self.q and other.n == self.n == len(other.vals) == len(self.n) and all((x-y) % self.q == 0 for x, y in zip(self.vals, other.vals))

    def __add__(self, other):
        if isinstance(other, PolyCoefRep) and other.q == self.q and other.n == self.n == len(other.vals) == len(self.vals):
            if len(self.vals) >= len(other.vals):
                result: PolyCoefRep = deepcopy(self)
                for i, x in enumerate(other.vals):
                    result.vals[i] = reduce(x=result.vals[i] + x, q=self.q)
            else:
                result: PolyCoefRep = deepcopy(other)
                for i, x in enumerate(self.vals):
                    result.vals[i] = reduce(x=result.vals[i] + x, q = self.q)
            return result
        elif not isinstance(other, PolyCoefRep):
            raise NotImplementedError(f'Addition of PolyCoefRep objects with non-PolyCoefRep objects is not implemented.')
        elif other.q != self.q:
            raise NotImplementedError(f'Addition of PolyCoefRep with two different moduli is not implemented.')
        elif other.n != self.n:
            raise NotImplementedError(f'Addition of PolyCoefRep with two different maximum degrees is not implemented.')
        elif self.n != len(self.vals):
            raise ValueError(f'Cannot add PolyCoefRep unless both objects are well-formed with the same number of values as the maximum degree, but had (self.n, len(self.vals))={(self.n, len(self.vals))}.')
        raise ValueError(f'Cannot add PolyCoefRep unless both objects are well-formed with the same number of values as the maximum degree, but had (other.n, len(other.vals))={(other.n, len(other.vals))}.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [-x for x in negative_other.vals]
        return self.__add__(other=other)

    def __mul__(self, other):
        if isinstance(other, PolyCoefRep) and other.q == self.q:
            # warn("Warning. You are attempting to multiply two polynomial coefficient representations. This is slow: in the worst-case, we will FOIL, and in the best case we will first convert to the NTT form and then back again. In either case, it would be better to perform arithmetic with the NTT form.")
            result = deepcopy(self)
            result.vals = [0 for _ in range(len(self.vals) + len(other.vals) - 1)]
            for i, x in enumerate(self.vals):
                for j in range(i+1):
                    result.vals[i] += self.vals[j] * other.vals[i - j]
                result.vals[i] = reduce(x=result.vals[i], q=self.q)
            for i in range(self.n - 1):
                result.vals[i] = reduce(x=result.vals[i] - result.vals[self.n+i], q=self.q)
            return result
        elif not isinstance(other, PolyCoefRep):
            raise NotImplementedError(f'Multiplication of PolyCoefRep against a non-PolyCoefRep is not implemented.')
        raise NotImplementedError(f'Multiplication of PolyCoefRep with two different moduli is not implemented.')


class PolyHalfNTTRep(object):
    # ``Half''-NTT representation of a polynomial modulo (q, X^n+1) for a prime q and a power-of-2 n
    q: int  # coefficient modulus
    half_q: int  # half the coefficient modulus
    n: int  # degree of polynomial modulus (and we must track n/2 linear polynomials)
    vals: list[PolyCoefRep]

    def __init__(self, q: int, n: int, vals: list[int]):
        q_is_prime: bool = is_prime(x=q)
        half_q: int = q//2
        if isinstance(q, int) and \
                q > 1 and \
                q_is_prime and \
                isinstance(n, int) and \
                n >= 1 and \
                is_pow_two(x=n) and \
                isinstance(vals, list) and \
                len(vals) == n and \
                all(isinstance(x, PolyCoefRep) for x in vals) and \
                all(x.q == q for x in vals) and \
                all(x.n == 2 for x in vals) and \
                all(-half_q <= y <= half_q for x in vals for y in x.vals):
            self.q = q
            self.half_q = half_q
            self.n = n
            self.vals = vals
        elif not isinstance(q, int):
            raise TypeError(f'Cannot instantiate a PolyHalfNTTRep with a non-integer modulus, but had type(coef_modulus)={type(q)}.')
        elif q <= 1:
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep with a modulus <=1, but had coef_modulus={q}.')
        elif not q_is_prime:
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep with a non-prime modulus, but had coef_modulus={q}.')
        elif not isinstance(n, int):
            raise TypeError(f'Cannot instantiate a PolyHalfNTTRep with a non-integer maximum degree n, but had n={n}.')
        elif n < 1:
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep with maximum degree n < 1, but had n={n}.')
        elif not is_pow_two(x=n):
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep with maximum degree n that is not a power of two, but had n={n}.')
        elif not isinstance(vals, list):
            raise TypeError(f'Cannot instantiate a PolyHalfNTTRep unless input vals is a list, but had type(vals)={type(vals)}.')
        elif len(vals) != n:
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep without n values, but had len(vals)={len(vals)}.')
        elif not all(isinstance(x, PolyCoefRep) for x in vals):
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep unless input vals is a list of PolyCoefRep objects, but had [type(x) for x in vals)]={[type(x) for x in vals]}.')
        elif not all(x.q == q for x in vals):
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep unless input vals all have the same coefficient modulus, but had [min(x.q for x in vals), max(x.q for x in vals)]={[min(x.q for x in vals), max(x.q for x in vals)]}.')
        elif not all(x.n == 2 for x in vals):
            raise ValueError(f'Cannot instantiate a PolyHalfNTTRep unless input vals all have the same degree of polynomial modulus, but had [min(x.n for x in vals), max(x.n for x in vals)]={[min(x.n for x in vals), max(x.n for x in vals)]}.')
        raise TypeError(f'Cannot instantiate a PolyHalfNTTRep unless input vals is a list of PolyCoefRep objects whose coefficients are all absolutely bounded by half_q, but had max(abs(y) for x in vals for y in x.vals)={max(abs(y) for x in vals for y in x.vals)}.')

    def __repr__(self):
        return f"PolyHalfNTTRep{(self.q, self.n, self.vals)}"

    def __eq__(self, other):
        return isinstance(other, PolyHalfNTTRep) and other.q == self.q and other.n == self.n == len(other.vals) == len(self.n) and all(x == y for x, y in zip(self.vals, other.vals))

    def __add__(self, other):
        if isinstance(other, PolyHalfNTTRep) and \
                other.q == self.q and \
                other.n == self.n == len(other.vals) == len(self.vals) and \
                all(x.n == 2 for x in other.vals) and \
                all(x.n == 2 for x in self.vals) and \
                all(x.q == self.q for x in other.vals) and \
                all(x.q == self.q for x in self.vals):
            result: PolyHalfNTTRep = deepcopy(self)
            for i, x in enumerate(other.vals):
                result.vals[i] += x
            return result
        elif not isinstance(other, PolyCoefRep):
            raise NotImplementedError(f'Addition of PolyHalfNTTRep objects with non-PolyHalfNTTRep objects is not implemented.')
        elif other.q != self.q:
            raise NotImplementedError(f'Addition of PolyHalfNTTRep with two different coefficient moduli is not implemented.')
        elif other.n != self.n:
            raise NotImplementedError(f'Addition of PolyHalfNTTRep with two different degrees of polynomial moduli is not implemented.')
        elif self.n != len(self.vals):
            raise ValueError(f'Cannot add PolyHalfNTTRep unless both objects are well-formed with the same number of values as the maximum degree, but had (self.n, len(self.vals))={(self.n, len(self.vals))}.')
        raise ValueError(f'Cannot add PolyHalfNTTRep unless both objects are well-formed with the same number of values as the maximum degree, but had (other.n, len(other.vals))={(other.n, len(other.vals))}.')

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other=other)

    def __sub__(self, other):
        negative_other = deepcopy(other)
        negative_other.vals = [-x for x in negative_other.vals]
        return self.__add__(other=other)

    def __mul__(self, other):
        if isinstance(other, PolyHalfNTTRep) and \
                other.q == self.q and \
                other.n == self.n == len(other.vals) == len(self.vals) and \
                all(x.n == 2 for x in self.vals) and \
                all(x.n == 2 for x in other.vals) and \
                all(x.q == self.q for x in self.vals) and \
                all(x.q == self.q for x in other.vals):
            result = deepcopy(self)
            for i, other_val in enumerate(other.vals):
                result.vals[i] = result.vals[i] * other_val
            return result
        elif not isinstance(other, PolyHalfNTTRep):
            raise NotImplementedError(f'Multiplication of PolyHalfNTTRep against non-PolyHalfNTTRep objects is not implemented.')
        elif not other.q == self.q:
            raise NotImplementedError(f'Multiplication of PolyHalfNTTRep objects with different coefficient moduli is not implemented.')
        elif not other.n == self.n:
            raise NotImplementedError(f'Mulitplication of PolyHalfNTTRep objects with different degrees of polynomial moduli is not implemented.')
        elif not other.n == len(other.vals):
            raise ValueError(f'Cannot multiply PolyHalfNTTRep objects unless those objects are well-formed, but had (other.n, len(other.vals)) == {(other.n, len(other.vals))}.')
        elif not self.n == len(self.vals):
            raise ValueError(f'Cannot multiply PolyHalfNTTRep objects unless those objects are well-formed, but had (self.n, len(self.vals)) == {(self.n, len(self.vals))}.')
        elif not all(x.n == 2 for x in other.vals + self.vals):
            raise ValueError(f'Cannot multiply PolyHalfNTTRep objects unless both objects have values that are linear polynomials.')
        raise NotImplementedError(f'Multiplication of PolyHalfNTTRep where the values have different degrees of polynomial moduli is not implemented.')


def get_prim_rou_and_rou_inv(modulus: int, rou_deg: int) -> None | tuple[int, int]:
    if (modulus - 1) % rou_deg != 0:
        raise ValueError(f'Cannot compute primitive {rou_deg}-th root of unity for modulus={modulus} since modulus - 1 % rou_deg = {(modulus-1) % rou_deg} != 0.')
    x: int = 2
    while x < modulus:
        if all(x ** k % modulus != 1 for k in range(1, rou_deg)) and x ** rou_deg % modulus == 1:
            break
        x += 1
    return x, ((x ** (rou_deg - 1)) % modulus)


already_computed_zetas_and_inverses: dict[tuple[int, int]] = dict()


def make_zetas_and_invs(modulus: int, rou_deg: int) -> tuple[list[int], list[int]]:
    if (modulus, rou_deg) not in already_computed_zetas_and_inverses:
        lg_rou_deg: int = ceil(log2(rou_deg))
        powers: list[int] = [rou_deg // (2 ** (s + 1)) for s in range(lg_rou_deg)]
        zeta, zeta_inv = get_prim_rou_and_rou_inv(modulus=modulus, rou_deg=rou_deg)
        rou_powers: list[int] = [int(zeta ** i) % modulus for i in powers]
        rou_powers = [i if i <= modulus // 2 else i - modulus for i in rou_powers]
        rou_inverse_powers: list[int] = [int(zeta_inv ** i) % modulus for i in powers]
        rou_inverse_powers = [i if i <= modulus // 2 else i - modulus for i in rou_inverse_powers]
        already_computed_zetas_and_inverses[(modulus, rou_deg)] = (rou_powers, rou_inverse_powers)
    return already_computed_zetas_and_inverses[(modulus, rou_deg)]


def coef_to_half_ntt(x: PolyCoefRep) -> PolyHalfNTTRep:
    zetas, zeta_inverses = make_zetas_and_invs(modulus=x.q, rou_deg=x.n)


def half_ntt_to_coef(x: PolyHalfNTTRep) -> PolyCoefRep:
    zetas, zeta_inverses = make_zetas_and_invs(modulus=x.q, rou_deg=x.n)


def half_ntt(x: PolyCoefRep | PolyHalfNTTRep) -> PolyHalfNTTRep | PolyCoefRep:
    if isinstance(x, PolyCoefRep):
        return coef_to_half_ntt(x=x)
    elif isinstance(x, PolyHalfNTTRep):
        return half_ntt_to_coef(x=x)
    raise ValueError(f'Cannot compute HalfNTT unless input x is a PolyCoefRep object or a PolyHalfNTTRep object, but had type(x)={type(x)}.')
