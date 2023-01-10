from math import pi, e, ceil, log2, sqrt, floor
from crystals.kyber import is_pow_two, is_prime


def log2binom(d: int, k: int) -> int:
    j = k
    if k > d//2:
        j = d-k
    result: float = 0.0
    i: int = 0
    while i < j:
        result += log2(d-i)
        result -= log2(i+1)
        i += 1
    return result

CRYSTALS_TARGET: float = 2.42
ALLOWABLE_SECPARS: list[int] = [128, 192, 256, 512, 1024]
secpar_input: str = input(f'Please input a security parameter. secpar = ')
secpar: int = int(secpar_input)
if secpar not in ALLOWABLE_SECPARS:
    raise ValueError

log2delta: float = log2(((2*secpar+9)/0.265) * (pi * ((2*secpar+9)/0.265))**(1/((2*secpar+9)/0.265)) / (2*pi*e))/(2*(((2*secpar+9)/0.265) - 1))
print(f'For this secpar, we have log2delta={log2delta}.')

aggregation_capacity_input: str = input(f'Please input an aggregation capacity. K = ')
capacity: int = int(aggregation_capacity_input)
if capacity < 2:
    raise ValueError

d: int = 128

found_pairs_of_beta_omega_ch: dict[tuple[int, int], float] = {}

# find omega_ch, beta_ch that minimizes the product subject to the constraint that
# -log2binom(d=d, k=omega_ch) - omega_ch*log2(2*beta_ch+1) < -secpar
minval: float = -1.0
argmin: tuple[int, int] = tuple()
for next_beta in range(1, 1024):
    for next_omega in range(1, d + 1):
        if -log2binom(d=d, k=next_omega) - next_omega * log2(2 * next_beta + 1) < -secpar:
            found_pairs_of_beta_omega_ch[(next_beta, next_omega)] = next_omega * next_beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ch[
                (next_beta, next_omega)] <= minval:
                argmin = (next_beta, next_omega)
                minval = found_pairs_of_beta_omega_ch[(next_beta, next_omega)]

beta_ch: int = argmin[0]
omega_ch: int = argmin[1]
log2_epsilon_ch: float = -log2binom(d=d, k=omega_ch) - omega_ch * log2(2 * beta_ch + 1)

# find omega_ag, beta_ag that minimizes the product subject to the constraints that
# -log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1) < -1 and
# -log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1)
# -log2(1-epsilon_ch) -log2(1-2**(-log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1))) < -(secpar+1)
found_pairs_of_beta_omega_ag: dict[tuple[int, int], float] = {}
minval: float = -1.0
argmin: tuple[int, int] = tuple()
for next_beta in range(1, 1024):
    for next_omega in range(1, d + 1):
        log2_epsilon_tmp: float = -log2binom(d=d, k=next_omega) - next_omega * log2(2 * next_beta + 1)
        if -log2binom(d=d, k=next_omega) - next_omega * log2(2 * next_beta + 1) < -1 and \
                log2_epsilon_tmp - log2(1 - 2 ** log2_epsilon_ch) - log2(1 - 2 ** log2_epsilon_tmp) < -(
                secpar + 1):
            found_pairs_of_beta_omega_ag[(next_beta, next_omega)] = next_omega * next_beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ag[
                (next_beta, next_omega)] <= minval:
                argmin = (next_beta, next_omega)
                minval = found_pairs_of_beta_omega_ag[(next_beta, next_omega)]

beta_ag: int = argmin[0]
omega_ag: int = argmin[1]
log2_epsilon_ag: float = -log2binom(d=d, k=omega_ag) - omega_ag * log2(2 * beta_ag + 1)


found_pairs_of_beta_omega_sk: dict[tuple[int, int], float] = {}
minval: float = -1.0
argmin: tuple[int, int] = tuple()
for next_beta in range(1, 2**12):
    for next_omega in range(1, d+1):
        omega_v_prime: int = min(d, (1 + omega_ch) * next_omega)
        beta_v_prime: int = next_beta*(1+min(d,omega_ch,next_omega)*beta_ch)
        val: float = 2 * log2binom(d=d, k=next_omega) + 2*next_omega * log2(2 * next_beta + 1) - log2binom(d=d, k=omega_v_prime) - omega_v_prime * log2(2 * beta_v_prime + 1)
        if val > 0:
            found_pairs_of_beta_omega_sk[(next_beta, next_omega)] = next_omega*next_beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_sk[(next_beta, next_omega)] <= minval:
                argmin = (next_beta, next_omega)
                minval = found_pairs_of_beta_omega_sk[(next_beta, next_omega)]

if len(argmin) >= 2:
    beta_sk: int = argmin[0]
    omega_sk: int = argmin[1]
    omega_v_prime: int = min(d, (1 + omega_ch) * omega_sk)
    beta_v_prime: int = beta_sk*(1+min(d,omega_ch,omega_sk)*beta_sk)
    omega_v: int = min(d, capacity * omega_ag * omega_v_prime)
    beta_v: int = capacity * min(d, omega_ag, omega_v_prime) * beta_ag * beta_v_prime
    omega: int = min(d, 2 * omega_v + 2 * omega_ag * omega_v_prime)
    beta: int = 2 * beta_v + 2 * min(d, 2 * omega_ag, omega_v_prime) * beta_ag * beta_v_prime

    max_d: int = 1024
    while d <= max_d:
        # p: int = 2*d + 1
        p: int = 2**29
        p_minus_one_mod_two_d: int = 2*d - ((p-1) % (2*d))
        p += p_minus_one_mod_two_d
        assert ((p-1) % (2*d)) == 0
        while not is_prime(x=p):
            p += 2*d
        while log2(p) < 31:
            cond_9: bool = 8 * min(d, 2 * omega_ag, 4 * omega_ch * omega_sk) * min(d, 2 * omega_ch, 2 * omega_sk) * beta_ag * beta_ch * beta_sk < (p - 1) / 2
            if cond_9:
                ell: int = ceil((secpar + 2*d*log2(p))/(2*log2binom(d=d, k=omega_sk) + 2*omega_sk*log2(2*beta_sk+1) - log2binom(d=d, k=omega_v_prime) - omega_v_prime*log2(2*beta_v_prime+1)))
                print(f'secpar={secpar}, capacity={capacity}, d={d}, p={p}, omega_ag={omega_ag}, beta_ag={beta_ag}, omega_ch={omega_ch}, beta_ch={beta_ch}, omega_sk={omega_sk}, beta_sk={beta_sk}, ell={ell}')
                security_val: float = 0.5 * log2(d) + log2(beta) - log2(p) / ell
                security_val = security_val / (ell * d - 1)
                print(f'security_val={security_val}, target={log2delta}')
                while ell <= 1024 and security_val > log2delta:
                    ell += 1
                    security_val: float = 0.5 * log2(d) + log2(beta) - log2(p) / ell
                    security_val = security_val / (ell * d - 1)
                    print(f'security_val={security_val}, target={log2delta}')
                if security_val <= log2delta:
                    print(f'WINNER: secpar={secpar}, capacity={capacity}, d={d}, p={p}, omega_ag={omega_ag}, beta_ag={beta_ag}, omega_ch={omega_ch}, beta_ch={beta_ch}, omega_sk={omega_sk}, beta_sk={beta_sk}, ell={ell}')
                    with open("some_winners.txt", "a") as wf:
                        wf.write(str((secpar, capacity, d, p, omega_ag, beta_ag, omega_ch, beta_ch, omega_sk, beta_sk, ell)) + "\n")
            p += 2*d
            while not is_prime(x=p):
                p += 2*d
        d *= 2
