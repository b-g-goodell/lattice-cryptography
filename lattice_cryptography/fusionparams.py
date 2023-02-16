from math import pi, e, ceil, log2, sqrt, floor
from crystals.kyber import is_pow_two, is_prime
from copy import deepcopy

CRYSTALS_TARGET: float = 2.42
FALCON_TARGET: float = 0.666
seen_log2binom: dict[tuple[int, int], float] = dict()
this_p: int = 2**31 - 2**9 + 1

def log2binom(d: int, k: int) -> float:
    if (d, k) not in seen_log2binom:
        if k == d or k == 0:
            seen_log2binom[(d, k)] = 0.0
        else:
            j = k
            if k > d // 2:
                j = d - k
            result: float = 0.0
            i: int = 0
            while i < j:
                result += log2(d - i)
                result -= log2(i + 1)
                i += 1
            seen_log2binom[(d, k)] = result
    return seen_log2binom[(d, k)]


seen_omega_ch_and_beta_ch: dict[tuple[int, int, int], tuple[int, int] | None] = dict()


def find_omega_ch_and_beta_ch(secpar: int, d: int, max_beta: int) -> tuple[int, int] | None:
    if (secpar, d, max_beta) not in seen_omega_ch_and_beta_ch:
        minprod: int = -1
        argmin: tuple[int, int] = tuple()
        for next_beta in range(1, max_beta+1):
            for next_omega in range(1, d+1):
                if (minprod == -1 or next_beta * next_omega**2 <= minprod) and log2binom(d=d, k=next_omega) + next_omega*log2(2*next_beta) > secpar:
                    argmin = (next_omega, next_beta)
                    minprod = next_beta * next_omega**2
        if minprod == -1:
            seen_omega_ch_and_beta_ch[(secpar, d, max_beta)] = None
        else:
            seen_omega_ch_and_beta_ch[(secpar, d, max_beta)] = argmin
    return seen_omega_ch_and_beta_ch[(secpar, d, max_beta)]


seen_omega_ag_and_beta_ag: dict[tuple[int, int, int, int, int], tuple[int, int] | None] = dict()


def find_omega_ag_and_beta_ag(secpar: int, d: int, max_beta: int, omega_ch: int, beta_ch: int) -> tuple[int, int] | None:
    if (secpar, d, max_beta, omega_ch, beta_ch) not in seen_omega_ag_and_beta_ag:
        minprod: int = -1
        argmin: tuple[int, int] = tuple()
        for next_beta in range(1, max_beta+1):
            for next_omega in range(1, d+1):
                log2_of_one_minus_epsilon_ch: float = log2(1 - 2**(-log2binom(d=d, k=omega_ch) - omega_ch*log2(2*beta_ch)))
                log2_of_one_minus_epsilon_ag: float = log2(1 - 2**(-log2binom(d=d, k=next_omega) - next_omega*log2(2*next_beta)))
                log2_epsilon_ag: float = -log2binom(d=d, k=next_omega) - next_omega*log2(2*next_beta)
                if (minprod < 0 or next_beta * next_omega <= minprod) and (log2_epsilon_ag < -1 and log2_epsilon_ag - log2_of_one_minus_epsilon_ag - log2_of_one_minus_epsilon_ch < -(secpar+1)):
                    argmin = (next_omega, next_beta)
                    minprod = next_beta * next_omega
        if minprod == -1:
            seen_omega_ag_and_beta_ag[(secpar, d, max_beta, omega_ch, beta_ch)] = None
        else:
            seen_omega_ag_and_beta_ag[(secpar, d, max_beta, omega_ch, beta_ch)] = argmin
    return seen_omega_ag_and_beta_ag[(secpar, d, max_beta, omega_ch, beta_ch)]


seen_omega_sk_and_beta_sk: dict[tuple[int, int, int, int, int, int, int], tuple[int, int] | None] = dict()


def find_omega_sk_and_beta_sk(secpar: int, d: int, omega_ch: int, beta_ch: int, omega_ag: int, beta_ag: int, capacity: int) -> tuple[int, int] | None:
    if (secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity) not in seen_omega_sk_and_beta_sk:
        minprod: int = -1
        argmin: tuple[int, int] = tuple()
        omega_range: list[int] = list(range(d, 0, -1))
        for next_omega in omega_range:
            new_max_beta: int = floor(((2**31 - 1) / 2) / (8 * min(d, 2 * omega_ag, 4 * omega_ch * next_omega) * min(d, 2 * omega_ch,2 * next_omega) * beta_ag * beta_ch))
            another_max_beta: int = floor(((2**31 - 1) / 2) / (2 * (capacity * min(d, omega_ag, min(d, (1 + omega_ch) * next_omega)) + min(d, 2 * omega_ag, min(d, (1 + omega_ch) * next_omega))) * beta_ag * (1 + min(d, omega_ch, next_omega) * beta_ch)))
            max_beta: int = min(new_max_beta, another_max_beta)
            if 1 < max_beta:
                for next_beta in range(1, max_beta):
                    omega_v_prime: int = next_omega*(1+omega_ch)
                    beta_v_prime: int = next_beta*(1+min(d, omega_ch, next_omega)*beta_ch)
                    condition: bool = 2*log2binom(d=d, k=next_omega) + 2*next_omega * log2(2*next_beta) > max([log2binom(d=d, k=i) + i*log2(2*beta_v_prime) for i in range(omega_v_prime+1)])
                    if condition and (minprod < 0 or minprod > next_beta * next_omega**2):
                        minprod = next_beta * next_omega**2
                        argmin = (next_omega, next_beta)
        if minprod == -1:
            seen_omega_sk_and_beta_sk[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity)] = None
        else:
            seen_omega_sk_and_beta_sk[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity)] = argmin
    return seen_omega_sk_and_beta_sk[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity)]


seen_primes: dict[tuple[int, int, int, int, int, int, int, int, int], int | None] = dict()


def find_prime(secpar: int, d: int, omega_ch: int, beta_ch: int, omega_ag: int, beta_ag: int, capacity: int, omega_sk: int, beta_sk: int) -> int | None:
    if (secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk) not in seen_primes:
        one_lower_bd: int = 16*min(d, 2*omega_ag, 4*omega_ch*omega_sk)*min(d, 2*omega_ch,2*omega_sk)*beta_ag*beta_ch*beta_sk + 1
        omega_v_prime: int = omega_sk * (1 + omega_ch)
        beta_v_prime: int = beta_sk * (1 + min(d, omega_ch, omega_sk) * beta_ch)
        omega_v: int = min(d, capacity*omega_ag*omega_v_prime)
        beta_v: int = capacity*min(d,omega_ag,omega_v_prime)*beta_ag*beta_v_prime
        omega: int = min(d, 2*omega_v + 2*omega_ag*omega_v_prime)
        beta: int = 2*beta_v + 2*min(2, 2*omega_ag, omega_v_prime)*beta_ag*beta_v_prime
        another_lower_bd: int = 2*beta+1
        starting_p: int = max(one_lower_bd, another_lower_bd) + 1
        starting_p += 2*d - ((starting_p - 1) % (2*d))
        p: int = starting_p
        while not is_prime(x=p) and log2(p) <= 31:
            p += 2*d
        seen_primes[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk)] = p
        if not is_prime(x=p) or log2(p) > 31:
            seen_primes[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk)] = None
    return seen_primes[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk)]



seen_ell: dict[tuple[int, int, int, int, int, int, int, int, int, int], int | None] = dict()


def find_ell(secpar: int, d: int, omega_ch: int, beta_ch: int, omega_ag: int, beta_ag: int, capacity: int, omega_sk: int, beta_sk: int, p: int) -> int|None:
    if (secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk, p) not in seen_ell:
        omega_v_prime: int = omega_sk * (1 + omega_ch)
        beta_v_prime: int = beta_sk * (1 + min(d, omega_ch, omega_sk) * beta_ch)
        omega_v: int = min(d, capacity*omega_ag*omega_v_prime)
        beta_v: int = capacity*min(d,omega_ag,omega_v_prime)*beta_ag*beta_v_prime
        omega: int = min(d, 2*omega_v + 2*omega_ag*omega_v_prime)
        beta: int = 2*beta_v + 2*min(2, 2*omega_ag, omega_v_prime)*beta_ag*beta_v_prime
        existence_val: float = 2 * log2binom(d=d, k=omega_sk) + 2 * omega_sk * log2(2 * beta_sk) - d * log2(2 * beta_v_prime + 1)
        existence_condition: bool = existence_val > 0
        if not existence_condition:
            seen_ell[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk, p)] = None
        else:
            log2delta: float = (log2(2*secpar+9)-1-log2(pi)-log2(e)-log2(0.265)+(0.265/(2*secpar+9))*(log2(pi)+log2(2*secpar+9)-log2(0.265)))/(2*((2*secpar+9)/0.265 - 1))
            ell: int = 1
            while existence_val * ell < secpar + 2*d*log2(p) or (log2(d)/2 + log2(beta) - log2(p)/ell)/(ell*d-1) >= log2delta:
                ell += 1
            assert verify(secpar=secpar, p=p, d=d, capacity=capacity, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, omega_sk=omega_sk, beta_sk=beta_sk, ell=ell)
            seen_ell[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk, p)] = ell
    return seen_ell[(secpar, d, omega_ch, beta_ch, omega_ag, beta_ag, capacity, omega_sk, beta_sk, p)]


def verify(secpar: int, p: int, d: int, capacity: int, omega_ch: int, beta_ch: int, omega_ag: int, beta_ag: int, omega_sk: int, beta_sk: int, ell: int) -> bool:
    log2delta: float = (log2(2*secpar+9)-1-log2(pi)-log2(e)-log2(0.265)+(0.265/(2*secpar+9))*(log2(pi)+log2(2*secpar+9)-log2(0.265)))/(2*((2*secpar+9)/0.265 - 1))
    omega_v_prime: int = min(d, (1 + omega_ch) * omega_sk)
    beta_v_prime: int = beta_sk * (1 + min(d, omega_ch, omega_sk) * beta_ch)
    omega_v: int = min(d, capacity * omega_ag * omega_v_prime)
    beta_v: int = capacity * min(d, omega_ag, omega_v_prime) * beta_ag * beta_v_prime
    omega: int = min(d, 2*omega_v + 2*omega_ag*omega_v_prime)
    beta: int = 2*beta_v+2*min(d,2*omega_ag,omega_v_prime)*beta_ag*beta_v_prime
    log2_epsilon_ch: float = -log2binom(d=d,k=omega_ch) - omega_ch*log2(2*beta_ch+1)
    log2_epsilon_ag: float = -log2binom(d=d,k=omega_ag) - omega_ag*log2(2*beta_ag+1)
    one_minus_epsilon_ch: float = 1 - 2**log2_epsilon_ch
    one_minus_epsilon_ag: float = 1 - 2**log2_epsilon_ag
    log2_adjusted_epsilon_ag: float = log2_epsilon_ag - log2(one_minus_epsilon_ch) - log2(one_minus_epsilon_ag)  # Improve this with taylor expansion
    ntt_friendly: bool = ((p-1) % (2*d) == 0)
    beta_small_enough: bool = beta < (p-1)/2
    combo_small_enough: bool = 8*min(d,2*omega_ag,4*omega_ch*omega_sk)*min(d,2*omega_ch,2*omega_sk)*beta_ag*beta_ch*beta_sk < (p-1)/2
    security_holds: bool = (log2(d)/2 + log2(beta) - log2(p)/ell)/(ell*d-1) < log2delta
    lemma_condition: bool = secpar + 2*d*log2(p) <= ell*(2*log2binom(d=d, k=omega_sk) + 2*omega_sk*log2(2*beta_sk) - d*log2(2*beta_v_prime+1))
    eps_ch_cond: bool = log2_epsilon_ch < -secpar
    eps_ag_cond: bool = log2_epsilon_ag < -1 and log2_adjusted_epsilon_ag < -(secpar+1)
    return ntt_friendly and beta_small_enough and combo_small_enough and security_holds and lemma_condition and eps_ch_cond and eps_ag_cond


seen_parameters: dict[tuple[int, int, int, int], tuple | None] = dict()
MAX_BETA_FOR_CH: int = 2**9
MAX_BETA_FOR_AG: int = 2**9


def find_params(secpar: int, d: int, capacity: int):
    omega_ch_and_beta_ch = find_omega_ch_and_beta_ch(secpar=secpar, d=d, max_beta=2**9)
    if omega_ch_and_beta_ch is not None:
        omega_ch, beta_ch = omega_ch_and_beta_ch
        omega_ag_and_beta_ag = find_omega_ag_and_beta_ag(secpar=secpar, d=d, max_beta=2**9, omega_ch=omega_ch, beta_ch=beta_ch)
        if omega_ag_and_beta_ag is not None:
            omega_ag, beta_ag = omega_ag_and_beta_ag
            omega_sk_and_beta_sk = find_omega_sk_and_beta_sk(secpar=secpar, d=d, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, capacity=capacity)
            if omega_sk_and_beta_sk is not None:
                omega_sk, beta_sk = omega_sk_and_beta_sk
                p = find_prime(secpar=secpar, d=d, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, capacity=capacity, omega_sk=omega_sk, beta_sk=beta_sk)
                if p is not None:
                    ell = find_ell(secpar=secpar, d=d, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, capacity=capacity, omega_sk=omega_sk, beta_sk=beta_sk, p=p)
                    if ell is not None:
                        print(f'secpar={secpar},d={d},capacity={capacity},omega_ch={omega_ch},beta_ch={beta_ch},omega_ag={omega_ag},beta_ag={beta_ag},omega_sk={omega_sk},beta_sk={beta_sk},p={p},ell={ell}')


for i in range(7,10):
    for j in range(6,10):
        for k in range(1,10):
            find_params(secpar=2**i, d=2**j, capacity=2**k)
#
# seen_forced_parameters: dict = dict()
#
#
# def find_parameters_with_forced_choices(secpar: int, p: int, d: int, capacity: int, omega_ch: int | None = None, beta_ch: int | None = None, omega_ag: int | None = None, beta_ag: int | None = None, omega_sk: int | None = None) -> list | None:
#     if (secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag, omega_sk) not in seen_forced_parameters:
#         seen_forced_parameters[(secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag)] = None
#         omega_sk: int
#         beta_sk: int
#         sk_out: tuple[int, int] | None = find_omega_sk_and_beta_sk(secpar=secpar, d=d, omega_ch=omega_ch,
#                                                                    beta_ch=beta_ch, omega_ag=omega_ag,
#                                                                    beta_ag=beta_ag, capacity=capacity, p=p, omega_sk=omega_sk)
#         if sk_out is not None:
#             omega_sk_observed, beta_sk_observed = sk_out
#             if omega_sk is not None:
#                 assert omega_sk_observed == omega_sk
#             ell: int | None = find_ell(secpar=secpar, d=d, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag,
#                                        beta_ag=beta_ag, capacity=capacity, omega_sk=omega_sk_observed, beta_sk=beta_sk_observed, p=p)
#             if ell is not None and verify(secpar=secpar, p=p, d=d, capacity=capacity, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, omega_sk=omega_sk_observed, beta_sk=beta_sk_observed, ell=ell):
#                 omega_v_prime: int = min(d, (1 + omega_ch) * omega_sk)
#                 beta_v_prime: int = beta_sk_observed * (1 + min(d, omega_ch, omega_sk_observed) * beta_ch)
#                 omega_v: int = min(d, capacity * omega_ag * omega_v_prime)
#                 beta_v: int = capacity * min(d, omega_ag, omega_v_prime) * beta_ag * beta_v_prime
#                 vk_weight: float = 2 * d * log2(p)
#                 # agg_sig_weight: float = ell*(d + omega_v*ceil(log2(2*beta_v+1)))  # if we use a d-bit string
#                 # agg_sig_weight: float = ell*(omega_v*ceil(log2(d)) + omega_v*ceil(log2(2*beta_v+1)))  # if we use a list of omega_v d-bit indices
#                 # agg_sig_weight: float = ell*d*ceil(log2(2*beta_v+1))  # if we use a list of d coefficients and count the zeros manually
#                 agg_sig_weight: float = ell * min(d + omega_v * ceil(log2(2 * beta_v + 1)),
#                                                   omega_v * ceil(log2(d)) + omega_v * ceil(log2(2 * beta_v + 1)),
#                                                   d * ceil(log2(2 * beta_v + 1)))
#                 avg_weight: float = vk_weight + agg_sig_weight / capacity
#                 avg_weight_kb: float = avg_weight / 8000
#
#                 seen_forced_parameters[(secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag)] = (secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag, omega_sk, beta_sk_observed, ell, avg_weight_kb)
#     return seen_forced_parameters[(secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag)]
#
#
# def find_parameters(secpar: int, p: int, d: int, capacity: int) -> list | None:
#     if (secpar, p, d, capacity) not in seen_parameters:
#         omega_ch: int
#         beta_ch: int
#         ch_out: tuple[int, int] | None = find_omega_ch_and_beta_ch(secpar=secpar, d=d, max_beta=MAX_BETA_FOR_CH)
#         if ch_out is None:
#             seen_parameters[(secpar, p, d, capacity)] = None
#         else:
#             omega_ch, beta_ch = ch_out
#             omega_ag: int
#             beta_ag: int
#             ag_out: tuple[int, int] | None = find_omega_ag_and_beta_ag(secpar=secpar, d=d, max_beta=MAX_BETA_FOR_AG, omega_ch=omega_ch, beta_ch=beta_ch)
#             if ag_out is None:
#                 seen_parameters[(secpar, p, d, capacity)] = None
#             else:
#                 omega_ag, beta_ag = ag_out
#                 omega_sk: int
#                 beta_sk: int
#                 sk_out: tuple[int, int] | None = find_omega_sk_and_beta_sk(secpar=secpar, d=d, omega_ch=omega_ch,
#                                                                            beta_ch=beta_ch, omega_ag=omega_ag,
#                                                                            beta_ag=beta_ag, capacity=capacity, p=p)
#                 if sk_out is None:
#                     seen_parameters[(secpar,p,d,capacity)] = None
#                 else:
#                     omega_sk, beta_sk = sk_out
#                     ell: int | None = find_ell(secpar=secpar, d=d, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, capacity=capacity, omega_sk=omega_sk, beta_sk=beta_sk, p=p)
#                     if ell is None:
#                         seen_parameters[(secpar,p,d,capacity)] = None
#                     elif verify(secpar=secpar, p=p, d=d, capacity=capacity, omega_ch=omega_ch, beta_ch=beta_ch, omega_ag=omega_ag, beta_ag=beta_ag, omega_sk=omega_sk, beta_sk=beta_sk, ell=ell):
#                         omega_v_prime: int = min(d, (1 + omega_ch) * omega_sk)
#                         beta_v_prime: int = beta_sk * (1 + min(d, omega_ch, omega_sk) * beta_ch)
#                         omega_v: int = min(d, capacity*omega_ag*omega_v_prime)
#                         beta_v: int = capacity*min(d, omega_ag, omega_v_prime)*beta_ag*beta_v_prime
#                         vk_weight: float = 2*d*log2(p)
#                         # agg_sig_weight: float = ell*(d + omega_v*ceil(log2(2*beta_v+1)))  # if we use a d-bit string
#                         # agg_sig_weight: float = ell*(omega_v*ceil(log2(d)) + omega_v*ceil(log2(2*beta_v+1)))  # if we use a list of omega_v d-bit indices
#                         # agg_sig_weight: float = ell*d*ceil(log2(2*beta_v+1))  # if we use a list of d coefficients and count the zeros manually
#                         agg_sig_weight: float = ell*min(d + omega_v*ceil(log2(2*beta_v+1)), omega_v*ceil(log2(d)) + omega_v*ceil(log2(2*beta_v+1)), d*ceil(log2(2*beta_v+1)))
#                         avg_weight: float = vk_weight + agg_sig_weight/capacity
#                         avg_weight_kb: float = avg_weight/8000
#
#                         seen_parameters[(secpar, p, d, capacity)] = (secpar, p, d, capacity, omega_ch, beta_ch, omega_ag, beta_ag, omega_sk, beta_sk, ell, avg_weight_kb)
#     return seen_parameters[(secpar, p, d, capacity)]
#
#
# ALLOWABLE_SECPARS: list[int] = [192]
# ALLOWABLE_DEGREE: list[int] = [32, 64, 128]
# ALLOWABLE_CAPACITIES: list[int] = [2]  # list(range(2, 64))
#
#
# def find_all_parameters() -> list:
#     p: int = this_p
#     result: dict = dict()
#     winners: dict = dict()
#     for secpar in ALLOWABLE_SECPARS:
#         for d in ALLOWABLE_DEGREE:
#             min_avg_size: float = -1.0
#             for capacity in ALLOWABLE_CAPACITIES:
#                 print((secpar, d, capacity))
#                 winner: list[int] = find_parameters(secpar=secpar, d=d, p=p, capacity=capacity)
#                 result[(secpar, d, p, capacity)] = winner
#                 if winner is not None and (min_avg_size < 0.0 or min_avg_size > winner[-1]) and \
#                         verify(secpar=winner[0], d=winner[2], p=winner[1], capacity=winner[3], omega_ch=winner[4],
#                                beta_ch=winner[5], omega_ag=winner[6], beta_ag=winner[7], omega_sk=winner[8],
#                                beta_sk=winner[9], ell=winner[10]):
#                     winners[(secpar, d, p)] = winner
#                     min_avg_size = winner[-1]
#     for next_secpar_d_p in winners:
#         winner = winners[next_secpar_d_p]
#         if winner is not None and len(winner) > 0 and verify(secpar=winner[0], d=winner[2], p=winner[1], capacity=winner[3], omega_ch=winner[4], beta_ch=winner[5], omega_ag=winner[6], beta_ag=winner[7], omega_sk=winner[8], beta_sk=winner[9], ell=winner[10]):
#             file_output: str = f'We use secpar={winner[0]}, p={winner[1]}, and d={winner[2]}. \n'
#             file_output += f'We set K={winner[3]}. The pair (omega_ch, beta_ch) = {(winner[4], winner[5])} minimizes omega_ch*beta_ch subject to the constraint (12). \n'
#             file_output += f'The pair (omega_ag, beta_ag)={(winner[6], winner[7])} minimizes omega_ag*beta_ag subject to the constraints (13) and (14). \n'
#             file_output += f'The pair (omega_sk, beta_sk)={(winner[8], winner[9])} maximizes binom(d,omega_sk)^2(2beta_sk+1)^2omega_sk/binom(d, omega_v_prime)(2beta_v_prime+1)^omega_v_prime > 1 subject to (8) and (9). \n'
#             file_output += f'Lastly, ell = {winner[10]} is the minimal value satisfying (10) and (11).\n'
#             file_output += f'secpar =& {winner[0]}\n'
#             file_output += f'p =& {winner[1]}\n'
#             file_output += f'd =& {winner[2]}\n'
#             file_output += f'K =& {winner[3]}\n'
#             file_output += f'omega_ch =& {winner[4]}\n'
#             file_output += f'beta_ch =& {winner[5]}\n'
#             file_output += f'omega_ag =& {winner[6]}\n'
#             file_output += f'beta_ag =& {winner[7]}\n'
#             file_output += f'omega_sk =& {winner[8]}\n'
#             file_output += f'beta_sk =& {winner[9]}\n'
#             file_output += f'ell =& {winner[10]}\n'
#             omega_v_prime: int = min(winner[2], winner[8]*(1+winner[4]))
#             beta_v_prime: int = winner[9]*(1+min(winner[2], winner[8], winner[4])*winner[5])
#             file_output += f'omega_v_prime =& {omega_v_prime}\n'
#             file_output += f'beta_v_prime =& {beta_v_prime}\n'
#             omega_v: int = min(winner[2], winner[3]*winner[6]*omega_v_prime)
#             beta_v: int = winner[3]*min(winner[2], winner[6], omega_v_prime)*winner[7]*beta_v_prime
#             file_output += f'omega_v =& {omega_v}\n'
#             file_output += f'beta_v =& {beta_v}\n'
#             omega: int = min(d, 2*omega_v+2*winner[6]*omega_v_prime)
#             beta: int = 2*beta_v + 2*min(d, 2*winner[6], omega_v_prime)*winner[7]*beta_v_prime
#             file_output += f'omega =& {omega}\n'
#             file_output += f'beta =& {beta}\n'
#             file_output += f'We have verification keys of size 2dlog(p)={2*winner[2]*ceil(log2(winner[1]))} bits, which are {2*winner[2]*ceil(log2(winner[1]))/8000} kilobytes. \n'
#             file_output += f'We have signatures that are ell*d*ceil(log(2*beta_v_prime+1))={winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))} bits, which are {winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))/8000} kilobytes. \n'
#             file_output += f'However, K={winner[3]} of these aggregate together. An aggregate signature is ell*d*ceil(log(2*beta_v+1))={winner[10]*winner[2]*ceil(log2(2*beta_v+1))} bits, which is {winner[10]*winner[2]*ceil(log2(2*beta_v+1))/8000} kilobytes. \n'
#             file_output += f'Together, K verification keys and one aggregate signature takes {winner[3]*(2*winner[2]*log2(winner[1])/8000) + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/8000} kilobytes. \n'
#             file_output += f'This cost in space averages out to {2*winner[2]*log2(winner[1])/8000 + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/(winner[3]*8000)} kilobytes per signer.\n'
#             if winner[-1] <= CRYSTALS_TARGET:
#                 file_output += f'Together, Fusion keys and aggregate signatures have better space efficiency than the smallest CRYSTALS-Dilithium signatures, which are {CRYSTALS_TARGET} kilobytes. \n'
#                 ct: int = 2
#                 while 2*winner[2]*log2(winner[1])/8000 + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/(ct*8000) > CRYSTALS_TARGET and ct < winner[3]:
#                     ct += 1
#                 if ct < winner[3]:
#                     file_output += f'In fact, if we do not aggregate at full capacity, then aggregating only {ct} signatures (i.e. occupying {100*ct/winner[3]}\% of signature capacity) is sufficient to beat CRYSTALS-Dilithium efficiency.\n '
#             if winner[-1] <= FALCON_TARGET:
#                 file_output += f'Together, Fusion keys and aggregate signatures have better space efficiency than the smallest Falcon signatures, which are {FALCON_TARGET} kilobytes.\n '
#                 ct: int = 2
#                 while 2 * winner[2] * log2(winner[1]) / 8000 + winner[10] * winner[2] * ceil(log2(2 * beta_v + 1)) / (ct * 8000) > FALCON_TARGET and ct < winner[3]:
#                     ct += 1
#                 if ct < winner[3]:
#                     file_output += f'In fact, if we do not aggregate at full capacity, then aggregating only {ct} signatures (i.e. occupying {100 * ct / winner[3]}\% of signature capacity) is sufficient to beat FALCON efficiency.\n '
#             file_output += f'The cost of using this scheme occurs in communication overhead: using these parameters, an aggregator must download K={winner[3]} public verification keys, signatures, and messages. \n'
#             file_output += f'Disregarding the cost of messages, this means downloading {winner[3]*2*winner[2]*log2(winner[1])/8000 + winner[3]*winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))/8000} kilobytes.\n'
#
#             print(f'secpar={winner[0]}, p={winner[1]}, d={winner[2]}, capacity={winner[3]}, omega_ch={winner[4]}, beta_ch={winner[5]}, omega_ag={winner[6]}, beta_ag={winner[7]}, omega_sk={winner[8]}, beta_sk={winner[9]}, ell={winner[10]}, avg_size={winner[11]}')
#             next_secpar_d_p=(winner[0], winner[2], winner[1])
#             with open(f"results{next_secpar_d_p}.txt", "w") as wf:
#                 wf.write("secpar,p,d,capacity,omega_ch,beta_ch,omega_ag,beta_ag,omega_sk,beta_sk,ell,avg_weight_kb\n" + str(file_output) + "\n\n\n")
#
#
# # find_all_parameters()
#
# some_winners: list = [
#     [128, 2**31-2**17+1, 64, 256, 30, 2, 30, 2, 64],
#     [128, 2**31-2**17+1, 256, 67268, 20, 1, 20, 1, 256],
#     [192, 2**31-2**17+1, 256, 10744, 34, 1, 34, 1, 256],
#     [256, 2**31-2**17+1, 256, 4673, 50, 1, 50, 1, 256],
#     [384, 2**31-2**17+1, 256, 847, 93, 1, 93, 1, 256],
#     [512, 2**31-2**17+1, 256, 59, 114, 2, 114, 2, 256],
# ]
#
# assert verify(secpar=128, p=2**31-2**17+1, d=64, capacity=1766, omega_ch=30, beta_ch=2, omega_ag=30, beta_ag=2, omega_sk=64, beta_sk=81, ell=174)
#
# updated_winners = []
#
# other_updated_winners = []
# #
# # for next_winner in some_winners:
# #     # minimize p
# #     candidate = find_parameters_with_forced_choices(secpar=next_winner[0], p=next_winner[1], d=next_winner[2], capacity=next_winner[3], omega_ch=next_winner[4], beta_ch=next_winner[5], omega_ag=next_winner[6], beta_ag=next_winner[7], omega_sk=next_winner[8])
# #     landscape = [candidate]
# #     assert verify(secpar=candidate[0], p=candidate[1], d=candidate[2], capacity=candidate[3], omega_ch=candidate[4], beta_ch=candidate[5], omega_ag=candidate[6], beta_ag=candidate[7], omega_sk=candidate[8], beta_sk=candidate[9], ell=candidate[10])
# #     if candidate is not None:
# #         other_updated_winners += [candidate]
# #         this_cost = candidate[-1]
# #         next_prime = 2**floor(log2(candidate[1])) - 1
# #         next_prime -= ((next_prime - 1) % (2*candidate[2]))
# #         while not is_prime(x=next_prime):
# #             next_prime -= 2 * candidate[2]
# #         next_candidate = find_parameters_with_forced_choices(secpar=candidate[0], p=next_prime, d=candidate[2], capacity=candidate[3], omega_ch=candidate[4],beta_ch=candidate[5], omega_ag=candidate[6],beta_ag=candidate[7], omega_sk=candidate[8])
# #         if next_candidate is not None:
# #             next_cost = next_candidate[-1]
# #             landscape += [next_candidate]
# #             while candidate is not None and \
# #                     next_candidate is not None and \
# #                     verify(secpar=candidate[0], p=candidate[1], d=candidate[2], capacity=candidate[3], omega_ch=candidate[4], beta_ch=candidate[5], omega_ag=candidate[6], beta_ag=candidate[7], omega_sk=candidate[8], beta_sk=candidate[9], ell=candidate[10]) and \
# #                     verify(secpar=next_candidate[0], p=next_candidate[1], d=next_candidate[2], capacity=next_candidate[3], omega_ch=next_candidate[4], beta_ch=next_candidate[5], omega_ag=next_candidate[6], beta_ag=next_candidate[7], omega_sk=next_candidate[8], beta_sk=next_candidate[9], ell=next_candidate[10]):
# #                 candidate = next_candidate
# #                 this_cost = candidate[-1]
# #                 next_prime = 2 ** floor(log2(candidate[1])) - 1
# #                 next_prime -= ((next_prime - 1) % (2 * candidate[2]))
# #                 while not is_prime(x=next_prime):
# #                     next_prime -= 2 * candidate[2]
# #                 next_candidate = find_parameters_with_forced_choices(secpar=candidate[0], p=next_prime, d=candidate[2],
# #                                                                      capacity=candidate[3], omega_ch=candidate[4],
# #                                                                      beta_ch=candidate[5], omega_ag=candidate[6],
# #                                                                      beta_ag=candidate[7], omega_sk=candidate[8])
# #                 if next_candidate is not None:
# #                     next_cost = next_candidate[-1]
# #                     landscape += [next_candidate]
# #             landscape = sorted(landscape, key=lambda x: x[-1], reverse=False)
# #             other_updated_winners[-1] = landscape[0]
# #
#
# for next_winner in some_winners:
#     # maximize capacity
#     candidate = find_parameters_with_forced_choices(secpar=next_winner[0], p=next_winner[1], d=next_winner[2], capacity=next_winner[3], omega_ch=next_winner[4], beta_ch=next_winner[5], omega_ag=next_winner[6], beta_ag=next_winner[7], omega_sk=next_winner[8])
#     assert verify(secpar=candidate[0], p=candidate[1], d=candidate[2], capacity=candidate[3], omega_ch=candidate[4], beta_ch=candidate[5], omega_ag=candidate[6], beta_ag=candidate[7], omega_sk=candidate[8], beta_sk=candidate[9], ell=candidate[10])
#     if candidate is not None:
#         updated_winners += [candidate]
#         this_cost = candidate[-1]
#         next_capacity = candidate[3] + 1
#         next_candidate = find_parameters_with_forced_choices(secpar=candidate[0], p=candidate[1], d=candidate[2], capacity=next_capacity, omega_ch=candidate[4],beta_ch=candidate[5], omega_ag=candidate[6],beta_ag=candidate[7], omega_sk=candidate[8])
#         if next_candidate is not None:
#             next_cost = next_candidate[-1]
#             while candidate is not None and \
#                     next_candidate is not None and \
#                     next_cost < this_cost and \
#                     verify(secpar=candidate[0], p=candidate[1], d=candidate[2], capacity=candidate[3], omega_ch=candidate[4], beta_ch=candidate[5], omega_ag=candidate[6], beta_ag=candidate[7], omega_sk=candidate[8], beta_sk=candidate[9], ell=candidate[10]) and \
#                     verify(secpar=next_candidate[0], p=next_candidate[1], d=next_candidate[2], capacity=next_candidate[3], omega_ch=next_candidate[4], beta_ch=next_candidate[5], omega_ag=next_candidate[6], beta_ag=next_candidate[7], omega_sk=next_candidate[8], beta_sk=next_candidate[9], ell=next_candidate[10]):
#                 candidate = next_candidate
#                 this_cost = candidate[-1]
#                 next_capacity += 1
#                 next_candidate = find_parameters_with_forced_choices(secpar=candidate[0], p=candidate[1], d=candidate[2],
#                                                                      capacity=next_capacity, omega_ch=candidate[4],
#                                                                      beta_ch=candidate[5], omega_ag=candidate[6],
#                                                                      beta_ag=candidate[7], omega_sk=candidate[8])
#                 if next_candidate is not None:
#                     next_cost = next_candidate[-1]
#             updated_winners[-1] = candidate
#
#
# for next_winner in updated_winners:
#     winner = find_parameters_with_forced_choices(secpar=next_winner[0], p=next_winner[1], d=next_winner[2], capacity=next_winner[3], omega_ch=next_winner[4], beta_ch=next_winner[5], omega_ag=next_winner[6], beta_ag=next_winner[7], omega_sk=next_winner[8])
#     file_output: str = f'We use secpar={winner[0]}, p={winner[1]}, and d={winner[2]}. \n'
#     file_output += f'We set K={winner[3]}. The pair (omega_ch, beta_ch) = {(winner[4], winner[5])} minimizes omega_ch*beta_ch subject to the constraint (12). \n'
#     file_output += f'The pair (omega_ag, beta_ag)={(winner[6], winner[7])} minimizes omega_ag*beta_ag subject to the constraints (13) and (14). \n'
#     file_output += f'The pair (omega_sk, beta_sk)={(winner[8], winner[9])} maximizes binom(d,omega_sk)^2(2beta_sk+1)^2omega_sk/binom(d, omega_v_prime)(2beta_v_prime+1)^omega_v_prime > 1 subject to (8) and (9). \n'
#     file_output += f'Lastly, ell = {winner[10]} is the minimal value satisfying (10) and (11).\n'
#     file_output += f'secpar =& {winner[0]}\n'
#     file_output += f'p =& {winner[1]}\n'
#     file_output += f'd =& {winner[2]}\n'
#     file_output += f'K =& {winner[3]}\n'
#     file_output += f'omega_ch =& {winner[4]}\n'
#     file_output += f'beta_ch =& {winner[5]}\n'
#     file_output += f'omega_ag =& {winner[6]}\n'
#     file_output += f'beta_ag =& {winner[7]}\n'
#     file_output += f'omega_sk =& {winner[8]}\n'
#     file_output += f'beta_sk =& {winner[9]}\n'
#     file_output += f'ell =& {winner[10]}\n'
#     omega_v_prime: int = min(winner[2], winner[8]*(1+winner[4]))
#     beta_v_prime: int = winner[9]*(1+min(winner[2], winner[8], winner[4])*winner[5])
#     file_output += f'omega_v_prime =& {omega_v_prime}\n'
#     file_output += f'beta_v_prime =& {beta_v_prime}\n'
#     omega_v: int = min(winner[2], winner[3]*winner[6]*omega_v_prime)
#     beta_v: int = winner[3]*min(winner[2], winner[6], omega_v_prime)*winner[7]*beta_v_prime
#     file_output += f'omega_v =& {omega_v}\n'
#     file_output += f'beta_v =& {beta_v}\n'
#     omega: int = min(winner[2], 2*omega_v+2*winner[6]*omega_v_prime)
#     beta: int = 2*beta_v + 2*min(winner[2], 2*winner[6], omega_v_prime)*winner[7]*beta_v_prime
#     file_output += f'omega =& {omega}\n'
#     file_output += f'beta =& {beta}\n'
#     file_output += f'We have verification keys of size 2dlog(p)={2*winner[2]*ceil(log2(winner[1]))} bits, which are {2*winner[2]*ceil(log2(winner[1]))/8000} kilobytes. \n'
#     file_output += f'We have signatures that are ell*d*ceil(log(2*beta_v_prime+1))={winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))} bits, which are {winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))/8000} kilobytes. \n'
#     file_output += f'However, K={winner[3]} of these aggregate together. An aggregate signature is ell*d*ceil(log(2*beta_v+1))={winner[10]*winner[2]*ceil(log2(2*beta_v+1))} bits, which is {winner[10]*winner[2]*ceil(log2(2*beta_v+1))/8000} kilobytes. \n'
#     file_output += f'Together, K verification keys and one aggregate signature takes {winner[3]*(2*winner[2]*log2(winner[1])/8000) + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/8000} kilobytes. \n'
#     file_output += f'This cost in space averages out to {2*winner[2]*log2(winner[1])/8000 + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/(winner[3]*8000)} kilobytes per signer.\n'
#     if winner[-1] <= CRYSTALS_TARGET:
#         file_output += f'Together, Fusion keys and aggregate signatures have better space efficiency than the smallest CRYSTALS-Dilithium signatures, which are {CRYSTALS_TARGET} kilobytes. \n'
#         ct: int = 2
#         while 2*winner[2]*log2(winner[1])/8000 + winner[10]*winner[2]*ceil(log2(2*beta_v+1))/(ct*8000) > CRYSTALS_TARGET and ct < winner[3]:
#             ct += 1
#         if ct < winner[3]:
#             file_output += f'In fact, if we do not aggregate at full capacity, then aggregating only {ct} signatures (i.e. occupying {100*ct/winner[3]}\% of signature capacity) is sufficient to beat CRYSTALS-Dilithium efficiency.\n '
#     if winner[-1] <= FALCON_TARGET:
#         file_output += f'Together, Fusion keys and aggregate signatures have better space efficiency than the smallest Falcon signatures, which are {FALCON_TARGET} kilobytes.\n '
#         ct: int = 2
#         while 2 * winner[2] * log2(winner[1]) / 8000 + winner[10] * winner[2] * ceil(log2(2 * beta_v + 1)) / (ct * 8000) > FALCON_TARGET and ct < winner[3]:
#             ct += 1
#         if ct < winner[3]:
#             file_output += f'In fact, if we do not aggregate at full capacity, then aggregating only {ct} signatures (i.e. occupying {100 * ct / winner[3]}\% of signature capacity) is sufficient to beat FALCON efficiency.\n '
#     file_output += f'The cost of using this scheme occurs in communication overhead: using these parameters, an aggregator must download K={winner[3]} public verification keys, signatures, and messages. \n'
#     file_output += f'Disregarding the cost of messages, this means downloading {winner[3]*2*winner[2]*log2(winner[1])/8000 + winner[3]*winner[10]*winner[2]*ceil(log2(2*beta_v_prime+1))/8000} kilobytes.\n'
#
#     print(f'secpar={winner[0]}, p={winner[1]}, d={winner[2]}, capacity={winner[3]}, omega_ch={winner[4]}, beta_ch={winner[5]}, omega_ag={winner[6]}, beta_ag={winner[7]}, omega_sk={winner[8]}, beta_sk={winner[9]}, ell={winner[10]}, avg_size={winner[11]}')
#     next_secpar_d_p=(winner[0], winner[2], winner[1])
#     with open(f"results{next_secpar_d_p}.txt", "w") as wf:
#         wf.write("secpar,p,d,capacity,omega_ch,beta_ch,omega_ag,beta_ag,omega_sk,beta_sk,ell,avg_weight_kb\n" + str(file_output) + "\n\n\n")
#
