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
ALLOWABLE_SECPARS: list[int] = [128, 256, 512, 1024]
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

degree_input: str = input(f'Please input a polynomial degree. d = ')
d: int = int(degree_input)
if not is_pow_two(x=d):
    raise ValueError

# find omega_ch, beta_ch that minimizes the product subject to the constraint that
# -log2binom(d=d, k=omega_ch) - omega_ch*log2(2*beta_ch+1) < -secpar
print(f'Finding omega_ch, beta_ch.')
found_pairs_of_beta_omega_ch: dict[tuple[int, int], float] = {}
minval: float = -1.0
argmin: tuple[int, int] = tuple()
for next_beta in range(1, 1024):
    for next_omega in range(1, d+1):
        if -log2binom(d=d, k=next_omega) - next_omega*log2(2*next_beta+1) < -secpar:
            found_pairs_of_beta_omega_ch[(next_beta, next_omega)] = next_omega*next_beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ch[(next_beta, next_omega)] <= minval:
                argmin = (next_beta, next_omega)
                minval = found_pairs_of_beta_omega_ch[(next_beta, next_omega)]

beta_ch: int = argmin[0]
omega_ch: int = argmin[1]
log2_epsilon_ch: float = -log2binom(d=d, k=omega_ch) - omega_ch*log2(2*beta_ch+1)

print(f'beta_ch={beta_ch}')
print(f'omega_ch={omega_ch}')
print(f'epsilon_ch={log2_epsilon_ch}')
print(f'Finding omega_ag, beta_ag.')

# find omega_ag, beta_ag that minimizes the product subject to the constraints that
# -log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1) < -1 and
# -log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1)
# -log2(1-epsilon_ch) -log2(1-2**(-log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1))) < -(secpar+1)
found_pairs_of_beta_omega_ag: dict[tuple[int, int], float] = {}
minval: float = -1.0
argmin: tuple[int, int] = tuple()
for next_beta in range(1, 1024):
    for next_omega in range(1, d+1):
        log2_epsilon_tmp: float = -log2binom(d=d, k=next_omega) - next_omega * log2(2 * next_beta + 1)
        if -log2binom(d=d, k=next_omega) - next_omega*log2(2*next_beta+1) < -1 and \
                log2_epsilon_tmp - log2(1-2**log2_epsilon_ch) - log2(1-2**log2_epsilon_tmp) < -(secpar+1):
            found_pairs_of_beta_omega_ch[(next_beta, next_omega)] = next_omega*next_beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ch[(next_beta, next_omega)] <= minval:
                argmin = (next_beta, next_omega)
                minval = found_pairs_of_beta_omega_ch[(next_beta, next_omega)]

beta_ag: int = argmin[0]
omega_ag: int = argmin[1]
log2_epsilon_ag: float = -log2binom(d=d, k=omega_ag) - omega_ag*log2(2*beta_ag+1)

print(f'beta_ag={beta_ag}')
print(f'omega_ag={omega_ag}')
print(f'epsilon_ag={log2_epsilon_ag}')

beta_v_prime_over_beta_sk: int = (1+min(d, omega_ch)*beta_ch)
beta_v_over_beta_sk: int = capacity*min(d,omega_ag)*beta_ag*beta_v_prime_over_beta_sk
beta_over_beta_sk: int = 2*(capacity*min(d,omega_ag)+min(d,2*omega_ag))*beta_ag*beta_v_prime_over_beta_sk

lower_bound_on_log_p_from_beta_without_beta_sk: float = log2(beta_over_beta_sk)+1
lower_bound_on_log_p_from_highlander_without_beta_sk: float = log2(8*min(d,2*omega_ag)*min(d,2*omega_ch)*beta_ag*beta_ch)+1
lower_bound_on_log_beta_sk_from_enough_keys: float = log2((min(d, omega_ch) * beta_ch - 1) / 2)
weak_lower_bound_on_p: float = max(lower_bound_on_log_p_from_highlander_without_beta_sk, lower_bound_on_log_p_from_beta_without_beta_sk) + lower_bound_on_log_beta_sk_from_enough_keys
candidate_p: int = ceil(2**weak_lower_bound_on_p)
if ((candidate_p - 1) % (2*d)) != 0:
    candidate_p += 2*d - ((candidate_p - 1) % (2*d))
assert ((candidate_p - 1) % 2*d) == 0
while not is_prime(x=candidate_p):
    candidate_p += 2*d
assert is_prime(x=candidate_p)

max_ell: int = 2**10
discovered_param_sets: dict[tuple, float] = {}

print(f'Beginning main loop.')

winner_efficiency = -1.0
winning_params = tuple()
while log2(candidate_p) <= 31:
    print(candidate_p)
    while not is_prime(x=candidate_p):
        candidate_p += 2*d
        print(candidate_p)
    # print(f'First finding bounds on beta_sk.')
    upper_bound_on_log_beta_sk: float = log2((candidate_p - 1) / 2) - log2(beta_over_beta_sk)
    upper_bound_on_beta_sk: int = floor(2**upper_bound_on_log_beta_sk)
    lower_bound_on_beta_sk: int = ceil(2**lower_bound_on_log_beta_sk_from_enough_keys)
    while upper_bound_on_beta_sk <= lower_bound_on_beta_sk:
        candidate_p += 2*d
        while not is_prime(x=candidate_p):
            candidate_p += 2 * d
        upper_bound_on_log_beta_sk: float = log2((candidate_p - 1) / 2) - log2(beta_over_beta_sk)
        upper_bound_on_beta_sk: int = floor(2**upper_bound_on_log_beta_sk)
        print(candidate_p)

    # print(f'Now searching beta_sk range for all minimal ell that satisfies both constraints.')
    found_ell_beta_sk_pairs: list[tuple[int, int]] = []
    for next_beta_sk in range(lower_bound_on_beta_sk, upper_bound_on_beta_sk+1):
        # print(f'Searching for next_beta_sk={next_beta_sk}.')
        next_ell: int = 1
        next_beta: int = beta_over_beta_sk*next_beta_sk
        next_beta_v_prime: int = beta_v_prime_over_beta_sk*next_beta_sk

        # print(f'Checking next_ell={next_ell}.')
        security_constraint: float = log2delta - (log2(d)/2 + log2(next_beta) - log2(candidate_p)/next_ell)/(next_ell*d-1)
        enough_keys_constraint: float = (2*next_ell*d*log2(2*next_beta_sk+1)) - (secpar + 2*d*log2(candidate_p) + next_ell*d*log2(2*next_beta_v_prime+1))
        while (enough_keys_constraint <= 0 or security_constraint <= 0) and next_ell <= max_ell:
            next_ell += 1
            security_constraint: float = log2delta - (log2(d)/2 + log2(next_beta) - log2(candidate_p)/next_ell)/(next_ell*d-1)
            enough_keys_constraint: float = (2*next_ell*d*log2(2*next_beta_sk+1)) - (secpar + 2*d*log2(candidate_p) + next_ell*d*log2(2*next_beta_v_prime+1))
            # print(f'Checking next_ell={next_ell}.')

        if enough_keys_constraint > 0 and security_constraint > 0 and next_ell <= max_ell:
            found_ell_beta_sk_pairs += [(next_ell, next_beta_sk)]

    if len(found_ell_beta_sk_pairs) > 0:
        found_ell_beta_sk_pairs = sorted(found_ell_beta_sk_pairs, key=lambda x: x[0])

        ell = found_ell_beta_sk_pairs[0][0]
        beta_sk = found_ell_beta_sk_pairs[0][1]
        beta_v_prime: int = beta_v_prime_over_beta_sk*beta_sk
        beta_v: int = beta_v_over_beta_sk*beta_sk
        beta: int = beta_over_beta_sk*beta_sk

        # print(f'Checking for success.')
        success = (candidate_p-1) % (2*d) == 0
        success = success and beta < (candidate_p-1)/2
        success = success and 8*min(d,2*omega_ag)*min(d, 2*omega_ch)*beta_ag*beta_sk*beta_ch < (candidate_p-1)/2
        success = success and log2delta - (log2(sqrt(d)) + log2(beta) - log2(candidate_p)/ell)/(ell*d-1) > 0
        success = success and secpar + 2*d*log2(candidate_p) + ell*d*log2(2*beta_v_prime+1) <= 2*ell*d*log2(2*beta_sk+1)
        success = success and log2_epsilon_ch < -secpar
        success = success and log2_epsilon_ag < -1
        success = success and log2_epsilon_ag - log2(1-2**log2_epsilon_ag) - log2(1-2**log2_epsilon_ch) < -(secpar+1)
        success = success and 2*beta_sk+1 > min(d, omega_ch)*beta_ch

        if success:
            vk_size: int = 2 * d * ceil(log2(candidate_p))
            vk_size_in_kb: float = ceil(vk_size / 8) / 1000
            sig_size: int = ell * d * ceil(log2(2 * beta_v_prime + 1))
            sig_size_in_kb: float = ceil(sig_size / 8) / 1000
            agg_sig_size: int = ell * d * ceil(log2(2 * beta_v + 1))
            agg_sig_size_in_kb: float = ceil(agg_sig_size / 8) / 1000
            average_size_in_kb: float = agg_sig_size_in_kb / capacity + vk_size_in_kb

            if agg_sig_size_in_kb / capacity + vk_size_in_kb < CRYSTALS_TARGET:
                params = (secpar, capacity, d, ell, candidate_p, beta_sk, beta_ch, beta_ag, omega_ch, omega_ag, beta_v)
                discovered_param_sets[params] = average_size_in_kb
                print(f"Found a solution that beats CRYSTALS. (secpar, capacity, d, ell, p, beta_sk, beta_ch, beta_ag, omega_ch, omega_ag, beta_v) = {params} with efficiency {average_size_in_kb}")
                print("Total solutions found that beat CRYSTALS: " + str(len(discovered_param_sets)))
                with open("result.txt", "a") as wf:
                    wf.write(str(params) + "," + str(discovered_param_sets[params]) + "\n")
                if winner_efficiency < 0 or (winner_efficiency > 0 and discovered_param_sets[params] < winner_efficiency):
                    winning_params = params
                    winner_efficiency = discovered_param_sets[params]
                    with open("winnerd128.txt", "w") as wf:
                        wf.write(str(params) + "," + str(discovered_param_sets[params]))

    candidate_p += 2*d
