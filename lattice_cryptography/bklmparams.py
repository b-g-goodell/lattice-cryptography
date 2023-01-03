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

CRYSTALS_TARGET: float = 2.42/2
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
for beta in range(1, 1024):
    for omega in range(1, d+1):
        if -log2binom(d=d, k=omega) - omega*log2(2*beta+1) < -secpar:
            found_pairs_of_beta_omega_ch[(beta, omega)] = omega*beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ch[(beta, omega)] <= minval:
                argmin = (beta, omega)
                minval = found_pairs_of_beta_omega_ch[(beta, omega)]

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
for beta in range(1, 1024):
    for omega in range(1, d+1):
        log2_epsilon_tmp: float = -log2binom(d=d, k=omega) - omega * log2(2 * beta + 1)
        if -log2binom(d=d, k=omega) - omega*log2(2*beta+1) < -1 and \
                log2_epsilon_tmp - log2(1-2**log2_epsilon_ch) - log2(1-2**log2_epsilon_tmp) < -(secpar+1):
            found_pairs_of_beta_omega_ch[(beta, omega)] = omega*beta
            if minval == -1.0 or minval > 0.0 and found_pairs_of_beta_omega_ch[(beta, omega)] <= minval:
                argmin = (beta, omega)
                minval = found_pairs_of_beta_omega_ch[(beta, omega)]

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

print(f'Beginning main loop.')
start_with_prime_lower_bound: str = input(f'Prime lower bound? Hit enter to skip if no lower bound.')
if start_with_prime_lower_bound:
    candidate_p_lower_bound: int = int(start_with_prime_lower_bound)
    candidate_p: int = candidate_p_lower_bound
    while (candidate_p - 1) % 2*d != 0:
        candidate_p += 1
    while not is_prime(x=candidate_p):
        candidate_p += 2*d
else:
    candidate_p: int = 2*d + 1
    while log2(candidate_p) < weak_lower_bound_on_p or not is_prime(x=candidate_p):
        candidate_p += 2*d

success: bool = False
max_ell: int = 2**10
ct: int = 0
CT_MOD: int = 1000

while not success and log2(candidate_p) <= 31:
    while not is_prime(x=candidate_p):
        candidate_p += 2*d
    if ct % CT_MOD == 0:
        print(f'Trying candidate_p={candidate_p}.')
    ct += 1
    # print(f'First finding bounds on beta_sk.')
    upper_bound_on_log_beta_sk: float = log2((candidate_p - 1) / 2) - log2(beta_over_beta_sk)
    while upper_bound_on_log_beta_sk <= lower_bound_on_log_beta_sk_from_enough_keys or not is_prime(x=candidate_p):
        candidate_p += 2*d
        upper_bound_on_log_beta_sk: float = log2((candidate_p - 1) / 2) - log2(beta_over_beta_sk)
    upper_bound_on_beta_sk: int = floor(2**upper_bound_on_log_beta_sk)

    # print(f'Now searching beta_sk range for all minimal ell that satisfies both constraints.')
    found_ell_beta_sk_pairs: list[tuple[int, int]] = []
    for next_beta_sk in range(ceil(2**lower_bound_on_log_beta_sk_from_enough_keys), upper_bound_on_beta_sk+1):
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

    success = False
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

        vk_size: int = 2 * d * ceil(log2(candidate_p))
        vk_size_in_kb: float = ceil(vk_size / 8) / 1000
        sig_size: int = ell * d * ceil(log2(2 * beta_v_prime + 1))
        sig_size_in_kb: float = ceil(sig_size / 8) / 1000
        agg_sig_size: int = ell * d * ceil(log2(2 * beta_v + 1))
        agg_sig_size_in_kb: float = ceil(agg_sig_size / 8) / 1000
        average_size_in_kb: float = agg_sig_size_in_kb / capacity + vk_size_in_kb

        if success and agg_sig_size_in_kb / capacity + vk_size_in_kb < CRYSTALS_TARGET:
            print(f'Found sufficient parameters.\n secpar = {secpar}\n K = {capacity}\n d={d}\n beta_ch={beta_ch}\n beta_ag={beta_ag}\n omega_ag={omega_ag}\n omega_ch={omega_ch}\n p={candidate_p}\n ell={ell}\n beta_sk={beta_sk}\n beta_v_prime={beta_v_prime}\n beta_v={beta_v}\n beta={beta}.')
            print(f'For this parameter set, we have a public verification takes {vk_size} bits, or {vk_size_in_kb} kilobytes.')
            print(f'Non-aggregated signatures take {sig_size} bits, or {sig_size_in_kb} kilobytes.')
            print(f'Aggregated signatures cover K={capacity} signatures and take {agg_sig_size} bits, or {agg_sig_size_in_kb} kilobytes.')
            print(f'Hence, one aggregated signature and all K keys take {agg_sig_size_in_kb + capacity * vk_size_in_kb} kilobytes.')
            print(f'This averages out to {average_size_in_kb} kilobytes per signer.')
            print(f'This compares favorably to CD.')
        else:
            # print(f'Found sufficient parameters.\n secpar = {secpar}\n K = {capacity}\n d={d}\n beta_ch={beta_ch}\n beta_ag={beta_ag}\n omega_ag={omega_ag}\n omega_ch={omega_ch}\n p={candidate_p}\n ell={ell}\n beta_sk={beta_sk}\n beta_v_prime={beta_v_prime}\n beta_v={beta_v}\n beta={beta}.')
            # print(f'For this parameter set, we have a public verification takes {vk_size} bits, or {vk_size_in_kb} kilobytes.')
            # print(f'Non-aggregated signatures take {sig_size} bits, or {sig_size_in_kb} kilobytes.')
            # print(f'Aggregated signatures cover K={capacity} signatures and take {agg_sig_size} bits, or {agg_sig_size_in_kb} kilobytes.')
            # print(f'Hence, one aggregated signature and all K keys take {agg_sig_size_in_kb + capacity * vk_size_in_kb} kilobytes.')
            # print(f'This averages out to {average_size_in_kb} kilobytes per signer.')
            success = False

    if not success:
        # print(f'Failed, trying again.')
        candidate_p += 2*d

if not success:
    print('FAIL')

print(f'Found sufficient parameters.\n secpar = {secpar}\n K = {capacity}\n d={d}\n beta_ch={beta_ch}\n beta_ag={beta_ag}\n omega_ag={omega_ag}\n omega_ch={omega_ch}\n p={candidate_p}\n ell={ell}\n beta_sk={beta_sk}\n beta_v_prime={beta_v_prime}\n beta_v={beta_v}\n beta={beta}.')

vk_size: int = 2 * d * ceil(log2(candidate_p))
vk_size_in_kb: float = ceil(vk_size / 8) / 1000
sig_size: int = ell * d * ceil(log2(2 * beta_v_prime + 1))
sig_size_in_kb: float = ceil(sig_size / 8) / 1000
agg_sig_size: int = ell * d * ceil(log2(2 * beta_v + 1))
agg_sig_size_in_kb: float = ceil(agg_sig_size / 8) / 1000
average_size_in_kb: float = agg_sig_size_in_kb / capacity + vk_size_in_kb

print(f'For this parameter set, we have a public verification takes {vk_size} bits, or {vk_size_in_kb} kilobytes.')
print(f'Non-aggregated signatures take {sig_size} bits, or {sig_size_in_kb} kilobytes.')
print(f'Aggregated signatures cover K={capacity} signatures and take {agg_sig_size} bits, or {agg_sig_size_in_kb} kilobytes.')
print(f'Hence, one aggregated signature and all K keys take {agg_sig_size_in_kb + capacity * vk_size_in_kb} kilobytes.')
print(f'This averages out to {average_size_in_kb} kilobytes per signer.')
if agg_sig_size_in_kb / capacity + vk_size_in_kb < CRYSTALS_TARGET:
    print(f'This compares favorably to CD.')
else:
    print(f'This does not compare favorably to CD. Since this is not favorable compared to CD, we will move to the next prime and try again.')
    success = False