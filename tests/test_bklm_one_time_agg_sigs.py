"""
We test the lattice-crypto.bklm_one_time_agg_sigs module.
"""
import pytest
from secrets import randbits
from lattice_cryptography.bklm_one_time_agg_sigs import \
    MODULI, \
    DEGREES, \
    LENGTHS, \
    SALTs, \
    BDs, \
    WTs, \
    ALLOWABLE_PARAMETERS, \
    DISTRIBUTION as DSTN, \
    LatticeParameters, \
    SecretSeed, \
    OneTimeSigningKey, \
    OneTimeVerificationKey, \
    SchemeParameters as SchPar, \
    make_setup_parameters, \
    keygen, \
    sign, \
    aggregate, \
    verify

SAMPLE_SIZE: int = 2 ** 0
MULTIPROCESSING_YES_OR_NO: bool = False
SS = SecretSeed
OTSK = OneTimeSigningKey
OTVK = OneTimeVerificationKey
OTK = tuple[SS, OTSK, OTVK]
REQ_SETUP_PARS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'vf_bd', 'vf_wt']
BOUND_NAMES = ['sk_bd', 'ch_bd', 'vf_bd']
WEIGHT_NAMES = ['sk_wt', 'ch_wt', 'vf_wt']
SALT_NAMES = ['sk_salt', 'ch_salt']

TARGET_SECPAR: int = 2 ** 8
TARGET_CAPACITY: int = 2 ** 11

MAKE_SCHEME_PARAMETERS_CASES = [(
    secpar,
    capacity,
    SchPar(secpar=secpar, lp=LatticeParameters(modulus=MODULI[(secpar, capacity)], degree=DEGREES[(secpar, capacity)],
                                               length=LENGTHS[(secpar, capacity)]), distribution=DSTN)
) for secpar, capacity in ALLOWABLE_PARAMETERS if secpar == TARGET_SECPAR and capacity == TARGET_CAPACITY]
MAKE_SETUP_PARAMETERS_CASES = [
    i + tuple([{
        'scheme_parameters': i[-1],
        'sk_salt': SALTs[(i[0], i[1])]['sk_salt'],
        'sk_bd': BDs[(i[0], i[1])]['sk_bd'],
        'sk_wt': WTs[(i[0], i[1])]['sk_wt'],
        'ch_salt': SALTs[(i[0], i[1])]['ch_salt'],
        'ch_bd': BDs[(i[0], i[1])]['ch_bd'],
        'ch_wt': WTs[(i[0], i[1])]['ch_wt'],
        'vf_bd': BDs[(i[0], i[1])]['sk_bd'] * (
                    1 + min(i[-1].lp.degree, WTs[(i[0], i[1])]['sk_wt'], WTs[(i[0], i[1])]['ch_wt']) *
                    BDs[(i[0], i[1])]['ch_bd']),
        'vf_wt': min(i[-1].lp.degree, WTs[(i[0], i[1])]['sk_wt'] * (1 + WTs[(i[0], i[1])]['ch_wt'])),
        'ag_bd': BDs[(i[0], i[1])]['ag_bd'],
        'ag_wt': WTs[(i[0], i[1])]['ag_wt'],
        'ag_cap': i[1],
        'ag_salt': SALTs[(i[0], i[1])]['ag_salt'],
        'avf_bd': i[1] * min(i[-1].lp.degree, WTs[(i[0], i[1])]['ag_wt'],
                             min(i[-1].lp.degree, WTs[(i[0], i[1])]['sk_wt'] * (1 + WTs[(i[0], i[1])]['ch_wt']))) *
                  BDs[(i[0], i[1])]['ag_bd'] * BDs[(i[0], i[1])]['sk_bd'] * (
                              1 + min(i[-1].lp.degree, WTs[(i[0], i[1])]['sk_wt'], WTs[(i[0], i[1])]['ch_wt']) *
                              BDs[(i[0], i[1])]['ch_bd']),
        'avf_wt': min(i[-1].lp.degree, i[1] * WTs[(i[0], i[1])]['ag_wt'] * WTs[(i[0], i[1])]['sk_wt'] * (
                    1 + WTs[(i[0], i[1])]['ch_wt'])),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES if i[0] == TARGET_SECPAR and i[1] == TARGET_CAPACITY
]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,cap,sp,expected_pp", MAKE_SETUP_PARAMETERS_CASES)
def test_make_setup_parameters(mocker, secpar, cap, sp, expected_pp):
    mocker.patch('lattice_cryptography.one_time_keys.random_polynomialvector', return_value=sp.key_ch)
    assert expected_pp == make_setup_parameters(secpar=secpar, carrying_capacity=cap)



MAKE_TEST_CORRECT_MAIN_PARAMETERS_CASES = list()
for i in MAKE_SETUP_PARAMETERS_CASES:
    if i[0] == TARGET_SECPAR and i[1] == TARGET_CAPACITY:
        for l in range(SAMPLE_SIZE):
            j = i  # secpar, cap, sp, pp
            j += tuple([keygen(pp=j[3], num=j[1], multiprocessing=MULTIPROCESSING_YES_OR_NO)])  # some_keys
            j += tuple([[k[2] for k in j[-1]]])  # some_otvks
            j += tuple([[bin(randbits(j[0]))[2:].zfill(j[0]) for _ in range(j[1])]])  # some_msgs
            j += tuple([[sign(pp=j[3], otk=otk, msg=msg) for otk, msg in zip(j[-3], j[-1])]])  # some_sigs
            j += tuple([[verify(pp=j[3], otvks=[next_otvk], msgs=[next_msg], sig=next_sig) for next_otvk, next_msg, next_sig
                         in zip(j[-3], j[-2], j[-1])]])  # some_valid_bits
            j += tuple([aggregate(pp=j[3], otvks=j[-4], msgs=j[-3], sigs=j[-2])])  # agg_sig
            j += tuple([verify(pp=j[3], otvks=j[-5], msgs=j[-4], sig=j[-1])])  # agg_valid_bit
            MAKE_TEST_CORRECT_MAIN_PARAMETERS_CASES += [j]


@pytest.mark.parametrize(
    "secpar,cap,sp,pp,some_keys,some_otvks,some_msgs,some_sigs,some_valid_bits,agg_sig,agg_valid_bit",
    MAKE_TEST_CORRECT_MAIN_PARAMETERS_CASES)
def test_correct_main(secpar, cap, sp, pp, some_keys, some_otvks, some_msgs, some_sigs, some_valid_bits, agg_sig,
                      agg_valid_bit):
    assert all(b for b in some_valid_bits) and agg_valid_bit
    for next_otvk, next_msg, next_sig, next_valid_bit in zip(some_otvks, some_msgs, some_sigs, some_valid_bits):
        assert verify(pp=pp, otvks=[next_otvk], msgs=[next_msg], sig=next_sig) == next_valid_bit
    assert verify(pp=pp, otvks=some_otvks, msgs=some_msgs, sig=agg_sig) == agg_valid_bit
