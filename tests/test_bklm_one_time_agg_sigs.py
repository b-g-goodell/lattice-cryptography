"""
We test the lattice-crypto.bklm_one_time_agg_sigs module.
"""
from copy import deepcopy
from random import shuffle
from secrets import randbits
from typing import Any, List, Dict
import pytest
from lattice_algebra import random_polynomialvector, random_polynomial
from lattice_cryptography.one_time_keys import ALLOWABLE_SECPARS, SecretSeed, OneTimeSigningKey, OneTimeVerificationKey, \
    SchemeParameters, bits_per_coefficient, bits_per_index_set, challenge_core as make_signature_challenge
from lattice_cryptography.bklm_one_time_agg_sigs import make_setup_parameters, keygen, sign, verify, \
    prepare_make_agg_coefs, prepare_hash2polyinput, make_agg_coefs, prepare_aggregate_input, aggregate, \
    make_aggregate_verify_target, aggregate_verify, LPs, BDs, WTs, SALTs, CAPs, DISTRIBUTION

SAMPLE_SIZE: int = 2 ** 6
REQUIRED_SETUP_PARAMETERS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'vf_bd',
                             'vf_wt', 'ag_bd', 'ag_wt', 'ag_salt', 'avf_bd', 'avf_wt', 'ag_cap']
BOUND_NAMES = ['sk_bd', 'ch_bd', 'vf_bd', 'ag_bd', 'avf_bd']
WEIGHT_NAMES = ['sk_wt', 'ch_wt', 'vf_wt', 'ag_wt', 'avf_wt']
SALT_NAMES = ['sk_salt', 'ch_salt', 'ag_salt']

# i = 0 : secpar
# i = 1 : SchemeParameters
# i = 2 : PublicParameters
SOME_SECPARS_AND_SPS = [(secpar, SchemeParameters(secpar=secpar, lp=LPs[secpar], distribution=DISTRIBUTION)) for secpar
                        in ALLOWABLE_SECPARS]
SOME_PUBLIC_PARAMETERS = [i + tuple([{
    'scheme_parameters': i[-1],
    'sk_salt': SALTs[i[0]]['sk_salt'],
    'sk_bd': BDs[i[0]]['sk_bd'],
    'sk_wt': WTs[i[0]]['sk_wt'],
    'ch_salt': SALTs[i[0]]['ch_salt'],
    'ch_bd': BDs[i[0]]['ch_bd'],
    'ch_wt': WTs[i[0]]['ch_wt'],
    'vf_bd': max(1, min(i[-1].lp.modulus // 2,
                        BDs[i[0]]['sk_bd'] * (1 + BDs[i[0]]['ch_bd'] * min(WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt'])))),
    'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
    'ag_cap': CAPs[i[0]],
    'ag_salt': SALTs[i[0]]['ag_salt'],
    'ag_bd': BDs[i[0]]['ag_bd'],
    'ag_wt': WTs[i[0]]['ag_wt'],
    'avf_bd': max(1, min(i[-1].lp.modulus // 2, CAPs[i[0]] * min(WTs[i[0]]['ag_wt'], max(1, min(i[-1].lp.degree,
                                                                                                WTs[i[0]]['sk_wt'] * (
                                                                                                        1 +
                                                                                                        WTs[i[0]][
                                                                                                            'ch_wt'])))) *
                         BDs[i[0]]['ag_bd'] * max(1, min(i[-1].lp.modulus // 2, BDs[i[0]]['sk_bd'] * (
            1 + BDs[i[0]]['ch_bd'] * min(WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt'])))))),
    'avf_wt': max(1, min(i[-1].lp.degree, CAPs[i[0]] * WTs[i[0]]['ag_wt'] * max(1, min(i[-1].lp.degree,
                                                                                       WTs[i[0]]['sk_wt'] * (
                                                                                               1 + WTs[i[0]][
                                                                                           'ch_wt']))))),
}]) for i in SOME_SECPARS_AND_SPS]
MAKE_RANDOM_SEED_CASES = [
    i + tuple([j, bin(j)[2:].zfill(i[0]), SecretSeed(secpar=i[0], lp=i[1].lp, seed=bin(j)[2:].zfill(i[0]))]) for i in
    SOME_PUBLIC_PARAMETERS for j in range(2 * SAMPLE_SIZE + 1)
]
BD_MIN: int = 1
BD_MAX: int = 1
SOME_BDS: List[int] = list(range(BD_MIN, BD_MAX + 1))

WT_MIN: int = 1
WT_MAX: int = 1
SOME_WTS: List[int] = list(range(WT_MIN, WT_MAX + 1))
MAKE_ONE_KEY_CASES = [i + tuple([
    sk_bd,
    sk_wt,
    random_polynomialvector(
        secpar=i[0],
        lp=i[1].lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': sk_bd, 'wt': sk_wt},
        num_coefs=sk_wt,
        bti=bits_per_index_set(secpar=i[0], degree=i[1].lp.degree, wt=sk_wt),
        btd=bits_per_coefficient(secpar=i[0], bd=sk_bd),
        const_time_flag=False),
    random_polynomialvector(
        secpar=i[0],
        lp=i[1].lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': sk_bd, 'wt': sk_wt},
        num_coefs=sk_wt,
        bti=bits_per_index_set(secpar=i[0], degree=i[1].lp.degree, wt=sk_wt),
        btd=bits_per_coefficient(secpar=i[0], bd=sk_bd),
        const_time_flag=False),
]) for i in MAKE_RANDOM_SEED_CASES for sk_bd in SOME_BDS for sk_wt in SOME_WTS]
for i in MAKE_ONE_KEY_CASES:
    i[8].const_time_flag = False
    i[9].const_time_flag = False
MAKE_ONE_KEY_CASES = [i + tuple([
    OneTimeSigningKey(
        secpar=i[0],
        lp=i[1].lp,
        left_key=i[-2],
        right_key=i[-1],
        const_time_flag=False)
]) for i in MAKE_ONE_KEY_CASES]
for i in MAKE_ONE_KEY_CASES:
    i[-1][0].const_time_flag = False
    i[-1][1].const_time_flag = False
MAKE_ONE_KEY_CASES = [i + tuple([
    i[1].key_ch * i[-3],
    i[1].key_ch * i[-2]
]) for i in MAKE_ONE_KEY_CASES]
MAKE_ONE_KEY_CASES = [i + tuple([
    OneTimeVerificationKey(secpar=i[0], lp=i[1].lp, left_key=i[-2], right_key=i[-1])
]) for i in MAKE_ONE_KEY_CASES]
for i in MAKE_ONE_KEY_CASES:
    i[-1][0].const_time_flag = False
    i[-1][1].const_time_flag = False
MAKE_SIGNATURE_CHALLENGE_CASES = []
for i in MAKE_ONE_KEY_CASES:
    for j in SOME_BDS:
        for k in SOME_WTS:
            MAKE_SIGNATURE_CHALLENGE_CASES += [
                i + tuple([
                    j,
                    k,
                    bin(randbits(i[0]))[2:].zfill(i[0]),
                    random_polynomial(
                        secpar=i[0],
                        lp=i[1].lp,
                        distribution=DISTRIBUTION,
                        dist_pars={'bd': j, 'wt': k},
                        num_coefs=k,
                        bti=bits_per_index_set(secpar=i[0], degree=i[1].lp.degree, wt=k),
                        btd=bits_per_coefficient(secpar=i[0], bd=j),
                        const_time_flag=False)])]
MAKE_SIGN_CASES = [
    i + tuple([
        i[8] ** i[-1] + i[9],
        i[10][0] ** i[-1] + i[10][1],
        True,
    ]) for i in MAKE_SIGNATURE_CHALLENGE_CASES
]
A_SIZE: int = 32
LIST_OF_LISTS_OF_UNQ_STRINGS = []
while len(LIST_OF_LISTS_OF_UNQ_STRINGS) < SAMPLE_SIZE:
    next_key_list = []
    while len(next_key_list) < SAMPLE_SIZE:
        next_key = bin(randbits(A_SIZE))[2:].zfill(A_SIZE)
        if next_key not in next_key_list:
            next_key_list += [next_key]
    LIST_OF_LISTS_OF_UNQ_STRINGS += [next_key_list]
PREPARE_MAKE_AGG_COEFS_CASES = [
    (i, [bin(randbits(A_SIZE))[2:].zfill(A_SIZE) for _ in range(SAMPLE_SIZE)]) for i in LIST_OF_LISTS_OF_UNQ_STRINGS
]
PREPARE_MAKE_AGG_COEFS_CASES = [
    i + tuple([
        list(zip(i[-2], i[-1]))
    ]) for i in PREPARE_MAKE_AGG_COEFS_CASES
]
PREPARE_MAKE_AGG_COEFS_CASES = [
    i + tuple([
        sorted(i[-1], key=lambda x: str(x[0]))
    ]) for i in PREPARE_MAKE_AGG_COEFS_CASES
]
PREPARE_MAKE_AGG_COEFS_CASES = [
    i + tuple([
        [j[0] for j in i[-1]],
        [j[1] for j in i[-1]]
    ]) for i in PREPARE_MAKE_AGG_COEFS_CASES
]


@pytest.mark.parametrize("unq_keys,ran_msgs,zipped_keys_and_msgs,srt_zipped_keys_and_msgs,exp_srt_keys,exp_srt_msgs",
                         PREPARE_MAKE_AGG_COEFS_CASES)
def test_prepare_make_agg_coefs(unq_keys, ran_msgs, zipped_keys_and_msgs, srt_zipped_keys_and_msgs, exp_srt_keys,
                                exp_srt_msgs):
    obs_srt_keys, obs_srt_msgs = prepare_make_agg_coefs(otvks=unq_keys, msgs=ran_msgs)
    assert all(next_key in unq_keys for next_key in obs_srt_keys)
    assert all(next_key in unq_keys for next_key in exp_srt_keys)
    assert all(next_msg in ran_msgs for next_msg in obs_srt_msgs)
    assert all(next_msg in ran_msgs for next_msg in exp_srt_msgs)
    assert len(obs_srt_keys) == len(unq_keys)
    assert len(obs_srt_msgs) == len(ran_msgs)
    assert obs_srt_keys == exp_srt_keys
    assert obs_srt_msgs == exp_srt_msgs
    shuffled_keys, shuffled_msgs = deepcopy(unq_keys), deepcopy(ran_msgs)
    shuffle(shuffled_keys)
    shuffle(shuffled_msgs)
    obs_srt_keys, obs_srt_msgs = prepare_make_agg_coefs(otvks=unq_keys, msgs=ran_msgs)
    assert all(next_key in unq_keys for next_key in obs_srt_keys)
    assert all(next_key in unq_keys for next_key in exp_srt_keys)
    assert all(next_msg in ran_msgs for next_msg in obs_srt_msgs)
    assert all(next_msg in ran_msgs for next_msg in exp_srt_msgs)
    assert len(obs_srt_keys) == len(unq_keys)
    assert len(obs_srt_msgs) == len(ran_msgs)
    assert obs_srt_keys == exp_srt_keys
    assert obs_srt_msgs == exp_srt_msgs


PREPARE_HASH2POLYINPUT_CASES = [i + j for i in SOME_PUBLIC_PARAMETERS for j in PREPARE_MAKE_AGG_COEFS_CASES]
PREPARE_HASH2POLYINPUT_CASES = [
    i + tuple([
        bits_per_coefficient(secpar=i[0], bd=i[2]['ag_bd']),
        bits_per_index_set(secpar=i[0], degree=i[1].lp.degree, wt=i[2]['ag_wt']),
        str(zip(i[7], i[8])),
        {
            'secpar': i[1].secpar,
            'lp': i[1].lp,
            'distribution': i[1].distribution,
            'dist_pars': {
                'bd': i[2]['ag_bd'],
                'wt': i[2]['ag_wt']},
            'num_coefs': i[2]['ag_wt'],
            'bti': bits_per_index_set(secpar=i[0], degree=i[1].lp.degree, wt=i[2]['ag_wt']),
            'btd': bits_per_coefficient(secpar=i[0], bd=i[2]['ag_bd']),
            'msg': str(list(zip(i[7], i[8]))),
            'const_time_flag': False
        }
    ]) for i in PREPARE_HASH2POLYINPUT_CASES]


@pytest.mark.parametrize(
    "secpar,sp,pp,unq_keys,ran_msgs,zipped_keys_and_msgs,srt_zipped_keys_and_msgs,exp_srt_keys,exp_srt_msgs,btd,bti,msg,exp_h2pinput",
    PREPARE_HASH2POLYINPUT_CASES)
def test_prepare_hash2polyinput(secpar, sp, pp, unq_keys, ran_msgs, zipped_keys_and_msgs, srt_zipped_keys_and_msgs,
                                exp_srt_keys, exp_srt_msgs, btd, bti, msg, exp_h2pinput):
    obs_h2pinput = prepare_hash2polyinput(pp=pp, otvks=unq_keys, msgs=ran_msgs)
    assert obs_h2pinput == exp_h2pinput


@pytest.mark.parametrize("secpar,sp,pp", SOME_PUBLIC_PARAMETERS)
def test_all(secpar, sp, pp):
    for i in range(SAMPLE_SIZE):
        # Sample some new one-time keys
        some_signing_keys = keygen(pp=pp, num_keys_to_gen=pp['ag_cap'])
        some_msgs = [bin(randbits(A_SIZE))[2:].zfill(A_SIZE) for _ in some_signing_keys]
        some_sigs = [sign(pp=pp, otk=i, msg=j) for i, j in zip(some_signing_keys, some_msgs)]
        assert all(verify(pp=pp, otvk=i[2], msg=j, sig=k) for i, j, k in zip(some_signing_keys, some_msgs, some_sigs))
        ag_sig = aggregate(pp=pp, otvks=[i[2] for i in some_signing_keys], msgs=some_msgs, sigs=some_sigs)
        assert aggregate_verify(pp=pp, otvks=[i[2] for i in some_signing_keys], msgs=some_msgs, ag_sig=ag_sig)
