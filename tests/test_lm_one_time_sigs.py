"""
We test the lattice-crypto.lm_one_time_sigs module.
"""
import pytest
from lattice_algebra import random_polynomialvector, random_polynomial
from lattice_cryptography.one_time_keys import ALLOWABLE_SECPARS, SecretSeed, OneTimeSigningKey, OneTimeVerificationKey, SchemeParameters, bits_per_coefficient, bits_per_index_set, challenge_core as make_signature_challenge
from lattice_cryptography.lm_one_time_sigs import PublicParameters, Message, Challenge, Signature, OneTimeKeyTuple, make_setup_parameters, make_one_key, keygen, sign, verify, LPs, SALTs, BDs, WTs, DISTRIBUTION, make_random_seed
from secrets import randbits
from typing import Any, Dict, List


SAMPLE_SIZE: int = 2 ** 0
REQUIRED_SETUP_PARAMETERS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'vf_bd', 'vf_wt']
BOUND_NAMES = ['sk_bd', 'ch_bd', 'vf_bd']
WEIGHT_NAMES = ['sk_wt', 'ch_wt', 'vf_wt']
SALT_NAMES = ['sk_salt', 'ch_salt']
MAKE_SCHEME_PARAMETERS_CASES = [(
    secpar,
    SchemeParameters(secpar=secpar, lp=LPs[secpar], distribution=DISTRIBUTION)
) for secpar in ALLOWABLE_SECPARS]
MAKE_SETUP_PARAMETERS_CASES = [
    i + tuple([{
        'scheme_parameters': i[-1],
        'sk_salt': SALTs[i[0]]['sk_salt'],
        'sk_bd': BDs[i[0]]['sk_bd'],
        'sk_wt': WTs[i[0]]['sk_wt'],
        'ch_salt': SALTs[i[0]]['ch_salt'],
        'ch_bd': BDs[i[0]]['ch_bd'],
        'ch_wt': WTs[i[0]]['ch_wt'],
        'vf_bd': BDs[i[0]]['sk_bd'] * (1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']),
        'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES
]


@pytest.mark.parametrize("secpar,sp,expected_pp", MAKE_SETUP_PARAMETERS_CASES)
def test_make_setup_parameters(mocker, secpar, sp, expected_pp):
    mocker.patch('lattice_cryptography.one_time_keys.random_polynomialvector', return_value=sp.key_ch)
    assert expected_pp == make_setup_parameters(secpar)

MAKE_RANDOM_SEED_CASES = [
    i + tuple([j, bin(j)[2:].zfill(i[0])]) for i in MAKE_SETUP_PARAMETERS_CASES for j in range(2*SAMPLE_SIZE+1)
]
MAKE_RANDOM_SEED_CASES = [
    i + tuple([SecretSeed(secpar=i[0], lp=i[1].lp, seed=i[-1])]) for i in MAKE_RANDOM_SEED_CASES
]
BD_MIN: int = 1
BD_MAX: int = 2**1
SOME_BDS: List[int] = list(range(BD_MIN, BD_MAX+1))
WT_MIN: int = 1
WT_MAX: int = 2**1
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
        bti=bits_per_index_set(
            secpar=i[0],
            degree=i[1].lp.degree,
            wt=sk_wt),
        btd=bits_per_coefficient(
            secpar=i[0],
            bd=sk_bd),
        const_time_flag=False),
    random_polynomialvector(
        secpar=i[0],
        lp=i[1].lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': sk_bd, 'wt': sk_wt},
        num_coefs=sk_wt,
        bti=bits_per_index_set(
            secpar=i[0],
            degree=i[1].lp.degree,
            wt=sk_wt),
        btd=bits_per_coefficient(
            secpar=i[0],
            bd=sk_bd),
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
        const_time_flag=False
    )
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
    for ch_bd in SOME_BDS:
        for ch_wt in SOME_WTS:
            MAKE_SIGNATURE_CHALLENGE_CASES += [
                i + tuple([
                    ch_bd,
                    ch_wt,
                    bin(randbits(i[0]))[2:].zfill(i[0]),
                    random_polynomial(
                        secpar=i[0],
                        lp=i[1].lp,
                        distribution=DISTRIBUTION,
                        dist_pars={'bd': ch_bd, 'wt': ch_wt},
                        num_coefs=ch_wt,
                        bti=bits_per_index_set(
                            secpar=i[0],
                            degree=i[1].lp.degree,
                            wt=ch_wt),
                        btd=bits_per_coefficient(
                            secpar=i[0],
                            bd=ch_bd),
                        const_time_flag=False)])]
for i in MAKE_SIGNATURE_CHALLENGE_CASES:
    i[-1].const_time_flag = False
for i in MAKE_SIGNATURE_CHALLENGE_CASES:
    assert not i[8].const_time_flag
    assert not i[9].const_time_flag
    assert not i[10][0].const_time_flag
    assert not i[10][1].const_time_flag
    assert not i[10].left_key.const_time_flag
    assert not i[10].right_key.const_time_flag
    assert not i[11].const_time_flag
    assert not i[12].const_time_flag
    assert not i[13][0].const_time_flag
    assert not i[13][1].const_time_flag
    assert not i[13].left_key.const_time_flag
    assert not i[13].right_key.const_time_flag
    assert not i[17].const_time_flag
MAKE_SIGN_CASES = [
    i + tuple([
        i[8] ** i[-1] + i[9],
        True,
        i[10][0] ** i[-1] + i[10][1],
        True
    ]) for i in MAKE_SIGNATURE_CHALLENGE_CASES
]


@pytest.mark.parametrize("secpar,sp,pp", MAKE_SETUP_PARAMETERS_CASES)
def test_correct_general(secpar,sp,pp):
    some_keys = keygen(pp=pp, num_keys_to_gen=SAMPLE_SIZE, multiprocessing=True)
    some_msgs = [bin(randbits(secpar))[2:].zfill(secpar) for _ in range(SAMPLE_SIZE)]
    some_sigs = [sign(pp=pp, otk=otk, msg=msg) for otk, msg in zip(some_keys, some_msgs)]
    assert all(verify(pp=pp, otvk=otk[2], msg=msg, sig=sig) for otk, msg, sig in zip(some_keys, some_msgs, some_sigs))


@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch", MAKE_SIGNATURE_CHALLENGE_CASES)
def test_correct_specific(secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch):
    sig = sign(pp=pp, otk=(secret_seed, sk, vk), msg=msg)
    assert verify(pp=pp, otvk=vk, msg=msg, sig=sig)


@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch,sig_a,valid_a,sig_b,valid_b", MAKE_SIGN_CASES)
def test_verify(mocker,secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch,sig_a,valid_a,sig_b,valid_b):
    mocker.patch('lattice_cryptography.lm_one_time_sigs.make_challenge', return_value=sig_ch)
    assert sig_a == sig_b
    assert verify(pp=pp, otvk=vk, msg=msg, sig=sig_a)


