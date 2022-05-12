"""
We test the lattice-crypto.keys.one_time_keys module.
TODO: Fix this so it does not depend on lattice parameters and salts and stuff from lm_one_time_sigs
"""
import pytest
from lattice_algebra import is_ntt_friendly_prime, random_polynomial
from lattice_cryptography.one_time_keys import *
from lattice_cryptography.lm_one_time_sigs import LPs, SALTs, BDs, WTs
from secrets import randbits
from typing import List


SAMPLE_SIZE: int = 2 ** 0
MIN_LOG2D: int = 5
MAX_LOG2D: int = 7
SOME_LOG2D: List[int] = list(range(MIN_LOG2D, MAX_LOG2D+1))
MAX_Q: int = 2**(MAX_LOG2D+1)
SOME_NTT_FRIENDLY_PAIRS = [(2**log2d, q) for log2d in SOME_LOG2D for q in range(2**(log2d+1)+1, MAX_Q, 2**(log2d+1)) if is_ntt_friendly_prime(modulus=q, degree=2**log2d)]
MIN_LEN: int = 1
MAX_LEN: int = 2**2
SOME_LEN: List[int] = list(range(MIN_LEN, MAX_LEN))
SECRET_SEED_INIT_CASES = []


for secpar in ALLOWABLE_SECPARS:
    for (d, q) in SOME_NTT_FRIENDLY_PAIRS:
        for l in SOME_LEN:
            for _ in range(SAMPLE_SIZE):
                lp: LatticeParameters = LatticeParameters(modulus=q, degree=d, length=l)
                seed = bin(randbits(secpar))[2:].zfill(secpar)
                SECRET_SEED_INIT_CASES += [(secpar, lp, seed, SecretSeed(secpar=secpar, lp=lp, seed=seed)) ]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed", SECRET_SEED_INIT_CASES)
def test_secret_seed_init(secpar, lp, seed, secret_seed):
    assert secret_seed.secpar == secpar
    assert secret_seed.lp == lp
    assert secret_seed.seed == seed


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed", SECRET_SEED_INIT_CASES)
def test_secret_seed_eq(secpar, lp, seed, secret_seed):
    assert SecretSeed(secpar=secpar, lp=lp, seed=seed) == secret_seed


DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT
BD_MIN: int = 1
BD_MAX: int = 2**1
SOME_BDS: List[int] = list(range(BD_MIN, BD_MAX+1))
WT_MIN: int = 1
WT_MAX: int = 2**1
SOME_WTS: List[int] = list(range(WT_MIN, WT_MAX + 1))
ONE_TIME_SECRET_WITNESS_CASES = [
    i + tuple([bd, wt,
        OneTimeSecretWitness(
            secpar=i[0],
            lp=i[1],
            key=random_polynomialvector(
                secpar=i[0],
                lp=i[1],
                distribution=DISTRIBUTION,
                dist_pars={'bd': bd, 'wt': wt},
                num_coefs=wt,
                bti=bits_per_index_set(secpar=i[0], degree=i[1].degree, wt=wt),
                btd=bits_per_coefficient(secpar=i[0], bd=bd),
                const_time_flag=False,
            ))
    ]) for i in SECRET_SEED_INIT_CASES for bd in SOME_BDS for wt in SOME_WTS
]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,bd,wt,wit", ONE_TIME_SECRET_WITNESS_CASES)
def test_one_time_secret_witness_init(secpar, lp, seed, secret_seed, bd, wt, wit):
    assert isinstance(wit, OneTimeSecretWitness)
    assert wit.secpar == secret_seed.secpar == secpar
    assert wit.lp == secret_seed.lp == lp
    assert isinstance(wit.key, PolynomialVector)
    assert all(k.const_time_flag for k in wit.key.entries)
    cnw = wit.key.get_coef_rep()
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)
    assert n <= bd
    assert w <= wt


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,bd,wt,wit", ONE_TIME_SECRET_WITNESS_CASES)
def test_one_time_secret_witness_eq(secpar, lp, seed, secret_seed, bd, wt, wit):
    assert OneTimeSecretWitness(secpar=secpar, lp=lp, key=wit.key) == wit


KEY_CHALLENGES = [
    i + tuple([
        random_polynomialvector(
            secpar=i[0],
            lp=i[1],
            distribution=DISTRIBUTION,
            dist_pars={'bd': i[1].modulus//2, 'wt': i[1].degree},
            num_coefs=i[1].degree,
            bti=bits_per_index_set(
                secpar=i[0],
                degree=i[1].degree,
                wt=i[1].degree),
            btd=bits_per_coefficient(
                secpar=i[0],
                bd=i[1].modulus//2),
            const_time_flag=False)
    ]) for i in ONE_TIME_SECRET_WITNESS_CASES
]
ONE_TIME_PUBSTAT_CASES = [i + tuple([OneTimePublicStatement(secpar=i[0], lp=i[1], key=i[-1] * i[-2].key)]) for i in KEY_CHALLENGES]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,bd,wt,wit,key_ch,stat", ONE_TIME_PUBSTAT_CASES)
def test_one_time_pubstat_init(secpar, lp, seed, secret_seed, bd, wt, wit, key_ch, stat):
    assert isinstance(key_ch, PolynomialVector)
    assert isinstance(stat, OneTimePublicStatement)
    assert stat.secpar == secpar
    assert stat.lp == lp
    assert isinstance(stat.key, Polynomial)
    assert key_ch * wit.key == stat.key


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,bd,wt,wit,key_ch,stat", ONE_TIME_PUBSTAT_CASES)
def test_one_time_pubstat_eq(secpar, lp, seed, secret_seed, bd, wt, wit, key_ch, stat):
    assert stat == OneTimePublicStatement(secpar=secpar, lp=lp, key=key_ch * wit.key)


ONE_TIME_SIGNING_KEY_CASES = [i + tuple([
    bd,
    wt,
    random_polynomialvector(
        secpar=i[0],
        lp=i[1],
        distribution=DISTRIBUTION,
        dist_pars={'bd': i[1].modulus//2, 'wt': i[1].degree},
        num_coefs=i[1].degree,
        bti=bits_per_index_set(
            secpar=i[0],
            degree=i[1].degree,
            wt=i[1].degree),
        btd=bits_per_coefficient(
            secpar=i[0],
            bd=i[1].modulus//2),
        const_time_flag=False),
    random_polynomialvector(
        secpar=i[0],
        lp=i[1],
        distribution=DISTRIBUTION,
        dist_pars={'bd': i[1].modulus // 2, 'wt': i[1].degree},
        num_coefs=i[1].degree,
        bti=bits_per_index_set(
            secpar=i[0],
            degree=i[1].degree,
            wt=i[1].degree),
        btd=bits_per_coefficient(
            secpar=i[0],
            bd=i[1].modulus // 2),
        const_time_flag=False),
]) for i in ONE_TIME_PUBSTAT_CASES for bd in SOME_BDS for wt in SOME_WTS]
ONE_TIME_SIGNING_KEY_CASES = [i + tuple([
    OneTimeSigningKey(
        secpar=i[0],
        lp=i[1],
        left_key=i[-2],
        right_key=i[-1]
    )
]) for i in ONE_TIME_SIGNING_KEY_CASES]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk", ONE_TIME_SIGNING_KEY_CASES)
def test_one_time_signing_key_init(secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk):
    assert isinstance(sk, OneTimeSigningKey)
    assert sk.secpar == secpar
    assert sk.lp == lp
    assert sk.left_key == left_sk
    assert sk.right_key == right_sk


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk", ONE_TIME_SIGNING_KEY_CASES)
def test_one_time_signing_key_eq(secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk):
    x = OneTimeSigningKey(secpar=secpar, lp=lp, left_key=left_sk, right_key=right_sk)
    assert x == sk

VERIFICATION_KEY_CASES = [
    i + tuple([
        i[-7] * i[-1].left_key,
        i[-7] * i[-1].right_key
    ]) for i in ONE_TIME_SIGNING_KEY_CASES
]
VERIFICATION_KEY_CASES = [
    i + tuple([
        OneTimeVerificationKey(
            secpar=i[0],
            lp=i[1],
            left_key=i[-2],
            right_key=i[-1]
        )
    ]) for i in VERIFICATION_KEY_CASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk", VERIFICATION_KEY_CASES)
def test_one_time_vf_key_init(secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk):
    assert isinstance(vk, OneTimeVerificationKey)
    assert vk.secpar == secpar
    assert vk.lp == lp
    assert vk.left_key == left_vk
    assert vk.right_key == right_vk
    assert vk.left_key == key_ch * left_sk
    assert vk.right_key == key_ch * right_sk


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk", VERIFICATION_KEY_CASES)
def test_one_time_vf_key_getitem(secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk):
    assert vk[0] == left_vk == vk.left_key
    assert vk[1] == right_vk == vk.right_key


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk", VERIFICATION_KEY_CASES)
def test_one_time_vf_key_eq(secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk):
    assert vk == OneTimeVerificationKey(secpar=secpar, lp=lp, left_key=left_vk, right_key=right_vk)


SCHEME_PARAMETERS_CASES = [
    i + tuple([UNIFORM_INFINITY_WEIGHT]) for i in VERIFICATION_KEY_CASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,lp,seed,secret_seed,wit_bd,wit_wt,wit,key_ch,stat,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,distn", SCHEME_PARAMETERS_CASES)
def test_schemeparameters_init(mocker, secpar, lp, seed, secret_seed, wit_bd, wit_wt, wit, key_ch, stat, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk, distn):
    mocker.patch('lattice_cryptography.one_time_keys.random_polynomialvector', return_value=key_ch)
    sp = SchemeParameters(secpar=secpar, lp=lp, distribution=distn, key_ch=None)
    assert sp.secpar == secpar
    assert sp.lp == lp
    assert sp.distribution == distn
    assert sp.key_ch == key_ch


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
        'vf_bd': max(1, min(
            i[-1].lp.modulus // 2,
            BDs[i[0]]['sk_bd'] * (1 + min(WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']))),
        'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES
]
MAKE_RANDOM_SEED_CASES = [
    i + tuple([j, bin(j)[2:].zfill(i[0])]) for i in MAKE_SETUP_PARAMETERS_CASES for j in range(2*SAMPLE_SIZE+1)
]
MAKE_RANDOM_SEED_CASES = [
    i + tuple([SecretSeed(secpar=i[0], lp=i[1].lp, seed=i[-1])]) for i in MAKE_RANDOM_SEED_CASES
]


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,sp,pp,expected_int,expected_str,expected_seed", MAKE_RANDOM_SEED_CASES)
def test_make_random_seed(mocker, secpar, sp, pp, expected_int, expected_str, expected_seed):
    mocker.patch('lattice_cryptography.one_time_keys.randbelow', return_value=expected_int)
    observed_seed = make_random_seed(secpar=secpar, pp=pp)
    assert observed_seed == expected_seed, f'Needed observed_seed == expected_seed, but had observed_seed = {str(observed_seed)} and expected_seed = {str(expected_seed)}'
    assert observed_seed.seed == expected_str


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
        right_key=i[-1]
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


# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk", MAKE_ONE_KEY_CASES)
def test_make_one_key(mocker, secpar, sp, pp: Dict[str, Any], j, seed, secret_seed, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk):
    mocker.patch('lattice_cryptography.one_time_keys.make_random_seed', return_value=secret_seed)
    mocker.patch('lattice_cryptography.one_time_keys.hash2polynomialvector', side_effect=[left_sk, right_sk])

    observed_key_tuple = make_one_key(pp=pp)
    observed_seed, observed_sk, observed_vk = observed_key_tuple
    assert observed_seed == secret_seed
    assert observed_sk == sk
    assert sk[0] == left_sk
    assert sk[1] == right_sk
    assert observed_sk[0].const_time_flag
    assert observed_sk[1].const_time_flag

    assert not observed_vk[0].const_time_flag
    assert not observed_vk[1].const_time_flag
    assert observed_vk.left_key == pp['scheme_parameters'].key_ch * sk.left_key
    assert observed_vk.right_key == pp['scheme_parameters'].key_ch * sk.right_key

    assert vk == OneTimeVerificationKey(secpar=secpar, lp=pp['scheme_parameters'].lp, left_key=pp['scheme_parameters'].key_ch * sk.left_key, right_key=pp['scheme_parameters'].key_ch * sk.right_key)
    cnw = sk[0].get_coef_rep() + sk[1].get_coef_rep()
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)
    assert 1 <= n <= pp['sk_bd']
    assert 1 <= w <= pp['sk_wt']

# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk", MAKE_ONE_KEY_CASES)
def test_keygen_core(mocker, secpar, sp, pp: Dict[str, Any], j, seed, secret_seed, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk):
    mocker.patch('lattice_cryptography.one_time_keys.make_one_key', return_value=(secret_seed, sk, vk))
    key_tuples = keygen_core(pp=pp, num_keys_to_gen=1, seeds=[secret_seed])
    assert isinstance(key_tuples, list)
    assert len(key_tuples) == 1
    for next_key_tuple in key_tuples:
        assert isinstance(next_key_tuple, tuple)
        assert next_key_tuple[0] == secret_seed
        assert next_key_tuple[1] == sk
        assert next_key_tuple[2] == vk
        assert sp.key_ch == pp['scheme_parameters'].key_ch
        assert sp.key_ch * next_key_tuple[1][0] == vk[0] == vk.left_key
        assert sp.key_ch * next_key_tuple[1][1] == vk[1] == vk.right_key


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

# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch", MAKE_SIGNATURE_CHALLENGE_CASES)
def test_challenge_core(mocker, secpar, sp, pp, j, seed, secret_seed, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk, ch_bd, ch_wt, msg, sig_ch):
    mocker.patch('lattice_cryptography.one_time_keys.hash2polynomial', return_value=sig_ch)
    assert challenge_core(pp=pp, otvk=vk, msg=msg) == sig_ch


MAKE_SIGN_CASES = [
    i + tuple([
        i[8] ** i[-1] + i[9],
        True,
        i[10][0] ** i[-1] + i[10][1],
        True
    ]) for i in MAKE_SIGNATURE_CHALLENGE_CASES
]

# @pytest.mark.skip()
@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch,expected_sig_a,expected_valid_a,expected_sig_b,expected_valid_b", MAKE_SIGN_CASES)
def test_sign_core(mocker, secpar, sp, pp, j, seed, secret_seed, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk, ch_bd, ch_wt, msg, sig_ch, expected_sig_a, expected_valid_a, expected_sig_b, expected_valid_b):
    mocker.patch('lattice_cryptography.one_time_keys.challenge_core', return_value=sig_ch)
    observed_sig = sign_core(pp=pp, otk=(secret_seed, sk, vk), msg=msg)
    assert expected_sig_a == expected_sig_b
    assert sp.key_ch == pp['scheme_parameters'].key_ch
    assert sp.key_ch * sk[0] == sp.key_ch * sk.left_key == vk[0] == vk.left_key
    assert sp.key_ch * sk[1] == sp.key_ch * sk.right_key == vk[1] == vk.right_key
    assert sk[0] ** sig_ch + sk[1] == expected_sig_a
    assert observed_sig == sk[0] ** sig_ch + sk[1]
    assert sp.key_ch * observed_sig == vk[0] * sig_ch + vk[1]
    cnw = observed_sig.get_coef_rep()
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)
    assert 1 <= n <= pp['vf_bd']
    assert 1 <= w <= pp['vf_wt']


@pytest.mark.parametrize("secpar,sp,pp,j,seed,secret_seed,sk_bd,sk_wt,left_sk,right_sk,sk,left_vk,right_vk,vk,ch_bd,ch_wt,msg,sig_ch,sig_a,valid_a,sig_b,valid_b", MAKE_SIGN_CASES)
def test_vf(mocker, secpar, sp, pp, j, seed, secret_seed, sk_bd, sk_wt, left_sk, right_sk, sk, left_vk, right_vk, vk, ch_bd, ch_wt, msg, sig_ch, sig_a, valid_a, sig_b, valid_b):
    mocker.patch('lattice_cryptography.one_time_keys.challenge_core', return_value=sig_ch)
    assert sig_a == sig_b
    assert sk[0] ** sig_ch + sk[1] == sig_a
    assert sp.key_ch == pp['scheme_parameters'].key_ch
    assert sp.key_ch * sk[0] == sp.key_ch * sk.left_key == vk[0] == vk.left_key
    assert sp.key_ch * sk[1] == sp.key_ch * sk.right_key == vk[1] == vk.right_key
    assert sp.key_ch * sig_a == vk[0] * sig_ch + vk[1]
    cnw = sig_a.get_coef_rep()
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)
    assert 1 <= n <= pp['vf_bd']
    assert 1 <= w <= pp['vf_wt']
    assert verify_core(sig=sig_a, bd=pp['vf_bd'], wt=pp['vf_wt'], key_ch=pp['scheme_parameters'].key_ch, target=left_vk*sig_ch + right_vk) == valid_a
    assert verify_core(sig=sig_b, bd=pp['vf_bd'], wt=pp['vf_wt'], key_ch=pp['scheme_parameters'].key_ch, target=left_vk*sig_ch + right_vk) == valid_b
