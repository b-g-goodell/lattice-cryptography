"""
We test the lattice-crypto.lm_one_time_sigs module.
"""
import pytest
from lattice_cryptography.one_time_keys import \
    ALLOWABLE_SECPARS, \
    SchemeParameters
from lattice_cryptography.adaptor_sigs import \
    make_setup_parameters, \
    LPs, \
    SALTs, \
    BDs, \
    WTs, \
    DISTRIBUTION, \
    make_random_seed, \
    make_one_wit, \
    make_one_key, \
    witgen, \
    keygen, \
    make_signature_challenge, \
    presign, \
    preverify, \
    adapt, \
    extract, \
    witness_verify, \
    sign, \
    verify
from secrets import randbits

SAMPLE_SIZE: int = 2 ** 3
REQUIRED_SETUP_PARAMETERS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'wit_bd',
                             'wit_wt', 'wit_salt', 'vf_bd', 'vf_wt', 'pvf_bd', 'pvf_wt', 'ext_wit_bd', 'ext_wit_wt']
BOUND_NAMES = ['sk_bd', 'ch_bd', 'wit_bd', 'vf_bd', 'pvf_bd', 'ext_wit_bd']
WEIGHT_NAMES = ['sk_wt', 'ch_wt', 'wit_wt', 'vf_wt', 'pvf_wt', 'ext_wit_wt']
SALT_NAMES = ['sk_salt', 'ch_salt', 'wit_salt']

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
        'wit_salt': SALTs[i[0]]['wit_salt'],
        'wit_bd': BDs[i[0]]['wit_bd'],
        'wit_wt': WTs[i[0]]['wit_wt'],
        'pvf_bd': max(1, min((i[-1].lp.modulus - 1) // 2, BDs[i[0]]['sk_bd'] * (
                    1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']))),
        'pvf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
        'vf_bd': max(1, min((i[-1].lp.modulus - 1) // 2, BDs[i[0]]['sk_bd'] * (
                    1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']) + BDs[i[0]][
                                'wit_bd'])),
        'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']) + WTs[i[0]]['wit_wt'])),
        'ext_wit_bd': max(1, min((i[-1].lp.modulus - 1) // 2, BDs[i[0]]['sk_bd'] * (
                    1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']) + BDs[i[0]][
                                     'sk_bd'] * (1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) *
                                                 BDs[i[0]]['ch_bd']) + BDs[i[0]]['wit_bd'])),
        'ext_wit_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']) + WTs[i[0]]['sk_wt'] * (
                    1 + WTs[i[0]]['ch_wt']) + WTs[i[0]]['wit_wt'])),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES
]


@pytest.mark.parametrize("secpar,sp,expected_pp", MAKE_SETUP_PARAMETERS_CASES)
def test_make_setup_parameters(mocker, secpar, sp, expected_pp):
    mocker.patch('lattice_cryptography.one_time_keys.random_polynomialvector', return_value=sp.key_ch)
    assert expected_pp == make_setup_parameters(secpar)


@pytest.mark.parametrize("secpar,sp,pp", MAKE_SETUP_PARAMETERS_CASES)
def test_main(secpar, sp, pp):
    for j in ALLOWABLE_SECPARS:
        pp = make_setup_parameters(secpar=j)
        for i in range(SAMPLE_SIZE):
            some_key = keygen(pp=pp, num_keys_to_gen=1)[0]
            some_otvk = some_key[2]
            some_witpair = witgen(pp=pp, num_wits_to_gen=1)[0]
            some_wit = some_witpair[1]
            some_pubstat = some_witpair[2]
            some_msg = bin(randbits(j))[2:].zfill(j)
            some_presig = presign(pp=pp, otk=some_key, msg=some_msg, st=some_pubstat)
            assert preverify(pp=pp, otvk=some_otvk, msg=some_msg, st=some_pubstat, presig=some_presig)
            some_sig = adapt(presig=some_presig, wit=some_wit)
            assert verify(pp=pp, otvk=some_otvk, msg=some_msg, st=some_pubstat, sig=some_sig)
            ext_wit = extract(pp=pp, presig=some_presig, sig=some_sig)
            assert witness_verify(pp=pp, wit=ext_wit, st=some_pubstat)
