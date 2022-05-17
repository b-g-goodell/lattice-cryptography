"""
We test the lattice-crypto.bklm_one_time_agg_sigs module.
"""
import pytest
from secrets import randbits
from lattice_cryptography.one_time_keys import \
    ALLOWABLE_SECPARS, \
    SecretSeed, \
    OneTimeSigningKey, \
    OneTimeVerificationKey, \
    SchemeParameters as SchPar
from lattice_cryptography.bklm_one_time_agg_sigs import \
    LPs, \
    BDs, \
    WTs, \
    SALTs, \
    CAPs, \
    DISTRIBUTION as DSTN, \
    make_setup_parameters, \
    keygen, \
    sign, \
    aggregate, \
    verify

SAMPLE_SIZE: int = 2 ** 3
SS = SecretSeed
OTSK = OneTimeSigningKey
OTVK = OneTimeVerificationKey
OTK = tuple[SS, OTSK, OTVK]
REQ_SETUP_PARS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'vf_bd', 'vf_wt']
BOUND_NAMES = ['sk_bd', 'ch_bd', 'vf_bd']
WEIGHT_NAMES = ['sk_wt', 'ch_wt', 'vf_wt']
SALT_NAMES = ['sk_salt', 'ch_salt']
MAKE_SCHEME_PARAMETERS_CASES = [(
    secpar,
    SchPar(secpar=secpar, lp=LPs[secpar], distribution=DSTN)
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
        'ag_bd': BDs[i[0]]['ag_bd'],
        'ag_wt': WTs[i[0]]['ag_wt'],
        'ag_cap': CAPs[i[0]],
        'ag_salt': SALTs[i[0]]['ag_salt'],
        'avf_bd': CAPs[i[0]] * min(i[-1].lp.degree, WTs[i[0]]['ag_wt'], max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt'])))) * BDs[i[0]]['ag_bd'] * BDs[i[0]]['sk_bd'] * (1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']),
        'avf_wt': max(1, min(i[-1].lp.degree, CAPs[i[0]] * WTs[i[0]]['ag_wt'] * WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES
]


@pytest.mark.parametrize("secpar,sp,expected_pp", MAKE_SETUP_PARAMETERS_CASES)
def test_make_setup_parameters(mocker, secpar, sp, expected_pp):
    mocker.patch('lattice_cryptography.one_time_keys.random_polynomialvector', return_value=sp.key_ch)
    assert expected_pp == make_setup_parameters(secpar)


def test_correct_main():
    for j in ALLOWABLE_SECPARS:
        pp = make_setup_parameters(secpar=j)
        for i in range(SAMPLE_SIZE):
            some_keys = keygen(pp=pp, num=pp['ag_cap'])
            some_otvks = [k[2] for k in some_keys]
            some_msgs = [bin(randbits(j))[2:].zfill(j) for _ in range(pp['ag_cap'])]
            some_sigs = [sign(pp=pp, otk=otk, msg=msg) for otk, msg in zip(some_keys, some_msgs)]
            for next_otvk, next_msg, next_sig in zip(some_otvks, some_msgs, some_sigs):
                assert verify(pp=pp, otvks=[next_otvk], msgs=[next_msg], sig=next_sig)
            ag_sig = aggregate(pp=pp, otvks=some_otvks, msgs=some_msgs, sigs=some_sigs)
            assert verify(pp=pp, otvks=some_otvks, msgs=some_msgs, sig=ag_sig)

