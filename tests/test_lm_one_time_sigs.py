"""
We test the lattice-crypto.lm_one_time_sigs module.
"""
import pytest
from lattice_cryptography.one_time_keys import \
    ALLOWABLE_SECPARS, \
    SecretSeed, \
    OneTimeSigningKey, \
    OneTimeVerificationKey, \
    SchemeParameters as SchPar
from lattice_cryptography.lm_one_time_sigs import \
    LPs, \
    SALTs, \
    BDs, \
    WTs, \
    DSTN, \
    Msg, \
    Sig, \
    OTK, \
    make_setup_parameters, \
    keygen, \
    sign, \
    verify
from secrets import randbits

SAMPLE_SIZE: int = 2 ** 3
SS = SecretSeed
OTSK = OneTimeSigningKey
OTVK = OneTimeVerificationKey
OTK = tuple[SS, OTSK, OTVK]
REQUIRED_SETUP_PARAMETERS = ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt', 'ch_bd', 'ch_wt', 'ch_salt', 'vf_bd', 'vf_wt']
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
        'vf_bd': BDs[i[0]]['sk_bd'] * (
                    1 + min(i[-1].lp.degree, WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']),
        'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
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
            otk: OTK = keygen(pp=pp, num=1)[0]
            otvk: OTVK = otk[2]
            msg: Msg = bin(randbits(j))[2:].zfill(j)
            sig: Sig = sign(pp=pp, otk=otk, msg=msg)
            assert verify(pp=pp, otvk=otvk, msg=msg, sig=sig)
