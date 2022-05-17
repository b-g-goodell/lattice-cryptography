"""
We test the lattice-crypto.keys.one_time_keys module.
TODO: Fix this so it does not depend on lattice parameters and salts and stuff from lm_one_time_sigs
"""
import pytest
from lattice_cryptography.one_time_keys import *
from lattice_cryptography.lm_one_time_sigs import LPs, SALTs, BDs, WTs
from secrets import randbits

SAMPLE_SIZE: int = 2 ** 3
DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT
MAKE_SCHEME_PARAMETERS_CASES = [(
    i,
    SchemeParameters(secpar=i, lp=LPs[i], distribution=DISTRIBUTION)
) for i in ALLOWABLE_SECPARS]
MAKE_SETUP_PARAMETERS_CASES = [
    i + tuple([{
        'scheme_parameters': i[-1],
        'sk_salt': SALTs[i[0]]['sk_salt'],
        'sk_bd': BDs[i[0]]['sk_bd'],
        'sk_wt': WTs[i[0]]['sk_wt'],
        'ch_salt': SALTs[i[0]]['ch_salt'],
        'ch_bd': BDs[i[0]]['ch_bd'],
        'ch_wt': WTs[i[0]]['ch_wt'],
        'vf_bd': BDs[i[0]]['sk_bd'] * (1 + min(WTs[i[0]]['sk_wt'], WTs[i[0]]['ch_wt']) * BDs[i[0]]['ch_bd']),
        'vf_wt': max(1, min(i[-1].lp.degree, WTs[i[0]]['sk_wt'] * (1 + WTs[i[0]]['ch_wt']))),
    }]) for i in MAKE_SCHEME_PARAMETERS_CASES
]


@pytest.mark.parametrize("secpar,sp,pp", MAKE_SETUP_PARAMETERS_CASES)
def test_correct_main(secpar, sp, pp):
    for _ in range(SAMPLE_SIZE):
        otk: OTK = keygen_core(pp=pp, num=1)[0]
        otvk: OTVK = otk[2]
        msg: Msg = bin(randbits(secpar))[2:].zfill(secpar)
        sig: Sig = sign_core(pp=pp, otk=otk, msg=msg)
        t: Polynomial = otvk[0] * challenge_core(pp=pp, otvk=otvk, msg=msg) + otvk[1]
        assert verify_core(sig=sig, bd=pp['vf_bd'], wt=pp['vf_wt'], key_ch=pp['scheme_parameters'].key_ch, target=t)
