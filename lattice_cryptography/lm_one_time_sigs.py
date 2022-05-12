from lattice_algebra import LatticeParameters
from lattice_cryptography.one_time_keys import SecurityParameter, PublicParameters, SecretSeed, OneTimeSigningKey, OneTimeVerificationKey, OneTimeKeyTuple, Message, Challenge, Signature, ALLOWABLE_SECPARS, SchemeParameters, UNIFORM_INFINITY_WEIGHT, make_random_seed, make_one_key, keygen_core as keygen, sign_core as sign, verify_core as verify
from typing import Dict

# COMPARE THE PARAMETERS HERE WITH OUR PARAMETER ANALYSIS
LPs: Dict[int, LatticeParameters] = dict()
LPs[128] = LatticeParameters(modulus=78593, degree=2**6, length=143)
LPs[256] = LatticeParameters(modulus=355841, degree=2**8, length=90)

BDs: Dict[int, Dict[str, int]] = dict()
BDs[128] = {'sk_bd': 76, 'ch_bd': 2}
BDs[256] = {'sk_bd': 172, 'ch_bd': 1}

WTs: Dict[int, Dict[str, int]] = dict()
WTs[128] = {'sk_wt': LPs[128].degree, 'ch_wt': LPs[128].degree}
WTs[256] = {'sk_wt': LPs[256].degree, 'ch_wt': LPs[256].degree}

SALTs: Dict[int, Dict[str, str]] = {i: {'sk_salt': 'SK_SALT', 'ch_salt': 'CH_SALT'} for i in ALLOWABLE_SECPARS}

DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT


def make_setup_parameters(secpar: SecurityParameter) -> PublicParameters:
    result: PublicParameters = {}

    result['scheme_parameters']: SchemeParameters = SchemeParameters(
        secpar=secpar, lp=LPs[secpar],
        distribution=DISTRIBUTION
    )

    result['sk_salt']: Message = SALTs[secpar]['sk_salt']
    result['sk_bd']: int = BDs[secpar]['sk_bd']
    result['sk_wt']: int = WTs[secpar]['sk_wt']

    result['ch_salt']: Message = SALTs[secpar]['ch_salt']
    result['ch_bd']: int = BDs[secpar]['ch_bd']
    result['ch_wt']: int = WTs[secpar]['ch_wt']

    result['vf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['sk_wt'] * (1 + result['ch_wt'])))
    result['vf_bd']: int = result['sk_bd'] * (1 + min(result['scheme_parameters'].lp.degree, result['sk_wt'], result['ch_wt']) * result['ch_bd'])
    return result
