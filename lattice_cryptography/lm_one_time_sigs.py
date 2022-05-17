from lattice_algebra import \
    LatticeParameters, \
    Polynomial
from lattice_cryptography.one_time_keys import \
    SecPar, \
    PubPar, \
    SecretSeed, \
    OneTimeSigningKey, \
    OneTimeVerificationKey, \
    Msg, \
    Chall, \
    Sig, \
    ALLOWABLE_SECPARS, \
    SchemeParameters as SchPar, \
    UNIFORM_INFINITY_WEIGHT, \
    keygen_core as kg, \
    challenge_core as make_challenge, \
    sign_core as sg, \
    verify_core
from typing import Dict

# Typing
SS = SecretSeed
OTSK = OneTimeSigningKey
OTVK = OneTimeVerificationKey
OTK = tuple[SS, OTSK, OTVK]

# COMPARE THE PARAMETERS HERE WITH OUR PARAMETER ANALYSIS
LPs: Dict[int, LatticeParameters] = dict()
LPs[128] = LatticeParameters(
    modulus=78593,
    degree=2 ** 6,
    length=143)
LPs[256] = LatticeParameters(
    modulus=314113,
    degree=2 ** 7,
    length=156)

BDs: Dict[int, Dict[str, int]] = dict()
BDs[128] = {
    'sk_bd': 76,
    'ch_bd': 2}
BDs[256] = {
    'sk_bd': 152,
    'ch_bd': 2}

WTs: Dict[int, Dict[str, int]] = dict()
WTs[128] = {'sk_wt': LPs[128].degree, 'ch_wt': LPs[128].degree}
WTs[256] = {'sk_wt': LPs[256].degree, 'ch_wt': LPs[256].degree}

SALTs: Dict[int, Dict[str, str]] = {i: {'sk_salt': 'SK_SALT', 'ch_salt': 'CH_SALT'} for i in ALLOWABLE_SECPARS}

DSTN: str = UNIFORM_INFINITY_WEIGHT


def make_setup_parameters(secpar: SecPar) -> PubPar:
    result: PubPar = dict()
    result['scheme_parameters']: SchPar = SchPar(secpar=secpar, lp=LPs[secpar], distribution=DSTN)
    result['sk_salt']: Msg = SALTs[secpar]['sk_salt']
    result['sk_bd']: int = BDs[secpar]['sk_bd']
    result['sk_wt']: int = WTs[secpar]['sk_wt']
    result['ch_salt']: Msg = SALTs[secpar]['ch_salt']
    result['ch_bd']: int = BDs[secpar]['ch_bd']
    result['ch_wt']: int = WTs[secpar]['ch_wt']
    result['vf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['sk_wt'] * (1 + result['ch_wt'])))
    result['vf_bd']: int = result['sk_bd'] * (
                1 + min(result['scheme_parameters'].lp.degree, result['sk_wt'], result['ch_wt']) * result['ch_bd'])
    return result


def keygen(pp: PubPar, num: int = 1, seeds: list[SS] | None = None, multiprocessing: bool | None = None) -> list[OTK]:
    return kg(pp=pp, num=num, seeds=seeds, multiprocessing=multiprocessing)


def sign(pp: PubPar, otk: OTK, msg: Msg) -> Sig:
    return sg(pp=pp, otk=otk, msg=msg)


def make_verify_target(pp: PubPar, otvk: OTVK, msg: Msg) -> Polynomial:
    c: Chall = make_challenge(pp=pp, otvk=otvk, msg=msg, const_time_flag=otvk.left_key.const_time_flag)
    return otvk.left_key * c + otvk.right_key


def verify(pp: PubPar, otvk: OTVK, msg: Msg, sig: Sig) -> bool:
    if otvk.left_key.const_time_flag != otvk.right_key.const_time_flag:
        raise ValueError('Cannot verify without equal const-time-flags in the verification keys.')
    return verify_core(sig=sig, bd=pp['vf_bd'], wt=pp['vf_wt'], key_ch=pp['scheme_parameters'].key_ch,
                       target=make_verify_target(pp=pp, otvk=otvk, msg=msg))
