from lattice_algebra import \
    LatticeParameters, \
    Polynomial, \
    PolynomialVector, \
    is_bitstring, \
    bits_to_decode as bits_per_coefficient, \
    bits_to_indices as bits_per_index_set, \
    hash2polynomial
from lattice_cryptography.one_time_keys import \
    UNIFORM_INFINITY_WEIGHT, \
    SecPar, \
    PubPar, \
    SecretSeed, \
    OneTimeSigningKey, \
    OneTimeVerificationKey, \
    OneTimeKeyTuple, \
    Msg, \
    Chall, \
    Sig, \
    SchemeParameters, \
    keygen_core as kg, \
    challenge_core as make_challenge, \
    sign_core as sg, \
    verify_core
from typing import Dict, List, Tuple, Any

# Typing
SS = SecretSeed
OTSK = OneTimeSigningKey
OTVK = OneTimeVerificationKey
OTK = OneTimeKeyTuple
AggCoef = Polynomial
AggSig = Sig

# COMPARE THE PARAMETERS HERE WITH OUR PARAMETER ANALYSIS
LPs: Dict[int, LatticeParameters] = dict()
LPs[128]: LatticeParameters = LatticeParameters(
    modulus=78593,
    degree=2 ** 6,
    length=143)
LPs[256]: LatticeParameters = LatticeParameters(
    modulus=355841,
    degree=2 ** 8,
    length=90)

BDs: Dict[int, Dict[str, int]] = dict()
BDs[128]: Dict[str, int] = {
    'sk_bd': 76,
    'ch_bd': 2,
    'ag_bd': 1}
BDs[256]: Dict[str, int] = {
    'sk_bd': 172,
    'ch_bd': 1,
    'ag_bd': 1}

WTs: Dict[int, Dict[str, int]] = dict()
WTs[128]: Dict[str, int] = {
    'sk_wt': LPs[128].degree,
    'ch_wt': LPs[128].degree,
    'ag_wt': LPs[128].degree}
WTs[256]: Dict[str, int] = {
    'sk_wt': LPs[256].degree,
    'ch_wt': LPs[256].degree,
    'ag_wt': LPs[256].degree}

SALTs: Dict[int, Dict[str, str]] = dict()
SALTs[128]: Dict[str, str] = {
    'sk_salt': 'SK_SALT',
    'ch_salt': 'CH_SALT',
    'ag_salt': 'AG_SALT'}
SALTs[256]: Dict[str, str] = {
    'sk_salt': 'SK_SALT',
    'ch_salt': 'CH_SALT',
    'ag_salt': 'AG_SALT'}

CAPs: Dict[int, int] = dict()
CAPs[128] = 2
CAPs[256] = 2

DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT


def make_setup_parameters(secpar: SecPar) -> PubPar:
    result: PubPar = dict()

    result['scheme_parameters']: SchemeParameters = SchemeParameters(
        secpar=secpar, lp=LPs[secpar], distribution=DISTRIBUTION
    )

    result['sk_salt']: Msg = SALTs[secpar]['sk_salt']
    result['sk_bd']: int = BDs[secpar]['sk_bd']
    result['sk_wt']: int = WTs[secpar]['sk_wt']

    result['ch_salt']: Msg = SALTs[secpar]['ch_salt']
    result['ch_bd']: int = BDs[secpar]['ch_bd']
    result['ch_wt']: int = WTs[secpar]['ch_wt']

    result['ag_salt']: Msg = SALTs[secpar]['ag_salt']
    result['ag_bd']: int = BDs[secpar]['ag_bd']
    result['ag_wt']: int = WTs[secpar]['ag_wt']
    result['ag_cap']: int = CAPs[secpar]

    result['vf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['sk_wt'] * (1 + result['ch_wt'])))
    result['vf_bd']: int = result['sk_bd'] * (
            1 + min(result['scheme_parameters'].lp.degree, result['sk_wt'], result['ch_wt']) * result['ch_bd'])

    result['avf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree,
                                       result['ag_cap'] * result['ag_wt'] * result['vf_wt']))
    result['avf_bd']: int = result['ag_cap'] * min(result['scheme_parameters'].lp.degree, result['ag_wt'],
                                                   result['vf_wt']) * result['ag_bd'] * result['vf_bd']

    return result


def keygen(pp: PubPar, num: int = 1, seeds: list[SS] | None = None, multiprocessing: bool | None = None) -> list[OTK]:
    return kg(pp=pp, num=num, seeds=seeds, multiprocessing=multiprocessing)


def sign(pp: PubPar, otk: OTK, msg: Msg) -> Sig:
    return sg(pp=pp, otk=otk, msg=msg)


def prepare_make_agg_coefs(otvks: List[OTVK], msgs: List[Msg]) -> tuple:
    if len(otvks) != len(msgs):
        raise ValueError("Cannot prepare_make_agg_coefs without two input vectors of equal length.")
    elif not all(is_bitstring(msg) for msg in msgs):
        raise ValueError("Input messages must be bitstrings.")
    elif any(not isinstance(otvk, OTVK) for otvk in otvks):
        raise ValueError(
            f"Input list of one-time verification keys must be List[OTVK], but had " +
            f"type(otvks) = {type(otvks)} and type(otvk for otvk in otvks) = {[type(otvk) for otvk in otvks]}")
    zipped_keys_and_msgs: List[Tuple[OTVK, Msg]] = list(zip(otvks, msgs))
    zipped_srtd_keys_and_msgs: List[Tuple[OTVK, Msg]] = sorted(zipped_keys_and_msgs, key=lambda x: str(x[0]))
    return tuple(([i[j] for i in zipped_srtd_keys_and_msgs] for j in range(len(zipped_srtd_keys_and_msgs[0]))))


def prepare_hash2polyinput(pp: PubPar, otvks: List[OTVK], msgs: List[Msg]) -> dict[str: Any]:
    srt_keys_and_msgs: Tuple[List[OTVK], List[Msg]] = prepare_make_agg_coefs(otvks=otvks, msgs=msgs)
    srt_keys: List[OTVK]
    srt_msgs: List[Msg]
    srt_keys, srt_msgs = srt_keys_and_msgs
    bpc: int = bits_per_coefficient(secpar=pp['scheme_parameters'].secpar, bd=pp['ag_bd'])
    bpis: int = bits_per_index_set(
        secpar=pp['scheme_parameters'].secpar, degree=pp['scheme_parameters'].lp.degree, wt=pp['ag_wt'])
    msg: Msg = str(list(zip(srt_keys, srt_msgs)))
    return {'secpar': pp['scheme_parameters'].secpar, 'lp': pp['scheme_parameters'].lp,
            'distribution': pp['scheme_parameters'].dstn, 'dist_pars': {'bd': pp['ag_bd'], 'wt': pp['ag_wt']},
            'num_coefs': pp['ag_wt'], 'bti': bpis, 'btd': bpc, 'msg': msg, 'const_time_flag': False}


def make_agg_coefs(pp: PubPar, otvks: List[OTVK], msgs: List[Msg]) -> List[AggCoef]:
    hash2polyinput: dict[str, Any] = prepare_hash2polyinput(pp=pp, otvks=otvks, msgs=msgs)
    return [hash2polynomial(**hash2polyinput, salt=pp['ag_salt'] + str(i)) for i in range(len(otvks))]


def prepare_aggregate_input(otvks: List[OTVK], msgs: List[Msg], sigs: List[Sig]) -> tuple:
    zipped_data: List[Tuple[OTVK, Msg, Sig]] = list(zip(otvks, msgs, sigs))
    zipped_data_srt_by_key: List[Tuple[OTVK, Msg, Sig]] = sorted(zipped_data, key=lambda x: str(x[0]))
    return tuple(([i[j] for i in zipped_data_srt_by_key] for j in range(len(zipped_data_srt_by_key[0]))))


def aggregate(pp: PubPar, otvks: List[OTVK], msgs: List[Msg], sigs: List[Sig]) -> Sig:
    srt_keys_and_msgs_and_sigs: Tuple[List[OTVK], List[Msg], List[Sig]] = prepare_aggregate_input(otvks=otvks,
                                                                                                  msgs=msgs, sigs=sigs)
    srt_keys: List[OTVK]
    srt_msgs: List[Msg]
    srt_sigs: List[Sig]
    srt_keys, srt_msgs, srt_sigs = srt_keys_and_msgs_and_sigs
    for srt_sig in srt_sigs:
        srt_sig.const_time_flag = False
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=otvks, msgs=msgs)
    return sum([sig ** ag_coef for sig, ag_coef in zip(srt_sigs, ag_coefs)])


def make_verify_target(pp: PubPar, otvk: OTVK, msg: Msg) -> Polynomial:
    c: Polynomial = make_challenge(pp=pp, otvk=otvk, msg=msg, const_time_flag=otvk.left_key.const_time_flag)
    return otvk.left_key * c + otvk.right_key


def make_aggregate_verify_target(pp: PubPar, otvks: List[OTVK], msgs: List[Msg]) -> Polynomial:
    if not all(otvk.left_key.const_time_flag == otvk.right_key.const_time_flag == otvks[0].left_key.const_time_flag for
               otvk in otvks):
        raise ValueError('Cannot verify without equal const-time-flags in all verification keys.')
    challs: List[Polynomial] = [make_challenge(pp=pp, otvk=otvk, msg=msg) for otvk, msg in zip(otvks, msgs)]
    zipped_keys_and_msgs_and_challs: List[Tuple[OTVK, Msg, Chall]] = list(
        zip(otvks, msgs, challs))
    srtd_keys_and_msgs_and_challs: List[Tuple[OTVK, Msg, Chall]] = sorted(
        zipped_keys_and_msgs_and_challs, key=lambda x: str(x[0]))
    srtd_keys: List[OTVK]
    srtd_msgs: List[Msg]
    srtd_challs: List[Chall]
    srtd_keys, srtd_msgs, srtd_challs = [[i[j] for i in srtd_keys_and_msgs_and_challs] for j in range(3)]
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=srtd_keys, msgs=srtd_msgs)
    return sum([(otvk[0] * c + otvk[1]) * ag_coef for ag_coef, c, otvk in zip(ag_coefs, srtd_challs, srtd_keys)])


def verify(pp: PubPar, otvks: List[OTVK], msgs: List[Msg],
           sig: Sig) -> bool:
    if len(otvks) < 1:
        raise ValueError('Cannot verify without at least one key.')
    elif len(otvks) > pp['ag_cap']:
        raise ValueError('Cannot verify with more keys than capacity.')
    elif len(otvks) != len(msgs):
        raise ValueError('Cannot verify without same number of keys as messages.')
    elif any(otvk.left_key.const_time_flag != otvk.right_key.const_time_flag for otvk in otvks) or any(
            otvk.left_key.const_time_flag != otvks[0].left_key.const_time_flag for otvk in otvks):
        raise ValueError('Cannot verify without equal const-time-flags in the verification keys.')
    bd: int
    wt: int
    key_ch: PolynomialVector = pp['scheme_parameters'].key_ch
    if len(otvks) == 1:
        bd: int = pp['vf_bd']
        wt: int = pp['vf_wt']
        t: Polynomial = make_verify_target(pp=pp, otvk=otvks[0], msg=msgs[0])
    else:
        bd: int = pp['avf_bd']
        wt: int = pp['avf_wt']
        t: Polynomial = make_aggregate_verify_target(pp=pp, otvks=otvks, msgs=msgs)
    return verify_core(sig=sig, bd=bd, wt=wt, key_ch=key_ch, target=t)
