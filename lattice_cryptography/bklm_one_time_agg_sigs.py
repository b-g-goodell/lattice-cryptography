from lattice_algebra import LatticeParameters, Polynomial, is_bitstring, bits_to_decode as bits_per_coefficient, bits_to_indices as bits_per_index_set, hash2polynomial
from lattice_cryptography.one_time_keys import SecurityParameter, PublicParameters, SecretSeed, OneTimeSigningKey, OneTimeVerificationKey, OneTimeKeyTuple, Message, Challenge, Signature, ALLOWABLE_SECPARS, SchemeParameters, UNIFORM_INFINITY_WEIGHT, make_random_seed, make_one_key, keygen_core as keygen, challenge_core as make_challenge, sign_core as sign, verify_core as verify
from typing import Dict, List, Tuple, Any

# Typing
AggCoef = Polynomial
AggSig = Signature

# COMPARE THE PARAMETERS HERE WITH OUR PARAMETER ANALYSIS
LPs: Dict[int, LatticeParameters] = dict()
LPs[128]: LatticeParameters = LatticeParameters(
    modulus=78593,
    degree=2**6,
    length=143)
LPs[256]: LatticeParameters = LatticeParameters(
    modulus=355841,
    degree=2**8,
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

DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT


def make_setup_parameters(secpar: SecurityParameter) -> PublicParameters:
    result: PublicParameters = dict()

    result['scheme_parameters']: SchemeParameters = SchemeParameters(
        secpar=secpar, lp=LPs[secpar], distribution=DISTRIBUTION
    )

    result['sk_salt']: Message = SALTs[secpar]['sk_salt']
    result['sk_bd']: int = BDs[secpar]['sk_bd']
    result['sk_wt']: int = WTs[secpar]['sk_wt']

    result['ch_salt']: Message = SALTs[secpar]['ch_salt']
    result['ch_bd']: int = BDs[secpar]['ch_bd']
    result['ch_wt']: int = WTs[secpar]['ch_wt']

    result['ag_salt']: Message = SALTs[secpar]['ag_salt']
    result['ag_bd']: int = BDs[secpar]['ag_bd']
    result['ag_wt']: int = WTs[secpar]['ag_wt']

    result['vf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['sk_wt'] * (1 + result['ch_wt'])))
    result['vf_bd']: int = result['sk_bd'] * (1 + min(result['scheme_parameters'].lp.degree, result['sk_wt'], result['ch_wt']) * result['ch_bd'])

    result['avf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['ag_cap'] * result['ag_wt'] * result['vf_wt']))
    result['avf_bd']: int = result['ag_cap'] * min(result['scheme_parameters'].lp.degree, result['ag_wt'], result['vf_wt']) * result['ag_bd'] * result['vf_bd']

    return result


def prepare_make_agg_coefs(otvks: List[OneTimeVerificationKey], msgs: List[Message]) -> tuple:
    # TODO: Refactor to verify types, and then modify test to take this into account
    # TODO: improve typing
    if len(otvks) != len(msgs):
        raise ValueError("Cannot prepare_make_agg_coefs without two input vectors of equal length.")
    elif not all(is_bitstring(msg) for msg in msgs):
        raise ValueError("Input messages must be bitstrings.")
    elif any(not isinstance(otvk, OneTimeVerificationKey) for otvk in otvks):
        raise ValueError(f"Input list ove one-time verification keys must be List[OneTimeVerificationKey], but had type(otvks) = {type(otvks)} and type(otvk for otvk in otvks) = {[type(otvk) for otvk in otvks]}")
    zipped_keys_and_msgs: List[Tuple[OneTimeVerificationKey, Message]] = list(zip(otvks, msgs))
    zipped_srtd_keys_and_msgs: List[Tuple[OneTimeVerificationKey, Message]] = sorted(zipped_keys_and_msgs, key=lambda x: str(x[0]))
    return tuple(([i[j] for i in zipped_srtd_keys_and_msgs] for j in range(len(zipped_srtd_keys_and_msgs[0]))))


def prepare_hash2polyinput(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message]) -> dict[str: Any]:
    srt_keys_and_msgs: Tuple[List[OneTimeVerificationKey], List[Message]] = prepare_make_agg_coefs(otvks=otvks, msgs=msgs)
    srt_keys: List[OneTimeVerificationKey]
    srt_msgs: List[Message]
    srt_keys, srt_msgs = srt_keys_and_msgs
    bpc: int = bits_per_coefficient(secpar=pp['scheme_parameters'].secpar, bd=pp['ag_bd'])
    bpis: int = bits_per_index_set(secpar=pp['scheme_parameters'].secpar, degree=pp['scheme_parameters'].lp.degree, wt=pp['ag_wt'])
    msg: Message = str(list(zip(srt_keys, srt_msgs)))
    return {'secpar': pp['scheme_parameters'].secpar, 'lp': pp['scheme_parameters'].lp,
            'distribution': pp['scheme_parameters'].distribution, 'dist_pars': {'bd': pp['ag_bd'], 'wt': pp['ag_wt']},
            'num_coefs': pp['ag_wt'], 'bti': bpis, 'btd': bpc, 'msg': msg, 'const_time_flag': False}


def make_agg_coefs(pp: PublicParameters, otvks: List[OneTimeVerificationKey],
                   msgs: List[Message]) -> List[AggCoef]:
    hash2polyinput: dict[str, Any] = prepare_hash2polyinput(pp=pp, otvks=otvks, msgs=msgs)
    return [hash2polynomial(**hash2polyinput, salt=pp['ag_salt'] + str(i)) for i in range(len(otvks))]


def prepare_aggregate_input(otvks: List[OneTimeVerificationKey], msgs: List[Message], sigs: List[Signature]) -> tuple:
    zipped_data: List[Tuple[OneTimeVerificationKey, Message, Signature]] = list(zip(otvks, msgs, sigs))
    zipped_data_srt_by_key: List[Tuple[OneTimeVerificationKey, Message, Signature]]  = sorted(zipped_data, key=lambda x: str(x[0]))
    return tuple(([i[j] for i in zipped_data_srt_by_key] for j in range(len(zipped_data_srt_by_key[0]))))
    # srt_keys, srt_msgs, srt_sigs = ([i[j] for i in zipped_data_srt_by_key] for j in range(len(zipped_data[0])))
    # return srt_keys, srt_msgs, srt_sigs


def aggregate(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message], sigs: List[Signature]) -> Signature:
    srt_keys_and_msgs_and_sigs: Tuple[List[OneTimeVerificationKey], List[Message], List[Signature]] = prepare_aggregate_input(otvks=otvks, msgs=msgs, sigs=sigs)
    srt_keys: List[OneTimeVerificationKey]
    srt_msgs: List[Message]
    srt_sigs: List[Signature]
    srt_keys, srt_msgs, srt_sigs = srt_keys_and_msgs_and_sigs
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=otvks, msgs=msgs)
    return sum([sig ** ag_coef for sig, ag_coef in zip(srt_sigs, ag_coefs)])


def aggregate_verify(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message],
                     ag_sig: Signature) -> bool:
    if len(otvks) < 1 or len(otvks) > pp['ag_cap'] or len(otvks) != len(msgs):
        return False

    cnw: List[Tuple[Dict[int, int], int, int]] = ag_sig.get_coef_rep()
    n: int
    w: int
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)

    if n < 1 or n > pp['avf_bd'] or w < 1 or w > pp['avf_wt']:
        return False

    challs: List[Polynomial] = [make_challenge(pp=pp, otvk=otvk, msg=msg) for otvk, msg in zip(otvks, msgs)]
    zipped_keys_and_msgs_and_challs: List[Tuple[OneTimeVerificationKey, Message, Challenge]] = list(zip(otvks, msgs, challs))
    srtd_keys_and_msgs_and_challs: List[Tuple[OneTimeVerificationKey, Message, Challenge]] = sorted(zipped_keys_and_msgs_and_challs, key=lambda x:str(x[0]))
    srtd_keys: List[OneTimeVerificationKey]
    srtd_msgs: List[Message]
    srtd_challs: List[Challenge]
    srtd_keys, srtd_msgs, srtd_challs = srtd_keys_and_msgs_and_challs
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=srtd_keys, msgs=srtd_msgs)
    sum_of_otvks: Polynomial = sum(
        [(otvk[0] * c + otvk[1]) * ag_coef for ag_coef, c, otvk in zip(ag_coefs, srtd_challs, srtd_keys)])
    return pp['scheme_parameters'].key_ch * ag_sig == sum_of_otvks
