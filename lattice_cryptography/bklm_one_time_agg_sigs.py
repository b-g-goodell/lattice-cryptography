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

# Allowable (security parameter, aggregation capacity) pairs
ALLOWABLE_PARAMETERS = [(2**j, 2**k) for j in [7, 8] for k in range(1, 12)]

# Moduli
MODULI: Dict[Tuple[int, int], int] = dict()
MODULI[(128, 2)] = 955649
MODULI[(128, 4)] = 955649
MODULI[(128, 8)] = 910849
MODULI[(128, 16)] = 1062913
MODULI[(128, 32)] = 2005249
MODULI[(128, 64)] = 3892481
MODULI[(128, 128)] = 7667713
MODULI[(128, 256)] = 15236609
MODULI[(128, 512)] = 30416641
MODULI[(128, 1024)] = 60639233
MODULI[(128, 2048)] = 122849537

MODULI[(256, 2)] = 7360513
MODULI[(256, 4)] = 7200769
MODULI[(256, 8)] = 7201793
MODULI[(256, 16)] = 8083457
MODULI[(256, 32)] = 15326209
MODULI[(256, 64)] = 29669377
MODULI[(256, 128)] = 58697729
MODULI[(256, 256)] = 115998721
MODULI[(256, 512)] = 236121089
MODULI[(256, 1024)] = 473172481
MODULI[(256, 2048)] = 1058768897

# Degrees
DEGREES: Dict[Tuple[int, int], int] = {(2**j, 2**k): 2**j for j in [7, 8] for k in range(1, 12)}

# Length/rank
LENGTHS: Dict[Tuple[int, int], int] = dict()
LENGTHS[(128, 2)] = 53
LENGTHS[(128, 4)] = 54
LENGTHS[(128, 8)] = 56
LENGTHS[(128, 16)] = 58
LENGTHS[(128, 32)] = 61
LENGTHS[(128, 64)] = 63
LENGTHS[(128, 128)] = 66
LENGTHS[(128, 256)] = 68
LENGTHS[(128, 512)] = 71
LENGTHS[(128, 1024)] = 74
LENGTHS[(128, 2048)] = 76

LENGTHS[(256, 2)] = 52
LENGTHS[(256, 4)] = 54
LENGTHS[(256, 8)] = 55
LENGTHS[(256, 16)] = 57
LENGTHS[(256, 32)] = 59
LENGTHS[(256, 64)] = 61
LENGTHS[(256, 128)] = 64
LENGTHS[(256, 256)] = 66
LENGTHS[(256, 512)] = 68
LENGTHS[(256, 1024)] = 70
LENGTHS[(256, 2048)] = 72

# Bounds
BDs: Dict[tuple[int, int], Dict[str, int]] = {x: {'ch_bd': 1, 'ag_bd': 1} for x in ALLOWABLE_PARAMETERS}
BDs[(128, 2)]['sk_bd'] = 22
BDs[(128, 4)]['sk_bd'] = 22
BDs[(128, 8)]['sk_bd'] = 21
BDs[(128, 16)]['sk_bd'] = 21
BDs[(128, 32)]['sk_bd'] = 21
BDs[(128, 64)]['sk_bd'] = 22
BDs[(128, 128)]['sk_bd'] = 21
BDs[(128, 256)]['sk_bd'] = 21
BDs[(128, 512)]['sk_bd'] = 21
BDs[(128, 1024)]['sk_bd'] = 21
BDs[(128, 2048)]['sk_bd'] = 21

BDs[(256, 2)]['sk_bd'] = 46
BDs[(256, 4)]['sk_bd'] = 45
BDs[(256, 8)]['sk_bd'] = 45
BDs[(256, 16)]['sk_bd'] = 44
BDs[(256, 32)]['sk_bd'] = 44
BDs[(256, 64)]['sk_bd'] = 44
BDs[(256, 128)]['sk_bd'] = 44
BDs[(256, 256)]['sk_bd'] = 44
BDs[(256, 512)]['sk_bd'] = 44
BDs[(256, 1024)]['sk_bd'] = 45
BDs[(256, 2048)]['sk_bd'] = 45

# Weights
# Note: cannot generally set these all to DEGREE!
WTs: Dict[tuple[int, int], Dict[str, int]] = {x: {'ch_wt': 26, 'ag_wt': 26} if x[0] == 128 else {'ch_wt': 50, 'ag_wt': 50} for x in ALLOWABLE_PARAMETERS}
WTs[(128, 2)]['sk_wt'] = 127
WTs[(128, 4)]['sk_wt'] = 122
WTs[(128, 8)]['sk_wt'] = 124
WTs[(128, 16)]['sk_wt'] = 123
WTs[(128, 32)]['sk_wt'] = 123
WTs[(128, 64)]['sk_wt'] = 123
WTs[(128, 128)]['sk_wt'] = 123
WTs[(128, 256)]['sk_wt'] = 123
WTs[(128, 512)]['sk_wt'] = 123
WTs[(128, 1024)]['sk_wt'] = 123
WTs[(128, 2048)]['sk_wt'] = 127

WTs[(256, 2)]['sk_wt'] = 254
WTs[(256, 4)]['sk_wt'] = 253
WTs[(256, 8)]['sk_wt'] = 251
WTs[(256, 16)]['sk_wt'] = 251
WTs[(256, 32)]['sk_wt'] = 252
WTs[(256, 64)]['sk_wt'] = 252
WTs[(256, 128)]['sk_wt'] = 251
WTs[(256, 256)]['sk_wt'] = 252
WTs[(256, 512)]['sk_wt'] = 253
WTs[(256, 1024)]['sk_wt'] = 255
WTs[(256, 2048)]['sk_wt'] = 251

# Salts
SALTs: Dict[tuple[int, int], Dict[str, str]] = {
    x: {'sk_salt': 'SK_SALT', 'ch_salt': 'CH_SALT', 'ag_salt': 'AG_SALT'} for x in ALLOWABLE_PARAMETERS
}

# Distribution
DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT

BAD_INPUT_ERR: str = 'Cannot make setup parameters without a security parameter in [128, 256] and carrying '
BAD_INPUT_ERR += 'capacity in [2, 4, 8, ..., 2048], but had (security parameter, carrying capacity) = '


def make_setup_parameters(secpar: SecPar, carrying_capacity: int) -> PubPar:
    # TODO: modify to test provable security conditions?
    x: tuple[int, int] = (secpar, carrying_capacity)
    if x not in ALLOWABLE_PARAMETERS:
        raise ValueError(BAD_INPUT_ERR + str(x))
    result: PubPar = dict()

    result['scheme_parameters']: SchemeParameters = SchemeParameters(
        secpar=secpar, lp=LatticeParameters(modulus=MODULI[x], degree=DEGREES[x], length=LENGTHS[x]),
        distribution=DISTRIBUTION
    )

    result['sk_salt']: Msg = SALTs[x]['sk_salt']
    result['sk_bd']: int = BDs[x]['sk_bd']
    result['sk_wt']: int = WTs[x]['sk_wt']

    result['ch_salt']: Msg = SALTs[x]['ch_salt']
    result['ch_bd']: int = BDs[x]['ch_bd']
    result['ch_wt']: int = WTs[x]['ch_wt']

    result['ag_salt']: Msg = SALTs[x]['ag_salt']
    result['ag_bd']: int = BDs[x]['ag_bd']
    result['ag_wt']: int = WTs[x]['ag_wt']
    result['ag_cap']: int = x[1]

    result['vf_wt']: int = max(1, min(result['scheme_parameters'].lp.degree, result['sk_wt'] * (1 + result['ch_wt'])))
    result['vf_bd']: int = min(result['scheme_parameters'].lp.degree, result['sk_wt'], result['ch_wt'])
    result['vf_bd'] *= result['ch_bd']
    result['vf_bd'] += 1
    result['vf_bd'] *= result['sk_bd']

    result['avf_wt']: int = result['ag_cap'] * result['ag_wt'] * result['vf_wt']
    result['avf_wt'] = min(result['scheme_parameters'].lp.degree, result['avf_wt'])
    result['avf_wt'] = max(1, result['avf_wt'])
    result['avf_bd']: int = min(result['scheme_parameters'].lp.degree, result['ag_wt'], result['vf_wt'])
    result['avf_bd'] *= result['ag_bd'] * result['vf_bd'] * result['ag_cap']

    return result


def keygen(pp: PubPar, num: int = 1, seeds: list[SS] | None = None, multiprocessing: bool | None = None) -> list[OTK]:
    return kg(pp=pp, num=num, seeds=seeds, multiprocessing=multiprocessing)


def sign(pp: PubPar, otk: OTK, msg: Msg) -> Sig:
    return sg(pp=pp, otk=otk, msg=msg)


BAD_PREP_INPUT_ERR: str = 'Input list of one-time verification keys must be List[OTVK], but had (type(otvks), '
BAD_PREP_INPUT_ERR += 'type(otvk for otvk in otvks)) = '


def prepare_make_agg_coefs(otvks: List[OTVK], msgs: List[Msg]) -> tuple:
    if len(otvks) != len(msgs):
        raise ValueError("Cannot prepare_make_agg_coefs without two input vectors of equal length.")
    elif not all(is_bitstring(msg) for msg in msgs):
        raise ValueError("Input messages must be bitstrings.")
    elif any(not isinstance(otvk, OTVK) for otvk in otvks):
        raise ValueError(BAD_PREP_INPUT_ERR + str((type(otvks), [type(otvk) for otvk in otvks])))
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


def prep_agg_input(otvks: List[OTVK], msgs: List[Msg], sigs: List[Sig]) -> tuple:
    zipped_data: List[Tuple[OTVK, Msg, Sig]] = list(zip(otvks, msgs, sigs))
    zipped_data_srt_by_key: List[Tuple[OTVK, Msg, Sig]] = sorted(zipped_data, key=lambda x: str(x[0]))
    return tuple(([i[j] for i in zipped_data_srt_by_key] for j in range(len(zipped_data_srt_by_key[0]))))


def aggregate(pp: PubPar, otvks: List[OTVK], msgs: List[Msg], sigs: List[Sig]) -> Sig:
    srt_keys_msgs_sigs: Tuple[List[OTVK], List[Msg], List[Sig]] = prep_agg_input(otvks=otvks, msgs=msgs, sigs=sigs)
    srt_keys: List[OTVK]
    srt_msgs: List[Msg]
    srt_sigs: List[Sig]
    srt_keys, srt_msgs, srt_sigs = srt_keys_msgs_sigs
    for srt_sig in srt_sigs:
        srt_sig.const_time_flag = False
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=otvks, msgs=msgs)
    return sum([sig ** ag_coef for sig, ag_coef in zip(srt_sigs, ag_coefs)])


def make_verify_target(pp: PubPar, otvk: OTVK, msg: Msg) -> Polynomial:
    c: Polynomial = make_challenge(pp=pp, otvk=otvk, msg=msg, const_time_flag=otvk.left_key.const_time_flag)
    return otvk.left_key * c + otvk.right_key


def make_aggregate_verify_target(pp: PubPar, otvks: List[OTVK], msgs: List[Msg]) -> Polynomial:
    if not all(otvk[0].const_time_flag == otvk[1].const_time_flag == otvks[0][0].const_time_flag for otvk in otvks):
        raise ValueError('Cannot verify without equal const-time-flags in all verification keys.')
    challs: List[Polynomial] = [make_challenge(pp=pp, otvk=otvk, msg=msg) for otvk, msg in zip(otvks, msgs)]
    zip_keys_msgs_chs: List[Tuple[OTVK, Msg, Chall]] = list(zip(otvks, msgs, challs))
    srtd_keys_and_msgs_and_challs: List[Tuple[OTVK, Msg, Chall]] = sorted(zip_keys_msgs_chs, key=lambda x: str(x[0]))
    srtd_keys: List[OTVK]
    srtd_msgs: List[Msg]
    srtd_challs: List[Chall]
    srtd_keys, srtd_msgs, srtd_challs = [[i[j] for i in srtd_keys_and_msgs_and_challs] for j in range(3)]
    ag_coefs: List[AggCoef] = make_agg_coefs(pp=pp, otvks=srtd_keys, msgs=srtd_msgs)
    return sum([(otvk[0] * c + otvk[1]) * ag_coef for ag_coef, c, otvk in zip(ag_coefs, srtd_challs, srtd_keys)])


def verify(pp: PubPar, otvks: List[OTVK], msgs: List[Msg], sig: Sig) -> bool:
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
