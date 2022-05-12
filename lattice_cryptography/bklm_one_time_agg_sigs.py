from lattice_algebra import Polynomial, PolynomialVector, LatticeParameters, hash2polynomialvector, hash2polynomial, is_bitstring
from lattice_cryptography.one_time_keys import SecretSeed, OneTimeSigningKey, OneTimeVerificationKey, ALLOWABLE_SECPARS, SchemeParameters, UNIFORM_INFINITY_WEIGHT, bits_per_index_set, bits_per_coefficient
from typing import Dict, List, Tuple, Any
from secrets import randbelow
from multiprocessing import Pool, cpu_count
from math import ceil

# Typing
SecurityParameter = int
PublicParameters = Dict[str, Any]
OneTimeKeyTuple = Tuple[SecretSeed, OneTimeSigningKey, OneTimeVerificationKey]
Message = str
Challenge = Polynomial
Signature = PolynomialVector
AggCoef = Polynomial
AggSig = Signature

# COMPARE THE PARAMETERS HERE WITH OUR PARAMETER ANALYSIS
LPs: Dict[int, LatticeParameters] = dict()
LPs[128]: LatticeParameters = LatticeParameters(modulus=78593, degree=2**6, length=143)
LPs[256]: LatticeParameters = LatticeParameters(modulus=355841, degree=2**8, length=90)

BDs: Dict[int, Dict[str, int]] = dict()
BDs[128]: Dict[str, int] = {'sk_bd': 76, 'ch_bd': 2, 'ag_bd': 1}
BDs[256]: Dict[str, int] = {'sk_bd': 172, 'ch_bd': 1, 'ag_bd': 1}

WTs: Dict[int, Dict[str, int]] = dict()
WTs[128]: Dict[str, int] = {'sk_wt': LPs[128].degree, 'ch_wt': LPs[128].degree, 'ag_wt': LPs[128].degree}
WTs[256]: Dict[str, int]= {'sk_wt': LPs[256].degree, 'ch_wt': LPs[256].degree, 'ag_wt': LPs[256].degree}

SALTs: Dict[int, Dict[str, str]] = dict()
SALTs[128]: Dict[str, str] = {'sk_salt': 'SK_SALT', 'ch_salt': 'CH_SALT', 'ag_salt': 'AG_SALT'}
SALTs[256]: Dict[str, str] = {'sk_salt': 'SK_SALT', 'ch_salt': 'CH_SALT', 'ag_salt': 'AG_SALT'}

DISTRIBUTION: str = UNIFORM_INFINITY_WEIGHT


def make_setup_parameters(secpar: SecurityParameter) -> PublicParameters:
    result: PublicParameters = {}

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


def make_random_seed(secpar: SecurityParameter, pp: PublicParameters) -> SecretSeed:
    # TODO: Move to one_time_keys.py
    seed: Message = bin(randbelow(2 ** secpar))[2:].zfill(secpar)
    return SecretSeed(secpar=secpar, lp=pp['scheme_parameters'].lp, seed=seed)


def make_one_key(pp: PublicParameters, seed: SecretSeed = None) -> OneTimeKeyTuple:
    # TODO: Move to one_time_keys.py
    secpar: SecurityParameter = pp['scheme_parameters'].secpar
    lp: LatticeParameters = pp['scheme_parameters'].lp
    x: SecretSeed = seed
    if not x:
        x = make_random_seed(secpar=secpar, pp=pp)
    left_signing_key: PolynomialVector = hash2polynomialvector(
        secpar=secpar, lp=lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': pp['sk_bd'], 'wt': pp['sk_wt']},
        num_coefs=pp['sk_wt'],
        bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=pp['sk_wt']),
        btd=bits_per_coefficient(secpar=secpar, bd=pp['sk_bd']),
        salt=pp['sk_salt'] + 'LEFT',
        msg=x.seed,
        const_time_flag=True
    )
    right_signing_key: PolynomialVector = hash2polynomialvector(
        secpar=secpar, lp=lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': pp['sk_bd'], 'wt': pp['sk_wt']},
        num_coefs=pp['sk_wt'],
        bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=pp['sk_wt']),
        btd=bits_per_coefficient(secpar=secpar, bd=pp['sk_bd']),
        salt=pp['sk_salt'] + 'RIGHT',
        msg=x.seed,
        const_time_flag=True
    )
    otsk: OneTimeSigningKey = OneTimeSigningKey(secpar=secpar, lp=lp, left_key=left_signing_key, right_key=right_signing_key)
    key_ch: PolynomialVector = pp['scheme_parameters'].key_ch
    key_ch.const_time_flag = True
    otvk: OneTimeVerificationKey = OneTimeVerificationKey(secpar=secpar, lp=lp, left_key=key_ch * left_signing_key, right_key=key_ch * right_signing_key)
    return x, otsk, otvk


def keygen(pp: PublicParameters, num_keys_to_gen: int = 1, seeds: List[SecretSeed] = None,
           multiprocessing: bool = None) -> List[OneTimeKeyTuple]:
    # TODO: Move to one_time_keys.py, import from one_time_keys, refactor so this is a wrapper
    """ Wraps keygen_core to handle workload distribution for batch generation """

    # Only default to parallelization if more than a few keys are needed (to avoid unnecessary overhead)
    if multiprocessing is None:
        multiprocessing: bool = num_keys_to_gen >= 16
    num_workers: int = cpu_count()

    # Pass straight through to keygen_core() if there is no reason or desire to parallelize (to avoid extra overhead)
    if (not multiprocessing) or (num_keys_to_gen == 1) or (num_workers == 1):
        return keygen_core(pp=pp, num_keys=num_keys_to_gen, seeds=seeds)

    # Prepare inputs for the pool
    if not seeds:
        iterable: List[Tuple[Dict[str, Any], int]] = [(pp, ceil(num_keys_to_gen / num_workers))] * num_workers
    else:
        seed_batches: List[List[Any]] = distribute_tasks(tasks=seeds)
        iterable: List[tuple] = list(zip([pp] * len(seed_batches), [len(x) for x in seed_batches], seed_batches))

    # Generate the keys and return the flattened results
    with Pool(num_workers) as pool:
        nested_keys: List[List[OneTimeKeyTuple]] = pool.starmap(func=keygen, iterable=iterable)
    return [item for sublist in nested_keys for item in sublist][:num_keys_to_gen]


def keygen_core(pp: PublicParameters, num_keys: int = 1, seeds: List[SecretSeed] = None) -> List[OneTimeKeyTuple]:
    # TODO: Move to one_time_keys.py
    if num_keys < 1:
        raise ValueError('Can only generate a natural number worth of keys.')
    elif seeds is not None and len(seeds) != num_keys:
        raise ValueError('Must either roll keys with no seeds, or with a seed for each key.')
    elif seeds is None and num_keys == 1:
        return [make_one_key(pp=pp)]
    elif seeds is not None and num_keys == 1:
        return [make_one_key(pp=pp, seed=seeds[0])]
    elif seeds is None:
        return [make_one_key(pp=pp) for _ in range(num_keys)]
    return [make_one_key(pp=pp, seed=next_seed) for next_seed in seeds]


def make_signature_challenge(pp: PublicParameters, otvk: OneTimeVerificationKey, msg: Message) -> Challenge:
    # TODO: Rename challenge_core, move to one_time_keys.py
    return hash2polynomial(
        secpar=pp['scheme_parameters'].secpar,
        lp=pp['scheme_parameters'].lp,
        distribution=DISTRIBUTION,
        dist_pars={'bd': pp['ch_bd'], 'wt': pp['ch_wt']},
        salt=pp['ch_salt'],
        msg=str(otvk) + ', ' + msg,
        num_coefs=pp['ch_wt'],
        bti=bits_per_index_set(
            secpar=pp['scheme_parameters'].secpar,
            degree=pp['scheme_parameters'].lp.degree,
            wt=pp['ch_wt']
        ),
        btd=bits_per_coefficient(
            secpar=pp['scheme_parameters'].secpar,
            bd=pp['ch_bd']
        ),
        const_time_flag=True
    )


def sign(pp: PublicParameters, otk: OneTimeKeyTuple, msg: Message) -> Signature:
    # TODO: Rename sign_core, move to one_time_keys.py, import sign_core, refactor so this is a wrapper
    c: Challenge = make_signature_challenge(pp=pp, otvk=otk[2], msg=msg)
    signature: Signature = otk[1][0] ** c + otk[1][1]
    return signature


def verify(pp: PublicParameters, otvk: OneTimeVerificationKey, msg: Message, sig: Signature) -> bool:
    # TODO: Rename verify_core, move to one_time_keys.py, import sign_core, refactor so this is a wrapper
    sig.const_time_flag = False
    cnws: List[Tuple[Dict[int, int], int, int]] = sig.get_coef_rep()
    n: int = max(i[1] for i in cnws)
    w: int = max(i[2] for i in cnws)

    if n > pp['vf_bd'] or w > pp['vf_wt']:
        return False

    key_ch: PolynomialVector = pp['scheme_parameters'].key_ch
    c: Challenge = make_signature_challenge(pp=pp, otvk=otvk, msg=msg)

    key_ch.const_time_flag = False
    c.const_time_flag = False
    otvk.left_key.const_time_flag = False
    otvk.right_key.const_time_flag = False

    lhs: Polynomial = key_ch * sig
    rhs: Polynomial = otvk[0] * c + otvk[1]

    return lhs == rhs


def distribute_tasks(tasks: List[Any], num_workers: int = None) -> List[List[Any]]:
    """
    Helper function that distributes a list of arbitrary tasks among a specific number of workers

    :param tasks: iterable containing list of tasks to carry out
    :param num_workers: number of workers available in the pool (usually = number of CPU cores)
    :return: task list broken up into num_workers segments
    TODO: move to one_time_keys.py
    """
    if not num_workers:
        num_workers = cpu_count()

    # Determine how the jobs should be split up per core
    r: int = len(tasks) % num_workers  # number of leftover jobs once all complete batches are processed
    job_counts: List[int] = r * [1 + (len(tasks) // num_workers)] + (num_workers - r) * [len(tasks) // num_workers]

    # Distribute the tasks accordingly
    i: int = 0
    task_list_all: List[List[Any]] = []
    for load_amount in job_counts:
        task_list_all.append(tasks[i:i + load_amount])
        i += load_amount
    return task_list_all


def prepare_make_agg_coefs(otvks: List[OneTimeVerificationKey], msgs: List[Message]) -> Tuple[List[OneTimeVerificationKey], List[Message]]:
    # TODO: Refactor to verify types, and then modify test to take this into account
    if len(otvks) != len(msgs):
        raise ValueError("Cannot prepare_make_agg_coefs without two input vectors of equal length.")
    elif not all(is_bitstring(msg) for msg in msgs):
        raise ValueError("Input messages must be bitstrings.")
    elif any(not isinstance(otvk, OneTimeVerificationKey) for otvk in otvks):
        raise ValueError(f"Input list ove one-time verification keys must be List[OneTimeVerificationKey], but had type(otvks) = {type(otvks)} and type(otvk for otvk in otvks) = {[type(otvk) for otvk in otvks]}")
    zipped_keys_and_msgs: List[Tuple[OneTimeVerificationKey, Message]] = list(zip(otvks, msgs))
    zipped_srtd_keys_and_msgs: List[Tuple[OneTimeVerificationKey, Message]] = sorted(zipped_keys_and_msgs, key=lambda x: str(x[0]))
    srt_keys, srt_msgs = ([i[j] for i in zipped_srtd_keys_and_msgs] for j in range(len(zipped_srtd_keys_and_msgs[0])))
    return srt_keys, srt_msgs


def prepare_hash2polyinput(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message]) -> dict[str: Any]:
    srt_keys, srt_msgs = prepare_make_agg_coefs(otvks=otvks, msgs=msgs)
    bpc: int = bits_per_coefficient(secpar=pp['scheme_parameters'].secpar, bd=pp['ag_bd'])
    bpis: int = bits_per_index_set(secpar=pp['scheme_parameters'].secpar, degree=pp['scheme_parameters'].lp.degree,
                             wt=pp['ag_wt'])
    msg: Message = str(list(zip(srt_keys, srt_msgs)))
    return {'secpar': pp['scheme_parameters'].secpar, 'lp': pp['scheme_parameters'].lp,
            'distribution': pp['scheme_parameters'].distribution, 'dist_pars': {'bd': pp['ag_bd'], 'wt': pp['ag_wt']},
            'num_coefs': pp['ag_wt'], 'bti': bpis, 'btd': bpc, 'msg': msg, 'const_time_flag': False}


def make_agg_coefs(pp: PublicParameters, otvks: List[OneTimeVerificationKey],
                   msgs: List[Message]) -> List[AggCoef]:
    hash2polyinput: dict[str, Any] = prepare_hash2polyinput(pp=pp, otvks=otvks, msgs=msgs)
    return [hash2polynomial(**hash2polyinput, salt=pp['ag_salt'] + str(i)) for i in range(len(otvks))]


def prepare_aggregate(otvks: List[OneTimeVerificationKey], msgs: List[Message], sigs: List[Signature]) -> Tuple[
    List[OneTimeVerificationKey], List[Message], List[Signature]]:
    zipped_data = list(zip(otvks, msgs, sigs))
    zipped_data_srt_by_key = sorted(zipped_data, key=lambda x: str(x[0]))
    srt_keys, srt_msgs, srt_sigs = ([i[j] for i in zipped_data_srt_by_key] for j in range(len(zipped_data[0])))
    return srt_keys, srt_msgs, srt_sigs


def aggregate(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message], sigs: List[Signature]) -> Signature:
    srt_keys, srt_msgs, srt_sigs = prepare_aggregate(otvks=otvks, msgs=msgs, sigs=sigs)
    ag_coefs = make_agg_coefs(pp=pp, otvks=otvks, msgs=msgs)
    return sum([sig ** ag_coef for sig, ag_coef in zip(srt_sigs, ag_coefs)])


def aggregate_verify(pp: PublicParameters, otvks: List[OneTimeVerificationKey], msgs: List[Message],
                     ag_sig: Signature) -> bool:
    cnw: List[Tuple[Dict[int, int], int, int]] = ag_sig.get_coef_rep()
    n, w = max(i[1] for i in cnw), max(i[2] for i in cnw)
    if n < 1 or n > pp['avf_bd'] or w < 1 or w > pp['avf_wt'] or len(otvks) < 1 or \
            len(otvks) > pp['ag_cap'] or len(otvks) != len(msgs):
        return False

    challenges: List[Polynomial] = [make_signature_challenge(pp=pp, otvk=otvk, msg=msg) for otvk, msg in
                                    zip(otvks, msgs)]
    zipped_keys_msgs_and_challs = list(zip(otvks, msgs, challenges))
    srt_otvks, sorted_msgs, sorted_challs = (
        [i[j] for i in sorted(zipped_keys_msgs_and_challs, key=lambda x: str(x[0]))]
        for j in range(len(zipped_keys_msgs_and_challs[0])))
    ag_coefs = make_agg_coefs(pp=pp, otvks=srt_otvks, msgs=sorted_msgs)
    sum_of_otvks: Polynomial = sum(
        [(otvk[0] * c + otvk[1]) * ag_coef for ag_coef, c, otvk in zip(ag_coefs, sorted_challs, srt_otvks)])
    return pp['scheme_parameters'].key_ch * ag_sig == sum_of_otvks
