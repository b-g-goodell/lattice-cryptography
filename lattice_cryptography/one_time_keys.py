"""
The lattice_cryptography.keys module handles keys and their generation.
"""
from lattice_algebra import Polynomial, PolynomialVector, LatticeParameters, random_polynomialvector, is_bitstring, UNIFORM_INFINITY_WEIGHT, hash2polynomial, hash2polynomialvector
from math import ceil, log2
from typing import Dict, Any, Tuple, List
from multiprocessing import Pool, cpu_count
from secrets import randbelow

# Typing
SecurityParameter = int
PublicParameters = Dict[str, Any]
Message = str
Challenge = Polynomial
Signature = PolynomialVector

# Allowed security parameters
ALLOWABLE_SECPARS = [128, 256]

# Error messages
GENERIC_ERR: str = 'Something went wrong.'
MISSING_DATA_ERR: str = 'Missing some required data.'
INCORRECT_DATA_TYPE_ERR: str = 'Required input data not the correct type.'
DATA_MISMATCH_ERR: str = 'Input data did not match.'
SEED_INST_ERR_NEED_BITS: str = INCORRECT_DATA_TYPE_ERR + ' Input must be a binary string.'
INVALID_DATA_VALUES_ERR: str = 'Required input data does not have valid values.'


class SecretSeed(object):
    secpar: int
    lp: LatticeParameters
    seed: str

    def __init__(self, seed: str, secpar: int, lp: LatticeParameters):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be an integer in' +
                             f' {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not is_bitstring(seed):
            raise ValueError(SEED_INST_ERR_NEED_BITS)
        elif not isinstance(lp, LatticeParameters):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be' +
                             f' LatticeParameters object.')
        elif len(seed) < secpar:
            raise ValueError(INVALID_DATA_VALUES_ERR + ' Input secret seed must have enough bits.')
        self.secpar = secpar
        self.lp = lp
        self.seed = seed

    def __eq__(self, other) -> bool:
        secpar_equal: bool = self.secpar == other.secpar
        lp_equal: bool = self.lp == other.lp
        seeds_equal: bool = self.seed == other.seed
        return secpar_equal and lp_equal and seeds_equal

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.seed)

    def __repr__(self):
        return str((self.secpar, self.lp, self.seed))


SECWIT_INST_ERR_NEED_VEC: str = MISSING_DATA_ERR
SECWIT_INST_ERR_NEED_VEC += ' Must instantiate SecretWitness with either a PolynomialVector.'
SECWIT_INST_ERR_NEED_POLYVEC: str = INCORRECT_DATA_TYPE_ERR
SECWIT_INST_ERR_NEED_POLYVEC += ' When instantiating SecretWitness with a PolynomialVector, input data must be a '
SECWIT_INST_ERR_NEED_POLYVEC += 'PolynomialVector!'
SECWIT_INST_ERR_NEED_SEED = INCORRECT_DATA_TYPE_ERR
SECWIT_INST_ERR_NEED_SEED += ' When instantiating SecretWitness with a SecretSeed, input data must be a SecretSeed!'
SECWIT_INST_ERR_LP_MISMATCH = DATA_MISMATCH_ERR
SECWIT_INST_ERR_LP_MISMATCH += ' Input LatticeParameters object does not match the LatticeParameters for the input '
SECWIT_INST_ERR_LP_MISMATCH += 'PolynomialVector.'
SECWIT_INST_ERR_SECPAR_MISMATCH = DATA_MISMATCH_ERR
SECWIT_INST_ERR_SECPAR_MISMATCH += ' The input LatticeParameters object or the input secpar integer does not match '
SECWIT_INST_ERR_SECPAR_MISMATCH += 'the input SecretSeed.'
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN = INVALID_DATA_VALUES_ERR
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN += ' Input secret witness has too large of a bound or weight, or has incorrect '
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN += 'length.'


class OneTimeSecretWitness(object):
    secpar: int
    lp: LatticeParameters
    key: PolynomialVector

    def __init__(self, secpar: int, lp: LatticeParameters, key: PolynomialVector):
        if secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.key = key
        for i in self.key.entries:
            i.const_time_flag = True

    def __eq__(self, other):
        same_secpar: bool = self.secpar == other.secpar
        same_lp: bool = self.lp == other.lp
        same_key: bool = self.key == other.key
        return same_secpar and same_lp and same_key

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.key)


PUBSTAT_INPUT_WITNESS_WRONG_TYPE_ERR = INCORRECT_DATA_TYPE_ERR
PUBSTAT_INPUT_WITNESS_WRONG_TYPE_ERR += ' Input witness must be a SecretWitness object.'
PUBSTAT_LP_OR_SECPAR_MISMATCH_ERR = DATA_MISMATCH_ERR + ' Input LatticeParameters or security parameters do not match.'
PUBSTAT_NEED_CHALL = MISSING_DATA_ERR + ' Need a key challenge if instantiating PublicStatement with a SecretWitness.'
PUBSTAT_NEED_POLYVEC_CHALL = INCORRECT_DATA_TYPE_ERR
PUBSTAT_NEED_POLYVEC_CHALL += ' If instantiating a PublicStatement with a SecretWitness, need a PolynomialVector key '
PUBSTAT_NEED_POLYVEC_CHALL += 'challenge.'
PUBSTAT_NEED_POLY = INCORRECT_DATA_TYPE_ERR + ' Input key must be a polynomial.'

PUBSTAT_NEED_KEY_ERR: str = 'Must instantiate a PublicStatement with a witness and a key challenge, or with a key.'


class OneTimePublicStatement(object):
    secpar: int
    lp: LatticeParameters
    key: Polynomial

    def __init__(self, secpar: int, lp: LatticeParameters, key: Polynomial):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatticeParameters):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatticeParameters' +
                             f' but had {type(lp)}.')
        elif not isinstance(key, Polynomial):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input key must be Polynomial but had {type(key)}.')
        elif key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.key = key
        self.key.const_time_flag = False

    def __eq__(self, other):
        secpars_match = self.secpar == other.secpar
        lps_match = self.lp == other.lp
        keys_match = self.key == other.key
        return secpars_match and lps_match and keys_match

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.key)


SK_NEED_SEED_OR_PAIR = MISSING_DATA_ERR + ' Need a SecretSeed or a left-and-right PolynomialVector pair.'
SK_NEED_SECRETSEED = INCORRECT_DATA_TYPE_ERR + ' If instantiating with a seed, need a SecretSeed object.'
SK_NEED_TWO_POLYVEC = INCORRECT_DATA_TYPE_ERR + ' If instantiating with a pair of polynomials, input them as a '
SK_NEED_TWO_POLYVEC += 'length-two list of PolynomialVectors.'
SK_KEY_LP_MISMATCH = DATA_MISMATCH_ERR + ' Input LatticeParameters objects or security parameter integers do not match.'


class OneTimeSigningKey(object):
    secpar: int
    lp: LatticeParameters
    left_key: PolynomialVector
    right_key: PolynomialVector

    def __init__(self, secpar: int, lp: LatticeParameters, left_key: PolynomialVector, right_key: PolynomialVector):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatticeParameters):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be a LatticeParameters' +
                             f' object, but had {type(lp)}.')
        elif not isinstance(left_key, PolynomialVector) or not isinstance(right_key, PolynomialVector):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Both input keys must be PolynomialVectors.')
        elif left_key.lp != lp or right_key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.left_key = left_key
        self.left_key.const_time_flag = True
        self.right_key = right_key
        self.right_key.const_time_flag = True

    def __getitem__(self, item: int):
        if item not in [0, 1]:
            raise ValueError('Can only get two items.')
        elif item:
            return self.right_key
        return self.left_key

    def __eq__(self, other):
        secpars_match = self.secpar == other.secpar
        lps_match = self.lp == other.lp
        left_match = self.left_key == other.left_key
        right_match = self.right_key == other.right_key
        return secpars_match and lps_match and left_match and right_match

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.left_key) and bool(self.right_key)


VK_NEED_SEED_OR_PAIR = MISSING_DATA_ERR + ' Need a SigningKey or a left-and-right Polynomial pair.'
VK_NEED_SK_AND_CH_OR_LEFT_AND_RIGHT = MISSING_DATA_ERR + ' Need a SigningKey and a key challenge, or need a left key '
VK_NEED_SK_AND_CH_OR_LEFT_AND_RIGHT += 'and a right key.'
VK_NEED_SK = INCORRECT_DATA_TYPE_ERR + ' Input sk is not a SigningKey.'
VK_CH_NEED_POLYVEC = INCORRECT_DATA_TYPE_ERR + ' Input key challenge is not a PolynomialVector.'
VK_LP_OR_SECPAR_MISMATCH = SK_KEY_LP_MISMATCH
VK_NEED_POLYVEC_LR = INCORRECT_DATA_TYPE_ERR + ' Left and right keys must be PolynomialVectors.'
VK_LENGTH_MISMATCH = DATA_MISMATCH_ERR + ' PolynomialVector length mismatch.'


class OneTimeVerificationKey(object):
    secpar: int
    lp: LatticeParameters
    left_key: Polynomial
    right_key: Polynomial

    def __init__(self, secpar: int, lp: LatticeParameters, left_key: Polynomial, right_key: Polynomial):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatticeParameters):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatticeParameters, but' +
                             f' had {type(lp)}.')
        elif not isinstance(left_key, Polynomial) or not isinstance(right_key, Polynomial):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Both input keys must be Polynomial but' +
                             f' had {type(left_key)} and {type(right_key)}.')
        elif left_key.lp != lp or right_key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.left_key = left_key
        self.left_key.const_time_flag = False
        self.right_key = right_key
        self.right_key.const_time_flag = False

    def __getitem__(self, item: int):
        if item not in [0, 1]:
            raise ValueError('Can only get two items.')
        elif item:
            return self.right_key
        return self.left_key

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.left_key) and bool(self.right_key)

    def __eq__(self, other):
        secpars_match = self.secpar == other.secpar
        lps_match = self.lp == other.lp
        left_keys_match = self.left_key == other.left_key
        right_keys_match = self.right_key == other.right_key
        return secpars_match and lps_match and left_keys_match and right_keys_match


OneTimeKeyTuple = Tuple[SecretSeed, OneTimeSigningKey, OneTimeVerificationKey]
ALLOWABLE_DISTRIBUTIONS = [UNIFORM_INFINITY_WEIGHT]


def bits_per_index_set(secpar: int, degree: int, wt: int) -> int:
    """Return number of bits required to sample wt indices uniformly without replacement from [0, 1, ..., degree - 1]
    for a power-of-two degree, with a bias that is O(2**-secpar).
    """
    return ceil(log2(degree)) + (wt - 1) * (ceil(log2(degree)) + secpar)


def bits_per_coefficient(secpar: int, bd: int) -> int:
    """Return number of bits required to sample a coefficient uniformly from [-bd, -bd + 1, ..., bd - 1, bd] with a
    bias that is O(2**-secpar).
    """
    if bd <= 0:
        raise ValueError('Cannot compute bits per coefficient for a non-positive bound bd.')
    return ceil(log2(bd)) + 1 + secpar


class SchemeParameters(object):
    secpar: int
    lp: LatticeParameters
    key_ch: PolynomialVector
    distribution: str

    def __init__(self, secpar: int, lp: LatticeParameters, distribution: str, key_ch: PolynomialVector = None):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatticeParameters):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatticeParameters.')
        elif key_ch is not None and not isinstance(key_ch, PolynomialVector):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input key challenge must be a PolynomialVector or None.')
        elif not isinstance(distribution, str) or distribution not in ALLOWABLE_DISTRIBUTIONS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input distribution must be a string code indicating' +
                             f' a supported distribution.')
        elif key_ch is not None and key_ch.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.distribution = distribution
        if key_ch is not None:
            self.key_ch = key_ch
            self.key_ch.const_time_flag = False
        elif distribution == UNIFORM_INFINITY_WEIGHT:
            self.key_ch = random_polynomialvector(
                secpar=secpar, lp=lp, distribution=distribution, dist_pars={'bd': lp.modulus//2, 'wt': lp.degree},
                bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=lp.degree),
                btd=bits_per_coefficient(secpar=secpar, bd=lp.modulus // 2),
                const_time_flag=True, num_coefs=lp.degree
            )
        else:
            raise ValueError('Unsupported distribution.')

    def __eq__(self, other) -> bool:
        same_secpar = self.secpar == other.secpar
        same_lp = self.lp == other.lp
        same_ch = self.key_ch == other.key_ch
        same_dist = self.distribution == other.distribution
        return same_secpar and same_lp and same_ch and same_dist


def make_random_seed(secpar: SecurityParameter, pp: PublicParameters) -> SecretSeed:
    seed = bin(randbelow(2 ** secpar))[2:].zfill(secpar)
    return SecretSeed(secpar=secpar, lp=pp['scheme_parameters'].lp, seed=seed)


def make_one_key(pp: PublicParameters, seed: SecretSeed = None) -> OneTimeKeyTuple:
    secpar = pp['scheme_parameters'].secpar
    distribution = pp['scheme_parameters'].distribution
    lp = pp['scheme_parameters'].lp
    x = seed
    if not x:
        x = make_random_seed(secpar=secpar, pp=pp)
    left_signing_key: PolynomialVector = hash2polynomialvector(
        secpar=secpar, lp=lp,
        distribution=distribution,
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
        distribution=distribution,
        dist_pars={'bd': pp['sk_bd'], 'wt': pp['sk_wt']},
        num_coefs=pp['sk_wt'],
        bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=pp['sk_wt']),
        btd=bits_per_coefficient(secpar=secpar, bd=pp['sk_bd']),
        salt=pp['sk_salt'] + 'RIGHT',
        msg=x.seed,
        const_time_flag=True
    )
    otsk = OneTimeSigningKey(secpar=secpar, lp=lp, left_key=left_signing_key, right_key=right_signing_key)
    key_ch = pp['scheme_parameters'].key_ch
    key_ch.const_time_flag = True
    otvk = OneTimeVerificationKey(secpar=secpar, lp=lp, left_key=key_ch * left_signing_key, right_key=key_ch * right_signing_key)
    return x, otsk, otvk


def make_one_key_wrapper(pp: PublicParameters, num_keys: int = 1, seeds: List[SecretSeed] = None) -> List[OneTimeKeyTuple]:
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


def distribute_tasks(tasks: List[Any], num_workers: int = None) -> List[List[Any]]:
    """
    Helper function that distributes a list of arbitrary tasks among a specific number of workers

    :param tasks: iterable containing list of tasks to carry out
    :param num_workers: number of workers available in the pool (usually = number of CPU cores)
    :return: task list broken up into num_workers segments
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


def keygen_core(pp: PublicParameters, num_keys_to_gen: int = 1, seeds: List[SecretSeed] = None,
                multiprocessing: bool = None) -> List[OneTimeKeyTuple]:
    # TODO: Move to one_time_keys.py, import from one_time_keys, refactor so this is a wrapper
    """ Wraps make_one_key_wrapper to handle workload distribution for batch generation """

    # Only default to parallelization if more than a few keys are needed (to avoid unnecessary overhead)
    if multiprocessing is None:
        multiprocessing: bool = num_keys_to_gen >= 16
    num_workers: int = cpu_count()

    # Pass straight through to make_one_key_wrapper() if there is no reason or desire to parallelize (to avoid extra overhead)
    if (not multiprocessing) or (num_keys_to_gen == 1) or (num_workers == 1):
        return make_one_key_wrapper(pp=pp, num_keys=num_keys_to_gen, seeds=seeds)

    # Prepare inputs for the pool
    if not seeds:
        iterable: List[Tuple[Dict[str, Any], int]] = [(pp, ceil(num_keys_to_gen / num_workers))] * num_workers
    else:
        seed_batches: List[List[Any]] = distribute_tasks(tasks=seeds)
        iterable: List[tuple] = list(zip([pp] * len(seed_batches), [len(x) for x in seed_batches], seed_batches))

    # Generate the keys and return the flattened results
    with Pool(num_workers) as pool:
        nested_keys: List[List[OneTimeKeyTuple]] = pool.starmap(func=keygen_core, iterable=iterable)
    return [item for sublist in nested_keys for item in sublist][:num_keys_to_gen]


def challenge_core(pp: PublicParameters, otvk: OneTimeVerificationKey, msg: Message) -> Challenge:
    return hash2polynomial(
        secpar=pp['scheme_parameters'].secpar,
        lp=pp['scheme_parameters'].lp,
        distribution=pp['scheme_parameters'].distribution,
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


def sign_core(pp: PublicParameters, otk: OneTimeKeyTuple, msg: Message) -> Signature:
    c: Challenge = challenge_core(pp=pp, otvk=otk[2], msg=msg)
    signature: Signature = otk[1][0] ** c + otk[1][1]
    return signature


def verify_core(pp: PublicParameters, otvk: OneTimeVerificationKey, msg: Message, sig: Signature) -> bool:
    sig.const_time_flag = False
    cnws = sig.get_coef_rep()
    n, w = max(i[1] for i in cnws), max(i[2] for i in cnws)
    if n > pp['vf_bd'] or w > pp['vf_wt']:
        return False

    key_ch = pp['scheme_parameters'].key_ch
    c: Challenge = challenge_core(pp=pp, otvk=otvk, msg=msg)

    key_ch.const_time_flag = False
    c.const_time_flag = False
    otvk.left_key.const_time_flag = False
    otvk.right_key.const_time_flag = False

    lhs = key_ch * sig
    rhs = otvk[0] * c + otvk[1]

    return lhs == rhs