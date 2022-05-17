"""
The lattice_cryptography.keys module handles keys and their generation.
"""
from lattice_algebra import \
    UNIFORM_INFINITY_WEIGHT, \
    Polynomial, \
    PolynomialVector, \
    LatticeParameters as LatPar, \
    random_polynomialvector, \
    is_bitstring, \
    hash2polynomial, \
    hash2polynomialvector
from math import ceil, log2
from typing import Dict, Any, Tuple, List
from multiprocessing import Pool, cpu_count
from secrets import randbelow

# Allowed security parameters
ALLOWABLE_SECPARS = [128, 256]

# Typing
SecPar = int
PubPar = Dict[str, Any]
Msg = str
Chall = Polynomial
Sig = PolynomialVector

# Error messages
GENERIC_ERR: str = 'Something went wrong.'
MISSING_DATA_ERR: str = 'Missing some required data.'
INCORRECT_DATA_TYPE_ERR: str = 'Required input data not the correct type.'
DATA_MISMATCH_ERR: str = 'Input data did not match.'
SEED_INST_ERR_NEED_BITS: str = INCORRECT_DATA_TYPE_ERR + ' Input must be a binary string.'
INVALID_DATA_VALUES_ERR: str = 'Required input data does not have valid values.'


class SecretSeed(object):
    secpar: SecPar
    lp: LatPar
    seed: Msg

    def __init__(self, seed: Msg, secpar: SecPar, lp: LatPar):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be an integer in' +
                             f' {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not is_bitstring(seed):
            raise ValueError(SEED_INST_ERR_NEED_BITS)
        elif not isinstance(lp, LatPar):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be' +
                             f' LatPar object.')
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
SECWIT_INST_ERR_LP_MISMATCH += ' Input LatPar object does not match the LatPar for the input '
SECWIT_INST_ERR_LP_MISMATCH += 'PolynomialVector.'
SECWIT_INST_ERR_SECPAR_MISMATCH = DATA_MISMATCH_ERR
SECWIT_INST_ERR_SECPAR_MISMATCH += ' The input LatPar object or the input secpar integer does not match '
SECWIT_INST_ERR_SECPAR_MISMATCH += 'the input SecretSeed.'
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN = INVALID_DATA_VALUES_ERR
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN += ' Input secret witness has too large of a bound or weight, or has incorrect '
SECWIT_INVALID_BOUND_WEIGHT_OR_LEN += 'length.'


class OneTimeSecretWitness(object):
    secpar: int
    lp: LatPar
    key: PolynomialVector

    def __init__(self, secpar: int, lp: LatPar, key: PolynomialVector, const_time_flag: bool = True):
        if secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.key = key
        for i in self.key.entries:
            i.const_time_flag = const_time_flag

    def __eq__(self, other):
        same_secpar: bool = self.secpar == other.secpar
        same_lp: bool = self.lp == other.lp
        same_key: bool = self.key == other.key
        return same_secpar and same_lp and same_key

    def __bool__(self):
        return bool(self.secpar) and bool(self.lp) and bool(self.key)


PUBSTAT_INPUT_WITNESS_WRONG_TYPE_ERR = INCORRECT_DATA_TYPE_ERR
PUBSTAT_INPUT_WITNESS_WRONG_TYPE_ERR += ' Input witness must be a SecretWitness object.'
PUBSTAT_LP_OR_SECPAR_MISMATCH_ERR = DATA_MISMATCH_ERR + ' Input LatPar or security parameters do not match.'
PUBSTAT_NEED_CHALL = MISSING_DATA_ERR + ' Need a key challenge if instantiating PublicStatement with a SecretWitness.'
PUBSTAT_NEED_POLYVEC_CHALL = INCORRECT_DATA_TYPE_ERR
PUBSTAT_NEED_POLYVEC_CHALL += ' If instantiating a PublicStatement with a SecretWitness, need a PolynomialVector key '
PUBSTAT_NEED_POLYVEC_CHALL += 'challenge.'
PUBSTAT_NEED_POLY = INCORRECT_DATA_TYPE_ERR + ' Input key must be a polynomial.'

PUBSTAT_NEED_KEY_ERR: str = 'Must instantiate a PublicStatement with a witness and a key challenge, or with a key.'


class OneTimePublicStatement(object):
    secpar: int
    lp: LatPar
    key: Polynomial

    def __init__(self, secpar: int, lp: LatPar, key: Polynomial):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatPar):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatPar' +
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
SK_KEY_LP_MISMATCH = DATA_MISMATCH_ERR + ' Input LatPar objects or security parameter integers do not match.'


class OneTimeSigningKey(object):
    secpar: int
    lp: LatPar
    left_key: PolynomialVector
    right_key: PolynomialVector

    def __init__(self, secpar: int, lp: LatPar, left_key: PolynomialVector, right_key: PolynomialVector,
                 const_time_flag: bool = True):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatPar):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be a LatPar' +
                             f' object, but had {type(lp)}.')
        elif not isinstance(left_key, PolynomialVector) or not isinstance(right_key, PolynomialVector):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Both input keys must be PolynomialVectors.')
        elif left_key.lp != lp or right_key.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        elif const_time_flag != left_key.const_time_flag or const_time_flag != right_key.const_time_flag:
            raise ValueError('Must __init__ with all-same const_time_flag.')
        self.secpar = secpar
        self.lp = lp
        self.left_key = left_key
        self.right_key = right_key

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
    lp: LatPar
    left_key: Polynomial
    right_key: Polynomial

    def __init__(self, secpar: int, lp: LatPar, left_key: Polynomial, right_key: Polynomial):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatPar):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatPar, but' +
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


SS = SecretSeed
OTSW = OneTimeSecretWitness
OTSK = OneTimeSigningKey
OTPS = OneTimePublicStatement
OTVK = OneTimeVerificationKey
OneTimeKeyTuple = Tuple[SS, OTSK, OTVK]
OTK = OneTimeKeyTuple
ALLOWABLE_DISTRIBUTIONS = [UNIFORM_INFINITY_WEIGHT]


def bits_per_index_set(secpar: int, degree: int, wt: int) -> int:
    """Return number of bits required to sample wt indices uniformly without replacement from [0, 1, ..., degree - 1]
    for a power-of-two degree, with a bias that is O(2**-secpar).
    """
    if min(degree, secpar, wt) < 1:
        raise ValueError('Cannot bits_per_index_set without strictly positive integer inputs.')
    return ceil(log2(degree)) + (wt - 1) * (ceil(log2(degree)) + secpar)


def bits_per_coefficient(secpar: int, bd: int) -> int:
    """Return number of bits required to sample a coefficient uniformly from [-bd, -bd + 1, ..., bd - 1, bd] with a
    bias that is O(2**-secpar).
    """
    if min(secpar, bd) < 1:
        raise ValueError('Cannot compute bits per coefficient for a non-positive secpar and bd.')
    return ceil(log2(bd)) + 1 + secpar


class SchemeParameters(object):
    secpar: int
    lp: LatPar
    key_ch: PolynomialVector
    dstn: str

    def __init__(self, secpar: int, lp: LatPar, distribution: str, key_ch: PolynomialVector = None):
        if not isinstance(secpar, int) or secpar not in ALLOWABLE_SECPARS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input security parameter must be' +
                             f' in {ALLOWABLE_SECPARS} but had {secpar}.')
        elif not isinstance(lp, LatPar):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input lattice parameters must be LatPar.')
        elif key_ch is not None and not isinstance(key_ch, PolynomialVector):
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input key challenge must be a PolynomialVector or None.')
        elif not isinstance(distribution, str) or distribution not in ALLOWABLE_DISTRIBUTIONS:
            raise ValueError(INVALID_DATA_VALUES_ERR + f' Input distribution must be a string code indicating' +
                             f' a supported distribution.')
        elif key_ch is not None and key_ch.lp != lp:
            raise ValueError(SECWIT_INST_ERR_LP_MISMATCH)
        self.secpar = secpar
        self.lp = lp
        self.dstn = distribution
        if key_ch is not None:
            self.key_ch = key_ch
            self.key_ch.const_time_flag = True
        elif distribution == UNIFORM_INFINITY_WEIGHT:
            self.key_ch = random_polynomialvector(
                secpar=secpar, lp=lp, distribution=distribution, dist_pars={'bd': lp.modulus // 2, 'wt': lp.degree},
                bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=lp.degree),
                btd=bits_per_coefficient(secpar=secpar, bd=lp.modulus // 2),
                const_time_flag=True, num_coefs=lp.degree
            )
        else:
            raise ValueError('Unsupported distribution.')

    def __eq__(self, other) -> bool:
        same_secpar: bool = self.secpar == other.secpar
        same_lp: bool = self.lp == other.lp
        same_ch: bool = self.key_ch == other.key_ch
        same_dist: bool = self.dstn == other.dstn
        return same_secpar and same_lp and same_ch and same_dist


def make_random_seed(secpar: SecPar, pp: PubPar) -> SS:
    if secpar not in ALLOWABLE_SECPARS:
        raise ValueError(f'Cannot make_random_seed without secpar in ALLOWABLE_SECPARS, but secpar={secpar}.')
    seed: Msg = bin(randbelow(2 ** secpar))[2:].zfill(secpar)
    return SS(secpar=secpar, lp=pp['scheme_parameters'].lp, seed=seed)


NEED_BD_WT_SALT_ERR = 'Need scheme_parameters, sk_bd, sk_wt, and sk_salt in input pp for make_one_key.'


def make_one_key(pp: PubPar, seed: SS = None, const_time_flag: bool = True) -> OTK:
    if any(x not in pp for x in ['scheme_parameters', 'sk_bd', 'sk_wt', 'sk_salt']):
        raise ValueError(NEED_BD_WT_SALT_ERR)
    secpar: SecPar = pp['scheme_parameters'].secpar
    dstn: str = pp['scheme_parameters'].dstn
    lp: LatPar = pp['scheme_parameters'].lp
    x: SS = seed
    if not x:
        x: SS = make_random_seed(secpar=secpar, pp=pp)
    left_signing_key: PolynomialVector = hash2polynomialvector(
        secpar=secpar, lp=lp,
        distribution=dstn,
        dist_pars={'bd': pp['sk_bd'], 'wt': pp['sk_wt']},
        num_coefs=pp['sk_wt'],
        bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=pp['sk_wt']),
        btd=bits_per_coefficient(secpar=secpar, bd=pp['sk_bd']),
        salt=pp['sk_salt'] + 'LEFT',
        msg=x.seed,
        const_time_flag=const_time_flag
    )
    right_signing_key: PolynomialVector = hash2polynomialvector(
        secpar=secpar, lp=lp,
        distribution=dstn,
        dist_pars={'bd': pp['sk_bd'], 'wt': pp['sk_wt']},
        num_coefs=pp['sk_wt'],
        bti=bits_per_index_set(secpar=secpar, degree=lp.degree, wt=pp['sk_wt']),
        btd=bits_per_coefficient(secpar=secpar, bd=pp['sk_bd']),
        salt=pp['sk_salt'] + 'RIGHT',
        msg=x.seed,
        const_time_flag=const_time_flag
    )
    otsk: OTSK = OTSK(secpar=secpar, lp=lp, left_key=left_signing_key, right_key=right_signing_key,
                      const_time_flag=const_time_flag)
    key_ch: PolynomialVector = pp['scheme_parameters'].key_ch
    key_ch.const_time_flag = const_time_flag
    otvk: OTVK = OTVK(secpar=secpar, lp=lp, left_key=key_ch * left_signing_key, right_key=key_ch * right_signing_key)
    return x, otsk, otvk


def make_one_key_wrapper(pp: PubPar, num: int = 1, seeds: List[SecretSeed] = None) -> List[OneTimeKeyTuple]:
    if num < 1:
        raise ValueError('Can only generate a natural number worth of keys.')
    elif seeds is not None and len(seeds) != num:
        raise ValueError('Must either roll keys with no seeds, or with a seed for each key.')
    elif seeds is None and num == 1:
        return [make_one_key(pp=pp)]
    elif seeds is not None and num == 1:
        return [make_one_key(pp=pp, seed=seeds[0])]
    elif seeds is None:
        return [make_one_key(pp=pp) for _ in range(num)]
    return [make_one_key(pp=pp, seed=next_seed) for next_seed in seeds]


def distribute_tasks(tasks: List[Any], num_workers: int = None) -> List[List[Any]]:
    """
    Helper function that distributes a list of arbitrary tasks among a specific number of workers.

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


def keygen_core(pp: PubPar, num: int = 1, seeds: List[SecretSeed] = None,
                multiprocessing: bool = None) -> List[OneTimeKeyTuple]:
    """ Wraps make_one_key_wrapper to handle workload distribution for batch generation """

    # Only default to parallelization if more than a few keys are needed (to avoid unnecessary overhead)
    if multiprocessing is None:
        multiprocessing: bool = num >= 16
    num_workers: int = cpu_count()

    # Pass to make_one_key_wrapper() to not parallelize to avoid extra overhead
    if (not multiprocessing) or (num == 1) or (num_workers == 1):
        return make_one_key_wrapper(pp=pp, num=num, seeds=seeds)

    # Prepare inputs for the pool
    if not seeds:
        iterable: List[Tuple[Dict[str, Any], int]] = [(pp, ceil(num / num_workers))] * num_workers
    else:
        seed_batches: List[List[Any]] = distribute_tasks(tasks=seeds)
        iterable: List[tuple] = list(zip([pp] * len(seed_batches), [len(x) for x in seed_batches], seed_batches))

    # Generate the keys and return the flattened results
    with Pool(num_workers) as pool:
        nested_keys: List[List[OneTimeKeyTuple]] = pool.starmap(func=keygen_core, iterable=iterable)
    return [item for sublist in nested_keys for item in sublist][:num]


def challenge_core(pp: PubPar, otvk: OneTimeVerificationKey, msg: Msg,
                   const_time_flag: bool = True) -> Chall:
    return hash2polynomial(
        secpar=pp['scheme_parameters'].secpar,
        lp=pp['scheme_parameters'].lp,
        distribution=pp['scheme_parameters'].dstn,
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
        const_time_flag=const_time_flag
    )


def sign_core(pp: PubPar, otk: OneTimeKeyTuple, msg: Msg) -> Sig:
    if otk[1].left_key.const_time_flag != otk[1].right_key.const_time_flag:
        raise ValueError('Must either do all constant time or not, no mixing.')
    c: Chall = challenge_core(pp=pp, otvk=otk[2], msg=msg, const_time_flag=otk[1].left_key.const_time_flag)
    signature: Sig = otk[1].left_key ** c + otk[1].right_key
    return signature


def verify_core(sig: Sig, bd: int, wt: int, key_ch: PolynomialVector, target: Polynomial) -> bool:
    for entry in sig.entries:
        entry.const_time_flag = False
    key_ch.const_time_flag = False

    norms_and_wts = sig.get_coef_rep()
    n, w = max(i[1] for i in norms_and_wts), max(i[2] for i in norms_and_wts)
    if n < 1 or n > bd or wt < 1 or w > wt:
        return False
    return key_ch * sig == target
