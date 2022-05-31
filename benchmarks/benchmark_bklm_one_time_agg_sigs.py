"""
We benchmark the lattice-crypto.lm_one_time_sigs.lm_one_time_sigs module.
"""
from timeit import default_timer as timer
from lattice_cryptography.bklm_one_time_agg_sigs import *
from multiprocessing import Pool, cpu_count
from math import ceil
from typing import List, Any
from secrets import randbits

# Benchmarking params
SAMPLE_SIZE: int = 2 ** 6
multiprocessing: bool = False  # << Set to False to disable multiprocessing (WIP)

# Parallelization params
num_cores: int = 64
if multiprocessing and num_cores > 1:
    num_workers = min((num_cores, cpu_count()))
    sample_size_per_core: int = ceil(SAMPLE_SIZE / num_workers)
    print(f"Beginning benchmarks (multiprocessing keygen with {num_workers} workers)")
else:
    num_workers = sample_size_per_core = 0  # N/A
    print(f"Beginning benchmarks (multiprocessing disabled)")


# Helper function
def flatten(some_nested_list: List[List[Any]]) -> List[Any]:
    return [item for sublist in some_nested_list for item in sublist]

time_per_key = 0.0
const_time = 0.0
# Benchmark parameter setup
for secpar, capacity in ALLOWABLE_PARAMETERS:
    # Make setup parameters
    print(f"Simple benchmarking for secpar, capacity = {secpar, capacity}.")
    print(f"\tParameter generation benchmarking for secpar, capacity = {secpar, capacity}.")
    print(f"\t\tGenerating parameters.  ")
    start = timer()
    pp = make_setup_parameters(secpar=secpar, carrying_capacity=capacity)
    end = timer()
    const_time += end - start
    print(f"\t\tElapsed time = {end - start}.")

    keygen_time: float = 0.0
    msg_sample_time: float = 0.0
    sign_time: float = 0.0
    vf_time: float = 0.0
    agg_time: float = 0.0
    agvf_time: float = 0.0

    some_key_samples: list[list[OneTimeKeyTuple]] = []
    some_msg_samples: list[list[Msg]] = []

    for i in range(SAMPLE_SIZE):
        print(f"\tKey generation for {i}-th sample with secpar, capacity = {secpar, capacity}.")
        print(f"\t\tGenerating {capacity} keys.")
        start = timer()
        some_key_samples += keygen(pp=pp, num=capacity, multiprocessing=multiprocessing)
        end = timer()
        keygen_time += end - start
        print(f"\t\tElapsed time = {end - start}, averaging {(end - start) / capacity} per item.")

        print(f"\tSampling messages to sign for {i}-th sample.")
        start = timer()
        some_msg_samples += [bin(randbits(256))[2:].zfill(256) for _ in range(capacity)]
        end = timer()
        msg_sample_time += end - start
        print(f"\t\tElapsed time = {end - start}, averaging {(end - start) / capacity} per item.")

        print(f"\tSignature benchmarking for {i}-th sample with secpar, capacity = {secpar, capacity}.")
        sign_multiprocessing: bool = False  # STILL A WIP - MAY NOT WORK
        sign_input_tuples = [[(pp, key, msg) for key, msg in zip(next_keys, next_msgs)] for next_keys, next_msgs in zip(some_key_samples, some_msg_samples)]
        if sign_multiprocessing and num_cores > 1:
            start = timer()
            with Pool(num_workers) as pool:
                some_sigs_for_keys_without_seeds = pool.starmap(func=sign, iterable=sign_input_tuples)
            end = timer()
        else:
            start = timer()
            some_sigs_for_keys_without_seeds = [[sign(*args) for args in sign_input_tuple] for sign_input_tuple in sign_input_tuples]
            end = timer()
        sign_time += end - start
        print(f"\t\tElapsed time = {end - start}, averaging {(end - start) / capacity} per item.")

        print(f"\tVerification benchmarking for {i}-th sample with secpar, capacity = {secpar, capacity}")
        verify_multiprocessing = False  # STILL A WIP - MAY NOT WORK
        verify_input_tuples = [[(pp, [otk[2]], [m], sig) for pp, otk, m, sig in zip(x, y, z, w)] for x, y, z, w in zip([pp] * len(some_key_samples), some_key_samples, some_msg_samples, some_sigs_for_keys_without_seeds)]
        if verify_multiprocessing and num_cores > 1:
            start = timer()
            with Pool(num_workers) as pool:
                results_without_seeds = pool.starmap(func=verify, iterable=verify_input_tuples)
            end = timer()
        else:
            start = timer()
            results_without_seeds = [[verify(*args) for args in verify_input_tuple] for verify_input_tuple in verify_input_tuples]
            end = timer()
        vf_time += end - start
        print(f"\t\tElapsed time = {end - start}, averaging {time_per_key} per item.")

        # TODO: Parallelize here, Mitchell?
        print(f"\tAggregation benchmarking for {i}-th sample with secpar, capacity = {secpar, capacity}")
        aggregate_input_tuples = [[pp]]
        start = timer()
        results_without_seeds = [[aggregate(*args) for args in verify_input_tuple] for verify_input_tuple in verify_input_tuples]
        end = timer()
        vf_time += end - start
        print(f"\t\tElapsed time = {end - start}, averaging {time_per_key} per item.")
