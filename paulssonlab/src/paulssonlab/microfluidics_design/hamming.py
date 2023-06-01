"""Implements Hamming encoding/decoding with or without an extra SECDED bit.

Inspired by https://github.com/dominiccarrano/hamming.
"""

from math import ceil, floor, log2

import numpy as np
from bitarray import bitarray


class DoubleError(Exception):
    pass


def encode(data: bitarray, secded=True):
    if not len(data):
        raise ValueError("data must be non-empty")
    if secded:
        prefix = 1
    else:
        prefix = 0
    num_parity_bits = _num_parity_bits_needed(len(data))
    encoded_length = len(data) + num_parity_bits + prefix
    encoded = bitarray(encoded_length)
    for parity_bit_index in _powers_of_two(num_parity_bits):
        encoded[parity_bit_index - 1 + prefix] = _calculate_parity(
            data, parity_bit_index
        )
    data_index = 0
    for encoded_index in range(prefix, len(encoded)):
        if not _is_power_of_two(encoded_index + 1 - prefix):
            encoded[encoded_index] = data[data_index]
            data_index += 1

    if secded:
        encoded[0] = encoded[prefix:].count() % 2
    return encoded


def decode(encoded: bitarray, secded=True):
    if not len(encoded):
        raise ValueError("encoded must be non-empty")
    if secded:
        prefix = 1
    else:
        prefix = 0
    encoded_length = len(encoded)
    num_parity_bits = int(floor(log2(encoded_length - prefix)) + 1)
    error_idx = 0
    decoded = _data_bits(encoded[prefix:])
    if secded:
        extra_parity_matches = encoded[0] == encoded[prefix:].count() % 2
    else:
        extra_parity_matches = True
    for parity_bit_index in _powers_of_two(num_parity_bits):
        expected = _calculate_parity(decoded, parity_bit_index)
        actual = encoded[parity_bit_index - 1 + prefix]
        if not expected == actual:
            error_idx += parity_bit_index
    if secded and error_idx and extra_parity_matches:
        raise DoubleError
    elif error_idx:
        # correct single errors
        encoded[error_idx - 1 + prefix] ^= 1
        decoded = _data_bits(encoded[prefix:])
    return decoded


def _num_parity_bits_needed(length: int):
    n = _next_power_of_two(length)
    lower_bin = floor(log2(n))
    upper_bin = lower_bin + 1
    data_bit_boundary = n - lower_bin - 1
    return lower_bin if length <= data_bit_boundary else upper_bin


def _calculate_parity(data: bitarray, parity: int):
    return (data & _parity_mask(parity, len(data))).count() % 2


def _parity_mask(parity: int, length: int):
    idxs = np.arange(length + _num_parity_bits_needed(length)) + 1
    mask = (idxs & parity != 0)[~_is_power_of_two(idxs)][:length]
    return bitarray(list(mask))


def _data_bits(encoded: bitarray):
    idxs = np.arange(len(encoded))
    ary = np.array(list(encoded))[~_is_power_of_two(idxs + 1)]
    return bitarray(list(ary))


def _next_power_of_two(n: int):
    if _is_power_of_two(n):
        return n << 1
    else:
        return 2 ** ceil(log2(n))


def _is_power_of_two(n: int):
    return np.logical_and(np.logical_not(n == 0), ((n & (n - 1)) == 0))


def _powers_of_two(n: int):
    return 1 << np.arange(n)
