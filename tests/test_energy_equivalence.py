"""
Test to check energy equivalence for two different energy calculators.
"""
import random
import struct
import pytest
from conformation import (linear_conformation, random_conformation, random_move, self_avoiding_walking_check)
from energy import (all_contacts, all_contacts_old, energy_old, energy)

SEQ_LENGTHS = list(range(4, 53, 4)) #10 sequences
PER_LENGTH = 25
SEED = 0


#To avoid 0.0 and -0.0 problem.
def bits(x):
    return struct.pack("<d", x)


#Random self-avoiding walks for compact and extended shapes (fresh start and chain of accepted moves).
def sample_conformations(n, rng, count):
    coords = random_conformation(n, rng)
    yield coords
    produced = 1
    while produced < count:
        candidate = random_move(coords, rng)
        if candidate is not None and self_avoiding_walking_check(candidate):
            coords = candidate
            yield coords
            produced += 1


@pytest.mark.parametrize("n", [1, 2, 3])
def test_degenerate_lengths_agree(n):
    # Too short for the random sampler: pivot_move needs an interior site.
    coords = linear_conformation(n)
    assert all_contacts(coords) == all_contacts_old(coords)
    assert bits(energy("H" * n, coords)) == bits(energy_old("H" * n, coords))


@pytest.mark.parametrize("n", SEQ_LENGTHS)
def test_contact_lists_agree(n):
    rng = random.Random(SEED + n)
    for coords in sample_conformations(n, rng, PER_LENGTH):
        assert all_contacts(coords) == all_contacts_old(coords), coords


@pytest.mark.parametrize("n", SEQ_LENGTHS)
def test_energies_agree_bit_for_bit(n):
    rng = random.Random(SEED + n)
    for coords in sample_conformations(n, rng, PER_LENGTH):
        sequences = (
            "".join(rng.choice("HP") for _ in range(n)),
            "H" * n,   # maximum number of contacts
            "P" * n,   # -0.0 would show up
        )
        for sequence in sequences:
            new_fast, old_slow = energy(sequence, coords), energy_old(sequence, coords)
            assert bits(new_fast) == bits(old_slow), (sequence, coords, new_fast, old_slow)