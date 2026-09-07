"""
Test to check energy equivalence for two different energy calculators
"""

import random 
import struct

from conformation import (linear_conformation, random_conformation, random_move, self_avoiding_walking_check)

from energy import (all_contacts, all_contacts_old, energy_old, energy)

SEQ_LENGTHS = list(range(4, 53, 4)) #10 sequences
PER_LENGTH = 10

#to avoid 0.0 and -0.0 problem
def bits(x):
    return struct.pack("<d", x)

#Random self-avoiding walks for compact and extended shapes (fresh start and chain of accepted moves).
def sample_conformations(rng):
    for n in SEQ_LENGTHS:
        coords = random_conformation(n, rng)
        yield n, coords
        produced = 1
        while produced < PER_LENGTH:
            candidate = random_move(coords, rng)
            if candidate is not None and self_avoiding_walking_check(candidate):
                coords = candidate
                yield n, coords
                produced += 1


def random_sequence(n, rng):
    return "".join(rng.choice("HP") for _ in range(n))

def main():
    rng = random.Random(0)
    checked = 0

    #Degenerate lengths cannot  be reach by random sampler
    for n in (1, 2, 3):
        coords = linear_conformation(n)
        assert all_contacts_old(coords) == all_contacts(coords)
        assert bits(energy("P" * n, coords)) == bits(energy_old("P" * n, coords))
        checked += 1

        for n, coords in sample_conformations(rng):
            assert all_contacts_old(coords) == all_contacts(coords), coords
            for sequence in (random_sequence(n, rng), "H" * n, "P" * n):
                new_fast, old_slow = energy(sequence, coords), energy_old(sequence, coords)
                assert bits(new_fast) == bits(old_slow), (sequence, coords, new_fast, old_slow)

            checked += 1

    print(f"Good, {checked} conformations, lengths {SEQ_LENGTHS[0]}-{SEQ_LENGTHS[-1]}, "
          f"Contacts and energies are identical!")

if __name__ == "__main__":
    main()