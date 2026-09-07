"""
This script is for fast time comparison 
between old and new energy calculations
"""

import random
import timeit
from conformation import random_conformation
from energy import energy, energy_old

SEQ_LENGTHS = (10, 30, 50, 100, 200)
REPEATS = 500

rng = random.Random(0)

print(f"{'n':>4} {'old_slow, us':>11} {'new_fast, us':>10} {'speedup':>9}")

for n in SEQ_LENGTHS:
    sequence = ("HP" * n)[:n]
    coords = random_conformation(n, rng)
    assert energy_old(sequence, coords) == energy(sequence, coords)

    old_slow = timeit.timeit(lambda: energy_old(sequence, coords), number=REPEATS)
    new_fast = timeit.timeit(lambda: energy(sequence, coords), number=REPEATS)
    print(f"{n:>4} {old_slow / REPEATS * 1e6:>11.1f} "
          f"{new_fast / REPEATS * 1e6:>10.1f} {old_slow / new_fast:>8.1f}x")