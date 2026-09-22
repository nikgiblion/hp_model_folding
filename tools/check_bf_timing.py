"""
Timing of the exact HP ground state search 
on a 20-mer, with and without pruning.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import time

from bf import exact_ground_states

SEQUENCE = "HPHPPHHPHPPHPHHPPHPH"

for prune in (False, True):
    start = time.perf_counter()
    e, optimal = exact_ground_states(SEQUENCE, prune=prune)
    elapsed = time.perf_counter() - start
    label = "branch and bound" if prune else "full enumeration"
    print(f"{label:>18}: E = {e}, {len(optimal)} ground states, {elapsed:.1f} s")