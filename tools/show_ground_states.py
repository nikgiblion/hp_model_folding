"""ASCII pictures of all exact ground states 
of a given HP sequence."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


from bf import exact_ground_states, exact_distinct_ground_states
from energy import energy
from visualization import ascii_vizualization

SEQUENCE = "HPHPPHHPHPPHPHHPPHPH"

e, optimal = exact_distinct_ground_states(SEQUENCE)
print(f"{SEQUENCE}: E = {e}, {len(optimal)} ground states up to D4 symmetry\n")
for i, coords in enumerate(optimal, start=1):
    print(f"--- ground state {i}: E = {energy(SEQUENCE, coords)} ---")
    print(ascii_vizualization(SEQUENCE, coords))
    print()