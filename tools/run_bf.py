"""
Exact ground-state search for an arbitrary HP sequence.

Usage from the repository root:
    python tools/run_bf.py --sequence HPHPPHHPHPPHPHHPPHPH
    python tools/run_bf.py --sequence "H2(P2H)3H" --no-prune
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import time

from bf import exact_distinct_ground_states, exact_ground_states
from conformation import sequence_check
from energy import energy
from visualization import ascii_vizualization

DEFAULT_SEQUENCE = "HPHPPHHPHPPHPHHPPHPH"  #known benchmark sequence with E = -9


def parse_args():
    parser = argparse.ArgumentParser(
        description="Exact HP ground states by brute-force enumeration")
    parser.add_argument("--sequence", default=DEFAULT_SEQUENCE,
                        help="HP sequence, plain or compressed, e.g. H2(P2H)7H")
    parser.add_argument("--no-prune", action="store_true",
                        help="disable branch and bound (full enumeration, much slower)")
    parser.add_argument("--keep-reversed", action="store_true",
                        help="report D4 classes only, without folding in chain reversal")
    parser.add_argument("--max-shown", type=int, default=5,
                        help="max ground-state structures printed (0 for none)")
    return parser.parse_args()


def main():
    args = parse_args()
    sequence = sequence_check(args.sequence)
    search = exact_ground_states if args.keep_reversed else exact_distinct_ground_states

    start = time.perf_counter()
    e, optimal = search(sequence, prune=not args.no_prune)
    elapsed = time.perf_counter() - start

    mode = "full enumeration" if args.no_prune else "branch and bound"
    symmetry = "D4 only" if args.keep_reversed else "D4 + chain reversal"
    print(f"sequence   : {sequence} (n = {len(sequence)})")
    print(f"method     : {mode}, {symmetry}")
    print(f"energy     : {e}")
    print(f"ground     : {len(optimal)} structures")
    print(f"runtime    : {elapsed:.2f} s")

    for i, coords in enumerate(optimal[:args.max_shown], start=1):
        print(f"\n--- ground state {i}: E = {energy(sequence, coords)} ---")
        print(ascii_vizualization(sequence, coords))
    if len(optimal) > args.max_shown:
        print(f"\n... {len(optimal) - args.max_shown} more not shown")


if __name__ == "__main__":
    main()