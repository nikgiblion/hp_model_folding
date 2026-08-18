import argparse
import time
from collections import Counter

from conformation import expand_sequence
from annealing import multi_start_annealing
from symmetry import unique_ground_states
from storage import make_run_id, save_run

"""
Were taken from refs: 
Task description : 6, 20,
https://doi.org/10.1063/1.2736681 : 24-48, 
https://doi.org/10.1186/1471-2105-8-342 : 50-100

"""

BENCHMARK_SEQUENCES = [
    {"id": "Short_Test", "length": 6, "best_known": -2,
     "sequence": "HPPHPH"},

    {"id": "Unger-Moult", "length": 20, "best_known": -9,
     "sequence": "HPHPPHHPHPPHPHHPPHPH"},

    {"id": "S1-2",  "length": 24,  "best_known": -9,
     "sequence": "H2(P2H)7H"},

    {"id": "S1-3",  "length": 25,  "best_known": -8,
     "sequence": "P2HP2(H2P4)3H2"},

    {"id": "S1-4",  "length": 36,  "best_known": -14,
     "sequence": "P3H2P2H2P5H7P2H2P4H2P2HP2"},

    {"id": "S1-5",  "length": 48,  "best_known": -23,
     "sequence": "P2H(P2H2)2P5H10P6(H2P2)2HP2H5"},

    {"id": "2D50", "length": 50, "best_known": -21,
     "sequence": "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"},

    {"id": "2D60", "length": 60, "best_known": -36,
     "sequence": "PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP"},

    {"id": "2D64", "length": 64, "best_known": -42,
     "sequence": "HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH"},

    {"id": "2D85", "length": 85, "best_known": -53,
     "sequence": "HHHHPPPPHHHHHHHHHHHHPPPPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHPPHHPPHHPPHPH"},

    {"id": "2D100a", "length": 100, "best_known": -48,
     "sequence": "PPPPPPHPHHPPPPPHHHPHHHHHPHHPPPPHHPPHHPHHHHHPHHHHHHHHHHPHHPHHHHHHHPPPPPPPPPPPHHHHHHHPPHPHHHPPPPPPHPHH"},

    {"id": "2D100b", "length": 100, "best_known": -50,
     "sequence": "PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPPHHHHHHPPHHHHHHPPPPPPPPPHPHHPHHHHHHHHHHHPPHHHPHHPHPPHPHHHPPPPPPHHH"},
 
]


def load_benchmarks():
    loaded = []
    for entry in BENCHMARK_SEQUENCES:
        seq = expand_sequence(entry["sequence"])
        if len(seq) != entry["length"]:
            raise ValueError(
                f"{entry['id']} : Not equal lengths: {len(seq)}, but expected is{entry['length']}")
        loaded.append({**entry, "sequence": seq})
    return loaded

def run_one(entry, n_starts=20, seed=0, **anneal_params):
    sequence = entry["sequence"]
    target = entry["best_known"]

    t0 = time.time()
    results = multi_start_annealing(sequence, n_starts=n_starts, seed=seed, **anneal_params)
    elapsed = time.time() - t0

    energies = []
    for e, _ in results:
        energies.append(e)
    found = min(energies)
    hits = sum(1 for e in energies if e == target)
    best_e, structures, counts = unique_ground_states(sequence, results)
    assert best_e == found, "unique_ground_states and min() disagree"

    if found < target:
        raise AssertionError(f"{entry['id']} was found {found} lower than {target}")
    return {

        "id": entry["id"],

        "length": entry["length"],

        "target": target,

        "found": found,

        "gap": found - target,

        "hit_rate": hits / n_starts,

        "mean_energy": sum(energies) / len(energies),

        "n_unique": len(structures),

        "seconds": elapsed,

        "distribution": dict(Counter(energies)),

    }

def print_table(rows):
    header = f"{'ID':<8}{'len':>5}{'target':>8}{'found':>7}{'gap':>5}{'hit%':>7}{'mean':>8}{'uniq':>6}{'sec':>7}"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(f"{r['id']:<8}{r['length']:>5}{r['target']:>8.0f}{r['found']:>7.0f}"
              f"{r['gap']:>5.0f}{100*r['hit_rate']:>6.0f}%{r['mean_energy']:>8.2f}"
              f"{r['n_unique']:>6}{r['seconds']:>7.1f}")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark: compare found energies with published values")
    parser.add_argument("--max-length", type=int, default=25,
                        help="skip sequences longer than this")
    parser.add_argument("--n-starts", type=int, default=20,
                        help="independent annealing runs per sequence")
    parser.add_argument("--seed", type=int, default=0,
                        help="same seed for all sequences (common random numbers)")
    parser.add_argument("--t-start", type=float, default=2.0)
    parser.add_argument("--t-end", type=float, default=0.05)
    parser.add_argument("--cooling", type=float, default=0.995)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--outdir", default="runs")
    parser.add_argument("--no-save", action="store_true")
    return parser.parse_args()

def main():
    args = parse_args()
    benchmarks = load_benchmarks()
    selected = [b for b in benchmarks if b["length"] <= args.max_length]

    rows = []
    for entry in selected:
        print(f"running {entry['id']} (n={entry['length']})...", flush=True)
        rows.append(run_one(entry, n_starts=args.n_starts, seed=args.seed))

    print()
    print_table(rows)

if __name__ == "__main__":
    main()


