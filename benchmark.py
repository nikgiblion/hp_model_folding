import argparse
import time
from datetime import datetime
from collections import Counter

from conformation import expand_sequence
from annealing import multi_start_annealing
from symmetry import unique_ground_states
from storage import save_report

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

def run_one(entry, nstarts=20, seed=0, **anneal_params):
    sequence = entry["sequence"]
    target = entry["best_known"]

    t0 = time.time()
    results, move_stats = multi_start_annealing(sequence, nstarts=nstarts, seed=seed, **anneal_params)
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

        "hit_rate": hits / nstarts,

        "mean_energy": sum(energies) / len(energies),

        "n_unique": len(structures),

        "seconds": elapsed,

        "distribution": dict(Counter(energies)),

        "sa_reject": move_stats["sa_reject_rate"],

        "accept": move_stats["accept_rate"],

    }

#Forming main information about benchmark
def format_table(rows):
    header = (f"{'ID':<8}{'len':>5}{'target':>8}{'found':>7}{'gap':>5}"
              f"{'hit%':>7}{'mean':>8}{'uniq':>6}{'saRej':>7}{'acc':>6}{'sec':>7}")
    lines = [header, "-" * len(header)]
    for r in rows:
        lines.append(
            f"{r['id']:<8}{r['length']:>5}{r['target']:>8.0f}{r['found']:>7.0f}"
            f"{r['gap']:>5.0f}{100*r['hit_rate']:>6.0f}%{r['mean_energy']:>8.2f}"
            f"{r['n_unique']:>6}{100*r['sa_reject']:>6.0f}%"
            f"{100*r['accept']:>5.0f}%{r['seconds']:>7.1f}")
    return "\n".join(lines)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark: compare found energies with published values")
    parser.add_argument("--max-length", type=int, default=25,
                        help="skip sequences longer than this")
    parser.add_argument("--nstarts", type=int, default=20,
                        help="independent annealing runs per sequence")
    parser.add_argument("--seed", type=int, default=0,
                        help="same seed for all sequences (common random numbers)")
    parser.add_argument("--t-start", type=float, default=2.0,
                        help="initial temperature")
    parser.add_argument("--t-end", type=float, default=0.05,
                        help="final temperature")
    parser.add_argument("--cooling", type=float, default=0.995,
                        help="geometric cooling factor per temperature step")
    parser.add_argument("--steps_per_temp", type=int, default=None,
                        help="moves per temperature (default: 10 x chain length)")
    parser.add_argument("--outdir", default="runs",
                        help="directory for saved reports")
    parser.add_argument("--no-save", action="store_true",
                        help="do not write a report file")
    return parser.parse_args()

def main():
    args = parse_args()
    benchmarks = load_benchmarks()
    selected = [b for b in benchmarks if b["length"] <= args.max_length]

    rows = []
    print("")
    print("Running Benchmark dataset")
    for entry in selected:
        print(f"running {entry['id']} (n={entry['length']})...", flush=True)
        rows.append(run_one(
            entry,
            nstarts=args.nstarts,
            seed=args.seed,
            t_start=args.t_start,
            t_end=args.t_end,
            cooling=args.cooling,
            steps_per_temp=args.steps_per_temp,
        ))

    table_text = format_table(rows)
    print()
    print(table_text)

    if not args.no_save:
        header = (f"benchmark  max_length={args.max_length} nstarts={args.nstarts} "
                  f"seed={args.seed} cooling={args.cooling} "
                  f"t_start={args.t_start} t_end={args.t_end}")
        run_id = "benchmark_" + datetime.now().strftime("%Y%m%d_%H%M%S")
        path = save_report(header + "\n\n" + table_text, run_id, outdir=args.outdir)
        print(f"Saved to : {path}")

if __name__ == "__main__":
    main()


