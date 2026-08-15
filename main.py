import argparse
from conformation import sequence_check, linear_conformation
from annealing import (annealing_process, multi_start_annealing)
from reporting import (report, multi_report)
from energy import energy

DEFAULT_SEQUENCE = "HPHPPHHPHPPHPHHPPHPH" #know benchmark sequence with E = -9

def parse_args():
    parser = argparse.ArgumentParser(description="HP model folding with annealing algorithm for 2D cases (ver.0.1)")
    parser.add_argument("--sequence", default=DEFAULT_SEQUENCE)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--t-start", type=float, default=2.0)
    parser.add_argument("--t-end", type=float, default=0.05)
    parser.add_argument("--cooling", type=float, default=0.995)
    parser.add_argument("--steps", type=int, default=200)  
    parser.add_argument("--n_starts", type=int, default=1) #number of independent annealing starts
    return parser.parse_args()

def main():
    args = parse_args()
    sequence = sequence_check(args.sequence)
    start = linear_conformation(len(sequence))

    print(report(sequence, start))

    if args.n_starts == 1:
        best, best_e, history = annealing_process(
        sequence, start,
        t_start=args.t_start, t_end=args.t_end,
        cooling=args.cooling, steps_per_temp=args.steps,
        seed=args.seed,
        )

        print(report(sequence, best))
    else:
        results = multi_start_annealing(
            sequence, n_starts=args.n_starts,
            seed=args.seed, t_start=args.t_start,
            t_end=args.t_end, cooling=args.cooling,
            steps_per_temp=args.steps,
            )
        print(multi_report(sequence, results))

        # print(f"Energy : {energy(sequence, start):.0f} -> {best_e:.0f}")
        # print(f"Temperature of steps : {len(history)}")

if __name__ == "__main__":
    main()