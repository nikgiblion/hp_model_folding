import argparse
from conformation import sequence_check, linear_conformation
from annealing import multi_start_annealing
from reporting import report, multi_report
from storage import make_run_id, save_report

DEFAULT_SEQUENCE = "HPHPPHHPHPPHPHHPPHPH" #know benchmark sequence with E = -9

def parse_args():
    parser = argparse.ArgumentParser(description="HP model folding with annealing algorithm for 2D cases (ver.0.1)")
    parser.add_argument("--sequence", default=DEFAULT_SEQUENCE)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--t-start", type=float, default=2.0)
    parser.add_argument("--t-end", type=float, default=0.05)
    parser.add_argument("--cooling", type=float, default=0.995)
    parser.add_argument("--steps_per_temp", type=int, default=None)  
    parser.add_argument("--nstarts", type=int, default=1) #number of independent annealing starts
    parser.add_argument("--outdir", default="runs")
    parser.add_argument("--no-save", action="store_true")
    parser.add_argument("--max-shown", type=int, default=5)
    parser.add_argument("--binding", action="store_true")
    parser.add_argument("--ligand", default="H", choices=["H", "P"])
    return parser.parse_args()

def main():
    args = parse_args()
    sequence = sequence_check(args.sequence)
    start = linear_conformation(len(sequence))

    print(report(sequence, start))

    results, move_stats = multi_start_annealing(
        sequence, nstarts=args.nstarts,
        seed=args.seed, t_start=args.t_start,
        t_end=args.t_end, cooling=args.cooling,
        steps_per_temp=args.steps_per_temp,
    )

    print(multi_report(sequence, results, max_shown=args.max_shown, binding=args.binding, ligand=args.ligand))
    print(f"Accepted moves: {100*move_stats['accept_rate']:.0f}%, "
          f"rejected by self-avoidance: {100*move_stats['sa_reject_rate']:.0f}%")

    if not args.no_save:
        run_id = make_run_id(sequence, args.seed)
        params = dict(vars(args))
        keys = ("seed", "nstarts", "t_start", "t_end", "cooling", "steps_per_temp")
        header = "params: " + " ".join(f"{k}={params[k]}" for k in keys)
        full_report = header + "\n\n" + multi_report(
            sequence, results, max_shown=None, per_level=5,
            binding=args.binding, ligand=args.ligand)
        path = save_report(full_report, run_id, outdir=args.outdir)
        print(f"Saved to : {path}")

if __name__ == "__main__":
    main()