from datetime import datetime
from collections import Counter
from energy import (energy, all_contacts, hh_contacts)
from symmetry import unique_ground_states
import json
import os

#ID creation
def make_run_id(sequence, seed):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{stamp}_n{len(sequence)}_s{seed}"

#Storage of data from run
def build_results_dict(sequence, results, params, run_id):
    best_e, structures, counts = unique_ground_states(sequence, results)
    distribution = Counter (e for e, _ in results)

    ground_states = []
    for i, (coords, count) in enumerate(zip(structures, counts), start=1):
        ground_states.append({
            "index": i,
            "energy": best_e,
            "found_times": count,
            "coords": [list(p) for p in coords],
            "HH_contacts": [list(c) for c in hh_contacts(sequence, coords)],
        })

    return {
        "run_id": run_id,
        "sequence": sequence,
        "length": len(sequence),
        "composition": {"H": sequence.count("H"), "P": sequence.count("P")},
        "params": params,
        "n_runs": len(results),
        "energy_distribution": {str(int(e)): c for e, c in sorted(distribution.items())},
        "all_energies": [int(e) for e, _ in results],
        "best_energy": int(best_e),
        "n_unique_ground_states": len(structures),
        "ground_states": ground_states,
        } 

#Writing JSON and txt files
def save_run(result_dict, report_text, outdir="runs"):
    run_dir = os.path.join(outdir, result_dict["run_id"])
    os.makedirs(run_dir, exist_ok=True)

    json_path = os.path.join(run_dir, "result.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result_dict, f, indent=2, ensure_ascii=False)

    report_path = os.path.join(run_dir, "report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_text)

    return run_dir