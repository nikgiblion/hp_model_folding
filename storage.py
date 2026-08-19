from datetime import datetime
from pathlib import Path

#ID creation
def make_run_id(sequence, seed):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{stamp}_n{len(sequence)}_s{seed}"

#Writing txt file
def save_report(report_text, run_id, outdir="runs"):
    run_dir = Path(outdir)
    run_dir.mkdir(parents=True, exist_ok=True)

    path = run_dir / f"{run_id}.txt"
    path.write_text(report_text, encoding="utf-8")
    return path