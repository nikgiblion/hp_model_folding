from energy import all_contacts, hh_contacts, energy
from visualization import ascii_vizualization
from symmetry import unique_ground_states

from collections import Counter

#Preparing 1st version of report, only main information.
def report(sequence, coords):
    all_pairs = all_contacts(coords)
    hh_pairs = hh_contacts(sequence, coords)
    lines = [
        # "Total information after global optimization run",
        "",
        "Conformation",
        ascii_vizualization(sequence, coords),
        "",
        f"Sequence : {sequence} (length : {len(sequence)})",
        f"Composition : H = {sequence.count('H')}, P = {sequence.count('P')}",
        f"Energy : {energy(sequence, coords):.0f}",
        f"All contacts : {len(all_pairs)}",
        f"Only H-H contacts : {len(hh_pairs)}",
        f"Contact pairs : {hh_pairs if hh_pairs else 'none'}", #indices for contact pairs in H-H contacts
        "",

    ]
    return "\n".join(lines)

def multi_report(sequence, results, max_shown=5):
    best_e, structures, counts = unique_ground_states(sequence, results)

    lines = []
    lines.append("=" * 50)
    lines.append(f"Runs : {len(results)}")
    lines.append("")
    energies = sorted(e for e, _ in results)
    distribution = Counter(energies)
    lines.append("Energy distribution : ")
    for e in sorted(distribution):
        bar = "#" * distribution[e]
        lines. append(f" E = {e:5.0f}  {bar}  ({distribution[e]})")

    lines.append("")
    lines.append(f"Energy minimum       : {best_e:.0f}")
    lines.append(f"Number of conformers : {len(structures)}")
    lines.append("")

    for i, (coords, count) in enumerate(zip(structures, counts), start=1):
        if i > max_shown:
            lines.append(f"... and also {len(structures) - max_shown} structures")
            break
        lines.append(f"--- structure {i} (was found {count} times) ---")
        lines.append(report(sequence,coords))
        lines.append("")
    return "\n".join(lines)
