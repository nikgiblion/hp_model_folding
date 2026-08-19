from energy import all_contacts, hh_contacts, energy
from visualization import ascii_vizualization
from symmetry import unique_ground_states, unique_by_energy
from binding import check_cavities
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

#Reporting about binding sites
def binding_report(sequence, coords, ligand="H"):
    sites = check_cavities(sequence, coords, ligand=ligand)
    if not sites:
        return "No inside cavities found."

    lines = [f"Cavities found : {len(sites)} (Ligand type : {ligand})", ""]
    for i, s in enumerate(sites, start=1):
        profile = "".join(s["profile"])
        lines.append(f"{i}. site {s['site']}, profile {profile} "
                     f"neighbors {s['n_neighbors']},  binding energy = {s['energy']:.0f}")
    best = sites[0]
    lines.append("")
    lines.append(f"Best binding site : {best['site']}, E_bind = {best['energy']:.0f}")
    lines.append("")
    lines.append("Best site only")
    lines.append(ascii_vizualization(sequence, coords, ligand_sites=[best["site"]]))
    lines.append("")
    lines.append("All binding sites")
    if len(sites) > 1:
        lines.append("")
        lines.append("All binding sites")
        all_sites = [s["site"] for s in sites]
        lines.append(ascii_vizualization(sequence, coords,
                                         ligand_sites=all_sites))
    return "\n".join(lines)

#Old one without json and txt
def multi_report(sequence, results, max_shown=None, per_level=None, binding=False, ligand="H"):
    best_e, structures, counts = unique_ground_states(sequence, results)

    lines = []
    lines.append("=" * 50)
    lines.append(f"Runs : {len(results)}")
    lines.append("")
    distribution = Counter(e for e, _ in results)
    lines.append("Energy distribution : ")
    max_count = max(distribution.values())
    for e in sorted(distribution):
        bar_len = round(30 * distribution[e] / max_count)
        lines.append(f" E = {e:5.0f}  {'#' * bar_len}  ({distribution[e]})")

    lines.append("")
    lines.append(f"Energy minimum       : {best_e:.0f}")
    lines.append(f"Number of conformers : {len(structures)}")
    lines.append("")

    for i, (coords, count) in enumerate(zip(structures, counts), start=1):
        if max_shown is not None and i > max_shown:
            lines.append(f"... and also {len(structures) - max_shown} structures")
            break
        lines.append(f"--- structure {i} (was found {count} times) ---")
        lines.append(report(sequence, coords))
        lines.append("")

    #This is for different energy levels
    if per_level is not None:
        lines.append("=" * 50)
        lines.append("All energy levels")
        lines.append("")
        for e, level_structures, level_counts in unique_by_energy(sequence, results):
            lines.append(f"### E = {e:.0f} : {len(level_structures)} unique, "
                         f"{sum(level_counts)} runs")
            lines.append("")
            for i, (coords, count) in enumerate(zip(level_structures, level_counts), start=1):
                if per_level > 0 and i > per_level:
                    lines.append(f"    ... and {len(level_structures) - per_level} more")
                    break
                lines.append(f"--- E={e:.0f}, structure {i} (found {count} times) ---")
                lines.append(report(sequence, coords))
                lines.append("")

    # Ligand binding analysis
    if binding:
        lines.append("=" * 50)
        lines.append("Ligand binding analysis")
        lines.append("")
        for i, coords in enumerate(structures, start=1):
            lines.append(f"--- ground state {i} ---")
            lines.append(binding_report(sequence, coords, ligand=ligand))
            lines.append("")

    return "\n".join(lines)

