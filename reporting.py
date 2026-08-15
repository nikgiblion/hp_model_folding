from energy import all_contacts, hh_contacts, energy
from visualization import ascii_vizualization

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