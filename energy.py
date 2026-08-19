#From conformation import linear_conformation

#Dictionary of all contact residues and their energies.
ENERGY_CONTACT_PAIRS = {
    ("H", "H"): -1.0,
    ("P", "H"): 0.0,
    ("H", "P"): 0.0,
    ("P", "P"): 0.0,
}

#Neighbor lattice positions for two residues (r1,r2).
def lattice_neighbors(r1,r2):
    (x1,y1), (x2,y2) = r1, r2
    return abs(x1 - x2) + abs(y1 - y2) == 1


#Combining all contact pairs in one list.
def all_contacts(coords):
    contacts = []
    n = len(coords)
    for i in range(n):
        for j in range(i + 3, n, 2):  #update from range(i + 2, n)
            if lattice_neighbors(coords[i], coords[j]):
                contacts.append((i,j))
    return contacts

#Finding only H-H contacts
def hh_contacts(sequence, coords):
    return [
        (i,j)
        for i, j in all_contacts(coords)
        if sequence[i] == "H" and sequence[j] == "H"
    ]

#Calculation of energy for current conformation.
def energy(sequence, coords):
    if len(sequence) != len(coords):
        raise ValueError(f"Different lengths: {len(sequence)} and {len(coords)}")
    etot = 0.0
    for i, j in all_contacts(coords):
        epair = ENERGY_CONTACT_PAIRS[(sequence[i], sequence[j])]
        etot += epair
    return etot

#Tests generated with Claude Opus 5

if __name__ == "__main__":
    from conformation import linear_conformation

    assert energy("HHHH", linear_conformation(4)) == 0
    square = [(0, 0), (1, 0), (1, 1), (0, 1)]
    assert energy("HHHH", square) == -1
    assert energy("HPPH", square) == -1
    assert energy("PHHP", square) == 0
    print("energy: ok")