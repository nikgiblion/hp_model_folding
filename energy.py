#Dictionary of all contact residues and their energies.
ENERGY_CONTACT_PAIRS = {
    ("H", "H"): -1.0,
    ("P", "H"): 0.0,
    ("H", "P"): 0.0,
    ("P", "P"): 0.0,
}

#Local constant: four lattice neighbours of a site.
NEIGHBOR_OFFSETS = ((1, 0), (-1, 0), (0, 1), (0, -1))


#
# Old implementation o(n^2), saved for tests and ref.
#

#Neighbor lattice positions for two residues (r1,r2).
def lattice_neighbors(r1,r2):
    (x1,y1), (x2,y2) = r1, r2
    return abs(x1 - x2) + abs(y1 - y2) == 1


#Combining all contact pairs in one list.
def all_contacts_old(coords):
    contacts = []
    n = len(coords)
    for i in range(n):
        for j in range(i + 3, n, 2):  #update from range(i + 2, n)
            if lattice_neighbors(coords[i], coords[j]):
                contacts.append((i,j))
    return contacts

#Finding only H-H contacts.
def hh_contacts(sequence, coords):
    return [
        (i,j)
        for i, j in all_contacts(coords)
        if sequence[i] == "H" and sequence[j] == "H"
    ]

#Calculation of energy for current conformation.
def energy_old(sequence, coords):
    if len(sequence) != len(coords):
        raise ValueError(f"Different lengths: {len(sequence)} and {len(coords)}")
    etot = 0.0
    for i, j in all_contacts_old(coords):
        epair = ENERGY_CONTACT_PAIRS[(sequence[i], sequence[j])]
        etot += epair
    return etot


#
# New implementation: faster than previous - o(n) using hash map of the lattice.
#

#All non-bonded contacts as (i, j) pairs with i < j, has the same order as all_contacts_old
#Assumes coords is SAW.
def all_contacts(coords):
    occupancy = {site: index for index, site in enumerate(coords)}
    contacts = []
    for i, (x, y) in enumerate(coords):
        partners = []
        for dx, dy in NEIGHBOR_OFFSETS:
            j = occupancy.get((x + dx, y + dy))
            if j is not None and j > i + 1 and (j - i) % 2 == 1:
                partners.append(j)
        partners.sort()
        for j in partners:
            contacts.append((i, j))
    return contacts

#Energy of current conformation.
#Only H residues enter the has map, only two of four directions are scanned, every pair is visited exactly once.
def energy(sequence, coords):
    if len(sequence) != len(coords):
        raise ValueError(f"Different lengths: {len(sequence)} and {len(coords)}")
    h_sites = {site: i for i, site in enumerate(coords) if sequence[i] == "H"}
    ncontacts = 0
    for (x, y), i in h_sites.items():
        for dx, dy in ((1, 0), (0, 1)):
            j = h_sites.get((x + dx, y + dy))
            if j is not None and abs(j - i) > 1 and (j - i) % 2 == 1:
                ncontacts += 1
    return float(-ncontacts)
        

    
     
# #Tests generated with Claude Opus 5.

# if __name__ == "__main__":
#     from conformation import linear_conformation

#     assert energy("HHHH", linear_conformation(4)) == 0
#     square = [(0, 0), (1, 0), (1, 1), (0, 1)]
#     assert energy("HHHH", square) == -1
#     assert energy("HPPH", square) == -1
#     assert energy("PHHP", square) == 0
#     assert all(
#     value == 0.0
#     for pair, value in ENERGY_CONTACT_PAIRS.items()
#     if pair != ("H", "H")), "New energy() assumes only H-H contributions. If table is changed, revise!"
#     assert energy_old("HHHH", square) == energy("HHHH", square)
#     print("energy: ok")