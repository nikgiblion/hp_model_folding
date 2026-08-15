sequence = "HPPHHPPH"

#Chekc if a sequence has only H and P residues and also normalized it.
def sequence_check(sequence):
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("The sequence is empty")
    wrong_symbols = set(sequence) - {"H", "P"}
    if wrong_symbols:
        raise ValueError(f"Forbidden symbols: {sorted(wrong_symbols)}")

    return sequence


#making linear conformation for sequence as an initial state on the lattice.
def linear_conformation(n):
    coords = []
    for i in range(n):
        point = (i,0)
        coords.append(point)
    return coords

#True if each pair of neighbors in sequence are neighbors in the lattice.
def connectivity_check(coords):
    for (x1,y1), (x2,y2) in zip(coords, coords[1:]):
        if abs(x1-x2) + abs(y1-y2) != 1:
            return False
    return True

#Check if two residues are in the same coordinate
def self_avoiding_walking_check(coords):
    return len(set(coords)) == len(coords)


#Two important properties to check
def final_validation_sequence(coords):
    return connectivity_check(coords) and self_avoiding_walking_check(coords)

