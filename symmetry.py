"""
To check if two conformations are unique or not (identical)
D_4 group is a group of quadratic symmetry

8 transformations : 4 rotations * 2 
(with reflection or without)
"""

SYMMETRIES = (
    lambda x, y: (x, y), #E 
    lambda x, y: (-y, x), #R_4 90 degrees
    lambda x, y: (-x, -y), #R_2 180 degrees
    lambda x, y: (y, -x), #270 degree rotation
    lambda x, y: (y, x), #reflection (diagonal)
    lambda x, y: (-x, y), #reflection (vertical)
    lambda x, y: (-y, -x), #reflection (2nd diagonal)
    lambda x, y: (x, -y), #reflection (horisontal)

)

#Applying symmetrical operation.
def apply_symmetry(coords, transform):
    new_coords = []
    for x, y in coords:
        new_coords.append(transform(x,y))

    return new_coords

#Normalization of translated conformation.
def normalize_translation(coords):
    min_x = min(x for x, _ in coords)
    min_y = min(y for _, y in coords)
    new_coords = []
    for x, y in coords:
        new_x, new_y = (x - min_x), (y - min_y)
        new_coords.append((new_x, new_y))

    return new_coords

#Generation of all symmetrical variants.
def canonical_form(coords, allow_reverse=True):
    variants = []
    chains = [list(coords)]
    if allow_reverse:
        chains.append(list(reversed(coords)))
    for chain in chains:
        for transform in SYMMETRIES:
            transformed = normalize_translation(apply_symmetry(chain, transform))
            variants.append(tuple(transformed))

    return min(variants)

#Generation all symmetrical variants (avoid non-palindrome).
def canonical_structure(sequence, coords, allow_reverse=True):
    variants = []
    chains = [(sequence, list(coords))]
    if allow_reverse:
        chains.append((sequence[::-1], list(reversed(coords))))
    for seq, chain in chains:
        for transform in SYMMETRIES:
            transformed = normalize_translation(apply_symmetry(chain, transform))
            variants.append((seq, tuple(transformed)))

    return min(variants)

#Chosing unique structures (count uniques and return their figures).
def unique_ground_states(sequence, results):
    best_e = min(e for e, _ in results)
    conformers = {} #coordinates that algorithm finds
    counts = {} #how many times was observed (runs to find this)
    for e, coords in results:
        if e != best_e:
            continue
        key = canonical_structure(sequence, coords)
        if key not in conformers:
            conformers[key] = coords
        counts[key] = counts.get(key, 0) + 1

    keys = list(conformers)
    return best_e, [conformers[k] for k in keys], [counts[k] for k in keys]

#Grouping by energy and chosing unique structures (count uniques and return their figures).
def unique_by_energy(sequence, results):
    by_energy = {}
    for e, coords in results:
        by_energy.setdefault(e, []).append(coords)

    levels = []
    for e in sorted(by_energy):
        conformers = {} #coordinates that algorithm finds
        counts = {} #how many times was observed (runs to find this)
        for coords in by_energy[e]:
            key = canonical_structure(sequence, coords)
            if key not in conformers:
                conformers[key] = coords
            counts[key] = counts.get(key, 0) + 1
        keys = list(conformers)
        levels.append((e, [conformers[k] for k in keys], [counts[k] for k in keys]))
    return levels

#Tests generated with Claude Opus 5
if __name__ == "__main__":
    square = [(0, 0), (1, 0), (1, 1), (0, 1)]
    assert canonical_form(square) == canonical_form([(x + 5, y + 3) for x, y in square])
    assert canonical_form(square) == canonical_form([(-y, x) for x, y in square])
    assert canonical_form(square) != canonical_form([(0,0),(1,0),(2,0),(2,1)])
    assert len({t(3, 7) for t in SYMMETRIES}) == 8
    print("symmetry: ok")