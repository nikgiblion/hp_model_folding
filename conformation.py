import re
# sequence = "HPPHHPPH"

#Operations for random choice. 
ROTATIONAL_OPERATIONS = (
    lambda dx, dy: (-dy, dx),  #90 degrees
    lambda dx, dy: (dy, -dx), #-90 degrees
    lambda dx, dy: (-dx, -dy), #180 degrees 
)


#Expand sequence : (HP)2 -> HPHP; H2P2 -> HHPP.
def expand_sequence(compressed):
    text = compressed.strip().upper().replace(" ", "")
    pattern = re.compile(r"\(([^()]+)\)(\d*)")
    while "(" in text:
        match = pattern.search(text)
        if match is None:
            raise ValueError(f"Unpaired brackets in : {compressed}")
        inner = _expand_flat(match.group(1))
        times = int(match.group(2)) if match.group(2) else 1
        text = text[:match.start()] + inner * times + text[match.end():]

    return _expand_flat(text)

#expand without brackets : H2P3 -> HHPPP.
def _expand_flat(text):
    leftover = re.sub(r"[HP\d]", "", text)
    if leftover:
        raise ValueError(f"Not allowed symbols : {leftover} in {text}")
    result = []
    for symbol, count in re.findall(r"([HP])(\d*)", text):
        result.append(symbol * (int(count) if count else 1))
    joined = "".join(result)
    return joined

#Chekc if a sequence has only H and P residues and also normalized it.
def sequence_check(sequence):
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("The sequence is empty")

    if any(c.isdigit() or c in "()" for c in sequence):
        sequence = expand_sequence(sequence)

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

#Rotational operations to perform for tail after rotational center.
def pivot_move(coords, rng):
    n = len(coords)
    pivot_r = rng.randrange(1, n-1) #edge rotations mean nothing
    rotate_side = rng.choice(ROTATIONAL_OPERATIONS)
    px, py = coords[pivot_r]
    new_coords = coords[:pivot_r + 1] #not a problem because of tuple

    for x, y in coords[pivot_r + 1:]:
        dx, dy = x - px, y - py #from pivot_r position
        rdx, rdy = rotate_side(dx, dy)
        new_coords.append((px + rdx, py + rdy))
    return new_coords

#Generation of random initial structures for future search.
def random_conformation(n, rng, n_moves=None):
    coords = linear_conformation(n)
    if n_moves is None:
        n_moves = 2 * n
    for _ in range(n_moves):
        candidate = pivot_move(coords, rng)
        if self_avoiding_walking_check(candidate):
            coords = candidate
    return coords

