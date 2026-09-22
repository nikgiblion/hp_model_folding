"""
This is a realization of brute-force approach
 for exact finidng of global minimum/a 
 in 2D lattice for protein folding.
"""

from energy import NEIGHBOR_OFFSETS
from symmetry import canonical_structure

def iter_saws(n):
    if n < 1:
        raise ValueError(f"Impossible sequence, must be at least 1")

    coords = [(0,0)]
    occupied = {(0,0)}

    def grow():
        if len(coords) == n:
            yield tuple(coords) #copy: the list keeps updated after yield
            return

        x, y = coords[-1]
        for dx, dy in NEIGHBOR_OFFSETS:
            site = (x + dx, y + dy)
            if site in occupied: #if self-intersection -> cut this branch
                continue
            coords.append(site)
            occupied.add(site)
            yield from grow()
            occupied.remove(site)
            coords.pop()

    yield from grow()

#Directions only afte 2nd residue: straight or up (+x or +y).
BEFORE_STEP = ((1,0), (0, 1))

def iter_symmetry_saws(n):
    if n < 1:
        raise ValueError(f"mpossible sequence, must be at least 1")
    if n == 1:
        yield ((0, 0),)
        return

    coords = [(0, 0), (1, 0)]
    occupied = {(0, 0), (1, 0)}

    def grow(turned):
        if len(coords) == n:
            yield tuple(coords)
            return
        x, y = coords[-1]
        for dx, dy in (NEIGHBOR_OFFSETS if turned else BEFORE_STEP):
            site = (x + dx, y + dy)
            if site in occupied:
                continue
            coords.append(site)
            occupied.add(site)
            yield from grow(turned or dy == 1) #1st +y step is 1st step
            occupied.remove(site)
            coords.pop()

    yield from grow(turned=False)


def exact_ground_states(sequence, prune=True):
    n = len(sequence)
    if n == 0:
        raise ValueError("The sequence is empty")
    if n <= 3:  #no contacts
        return 0.0, list(iter_symmetry_saws(n))

    coords = [(0, 0), (1, 0)]
    occupied = {(0, 0), (1, 0)}
    h_index = {site: i for i, site in enumerate(coords) if sequence[i] == "H"}
    best = -1  #largest number of H-H contacts found
    optimal = []
    #upper bound on the contacts that residues k...n-1 can still add
    free_valence = [(0 if aa != "H" else 3 if i in (0, n - 1) else 2) for i, aa in enumerate(sequence)]
    bound = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        bound[i] = bound[i + 1] + free_valence[i]

    def grow(turned, contacts):
        nonlocal best, optimal
        k = len(coords) #index of residue to place now
        if prune and contacts + bound[k] < best:
            return  #cut the subtree if not better than best     
        if k == n:
            if contacts > best:
                best, optimal = contacts, [tuple(coords)]
            elif contacts == best:
                optimal.append(tuple(coords))
            return
        x, y = coords[-1]
        is_h = sequence[k] == "H"
        for dx, dy in (NEIGHBOR_OFFSETS if turned else BEFORE_STEP):
            site = (x+ dx, y + dy)
            if site in occupied:
                continue
            gained = 0
            if is_h:
                for ox, oy in NEIGHBOR_OFFSETS:
                    j = h_index.get((site[0] + ox, site[1] + oy))
                    if j is not None and j != k - 1:  #the chain neighbor is a bond
                        gained += 1
            coords.append(site)
            occupied.add(site)
            if is_h:
                h_index[site] = k
            grow(turned or dy == 1, contacts + gained)
            if is_h:
                del h_index[site]
            occupied.remove(site)
            coords.pop()

    grow(False, 0)
    return float(-best), optimal

#Ground states up to lattice symmetry and chain reversal (palindrome check)
def exact_distinct_ground_states(sequence, prune=True):
    e, optimal = exact_ground_states(sequence, prune=prune)
    representatives = {}
    for coords in optimal:
        representatives.setdefault(canonical_structure(sequence, coords), coords)
    return e, list(representatives.values())