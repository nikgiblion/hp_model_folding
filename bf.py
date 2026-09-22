"""
This is a realization of brute-force approach
 for exact finidng of global minimum/a 
 in 2D lattice for protein folding.
"""

from energy import NEIGHBOR_OFFSETS

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