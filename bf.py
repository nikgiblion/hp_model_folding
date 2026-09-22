"""
This is a realization of brute-force approach
 for exact finidng of global minimum/a 
 in 2D lattice for protein folding.
"""

from energy import NEIGHBOR_OFFSETS

def iter_saws(n):
    if n < 1:
        raise ValueError(f"Impossible sequence")

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

