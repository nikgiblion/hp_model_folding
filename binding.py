"""
This module is extension to the task and is about finding binding sites
for ligands. Two types of ligand : H and P, based on it contact energy is different.
"""

from collections import deque
from conformation import DIRECTIONS

LIGAND_ENERGY = {
    ("H", "H"): -1.0, #for H-Ligand
    ("H", "P"): 0.0,
    ("P", "H"): 0.0,
    ("P", "P"): -1.0, #for P-Ligand
}

#on the edge cavities (3 member surrounding)
def edge_cavities(coords, margin=2):
    occupied = set(coords)
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    min_x, max_x = min(xs) - margin, max(xs) + margin
    min_y, max_y = min(ys) - margin, max(ys) + margin

    start = (min_x, min_y)
    seen = {start}
    queue = deque([start])

    while queue:
        x, y = queue.popleft()
        for dx, dy in DIRECTIONS:
            nx, ny = x + dx, y + dy
            if not (min_x <= nx <= max_x and min_y <= ny <= max_y):
                continue
            if (nx, ny) in occupied or (nx, ny) in seen:
                continue
            seen.add((nx, ny))
            queue.append((nx, ny))

    return seen

#inside cavities or holes (4 or more member surrounding)
def inside_cavities(coords):
    occupied = set(coords)
    outside = edge_cavities(coords)
    holes = set()
    for x, y in coords:
        for dx, dy in DIRECTIONS:
            site = (x + dx, y + dy)
            if site not in occupied and site not in outside:
                holes.add(site)
    return holes

def site_profile(site, coords, sequence):
    position_to_index = {point: i for i, point in enumerate(coords)}
    x, y = site
    types = []
    for dx, dy in DIRECTIONS:
        neighbor = (x + dx, y + dy)
        if neighbor in position_to_index:
            types.append(sequence[position_to_index[neighbor]])
    return tuple(sorted(types))

def check_cavities(sequence, coords, ligand="H", min_neighbors=2):
    sites = []
    for site in inside_cavities(coords):
        profile = site_profile(site, coords, sequence)
        if len(profile) < min_neighbors:
            continue
        e_bind = sum(LIGAND_ENERGY[(ligand, t)] for t in profile)
        sites.append({
            "site": site,
            "profile": profile,
            "n_neighbors": len(profile),
            "energy": e_bind,
        })
    return sorted(sites, key=lambda s: s["energy"])

#Tests generated with Claude Opus 5
if __name__ == "__main__":
    from conformation import linear_conformation

    ring = [(0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1)]
    assert inside_cavities(ring) == {(1, 1)}
    assert inside_cavities(linear_conformation(10)) == set()

    result = check_cavities("HHHHHHHH", ring, ligand="H")
    assert result[0]["energy"] == -4.0
    assert check_cavities("PPPPPPPP", ring, ligand="H")[0]["energy"] == 0.0
    print("binding: ok")
