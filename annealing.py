import random
import math

from conformation import self_avoiding_walking_check
from energy import energy

#Operations for random choice. 
ROTATIONAL_OPERATIONS = (
    lambda dx, dy: (-dy, dx),  #90 degrees
    lambda dx, dy: (dy, -dx), #-90 degrees
    lambda dx, dy: (-dx, -dy), #180 degrees 
)

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



#Metropolis Criterium to accept operations to go out of local minimum.
# dE <= 0 (always take it), else with probability p = exp(-dE/T)
def metropolis_crit(delta_e, temperature, rng):
    if delta_e <= 0:
        return True
    return rng.random() < math.exp(-delta_e / temperature) #Boltzmann distr

#Main function of annealing process.
def annealing_process(sequence, start_coords, t_start=2.0, t_end=0.05, cooling=0.995, steps_per_temp=200, seed=None):
    rng = random.Random(seed)
    current = list(start_coords)
    current_e = energy(sequence, current)
    best, best_e = list(current), current_e
    temperature = t_start
    history = []  #not to allocate too much, clear each activation
    while temperature > t_end:
        for _ in range(steps_per_temp): #only repeats
            candidate = pivot_move(current, rng)
            if not self_avoiding_walking_check(candidate):
                continue
            candidate_e = energy(sequence, candidate)
            diff_e = candidate_e - current_e
            if metropolis_crit(diff_e, temperature, rng):
                current, current_e = candidate, candidate_e
                if current_e < best_e:
                    best, best_e = list(current), current_e
        history.append((temperature, current_e, best_e))
        temperature *= cooling
    return best, best_e, history

