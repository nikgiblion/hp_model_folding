import random
import math

from energy import energy

# from conformation import self_avoiding_walking_check
# from conformation import pivot_move
# from conformation import random_conformation

from conformation import (
    self_avoiding_walking_check,
    pivot_move,
    random_conformation,
    linear_conformation,
)

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

#To do multiple starts with different initial states (linear + n_starts - 1 structures)
def multi_start_annealing(sequence, n_starts=10, seed=0, t_start=2.0, t_end=0.05, cooling=0.995, steps_per_temp=200):    
    rng = random.Random(seed)
    n = len(sequence)
    results = []

    for i in range(n_starts):
        if i == 0:
            start = linear_conformation(n)
        else:
            start = random_conformation(n,rng)

        run_seed = rng.randrange(2**32)
        best, best_e, history = annealing_process(sequence, start, seed=run_seed, t_start=t_start, t_end=t_end, cooling=cooling, steps_per_temp=steps_per_temp)
        results.append((best_e, best))

    return results

