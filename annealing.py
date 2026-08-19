import random
import math

from energy import energy

# from conformation import self_avoiding_walking_check
# from conformation import pivot_move
# from conformation import random_conformation

from conformation import (
    self_avoiding_walking_check,
    random_move,

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
def annealing_process(sequence, start_coords, t_start=2.0, t_end=0.05, cooling=0.995, steps_per_temp=None, seed=None):
    if steps_per_temp is None:
        steps_per_temp = 10 * len(sequence) #for n=20 -> 200, for 48 -> 480
    rng = random.Random(seed)
    current = list(start_coords)
    current_e = energy(sequence, current)
    best, best_e = list(current), current_e
    temperature = t_start
    history = []  #not to allocate too much, clear each activation

    attempted = 0
    rejected_no_move = 0
    rejected_sa = 0
    rejected_metropolis = 0
    while temperature > t_end:
        for _ in range(steps_per_temp): #only repeats
            attempted += 1 
            candidate = random_move(current, rng)
            if candidate is None:
                rejected_no_move += 1
                continue

            if not self_avoiding_walking_check(candidate):
                rejected_sa += 1
                continue
            
            candidate_e = energy(sequence, candidate)
            diff_e = candidate_e - current_e
            if metropolis_crit(diff_e, temperature, rng):
                current, current_e = candidate, candidate_e
                if current_e < best_e:
                    best, best_e = list(current), current_e
            else:
                rejected_metropolis += 1
        history.append((temperature, current_e, best_e))
        temperature *= cooling

    accepted = attempted - rejected_no_move - rejected_metropolis - rejected_sa
    stats_moves = {
        "attempted": attempted,
        "accepted": accepted,
        "sa_reject_rate": rejected_sa / attempted, 
        "no_move_rate": rejected_no_move / attempted,
        "metropolis_reject_rate": rejected_metropolis / attempted,
        "accept_rate": accepted / attempted,

    }
    return best, best_e, history, stats_moves

#To do multiple starts with different initial states (linear + nstarts - 1 structures)
def multi_start_annealing(sequence, nstarts=10, seed=0, t_start=2.0, t_end=0.05, cooling=0.995, steps_per_temp=None):    
    rng = random.Random(seed)
    n = len(sequence)
    results = []
    all_stats = []

    for i in range(nstarts):
        if i == 0:
            start = linear_conformation(n)
        else:
            start = random_conformation(n,rng)

        run_seed = rng.randrange(2**32)
        best, best_e, history, stats_moves = annealing_process(sequence, start, seed=run_seed, t_start=t_start, t_end=t_end, cooling=cooling, steps_per_temp=steps_per_temp)
        results.append((best_e, best))
        all_stats.append(stats_moves)

    mean_stats = {                                   
        key: sum(s[key] for s in all_stats) / len(all_stats)
        for key in all_stats[0]
    }

    return results, mean_stats

