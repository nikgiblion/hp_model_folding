## Results

The method was validated against reference sequences from the task description and from the literature (see README.md). For HPPHPH (6 residues) the search reaches E = −2 in every run, matching the stated optimum. For the Unger–Moult 20-mer the published ground-state energy of −9 is reached in 32% of 100 independent runs. Every run that reaches it converges to the same conformation up to the symmetries of the square lattice and chain reversal, so the ground-state appears to be non-degenerate for this sequence.

Benchmark against published values (100 independent runs per sequence, one shared seed. `Unique GS` = unique ground-states `saRej` = moves rejected for self-intersection, `acc` = moves accepted):

| ID | Length | Published | Found | Gap | hit % | Mean E | Unique GS | saRej | acc |
|----|--------|-----------|-------|-----|-------|--------|-----------|-------|-----|
| HPPHPH | 6 | −2 | −2 | 0 | 100% | −2.00 | 1 | 44% | 18% |
| Unger–Moult | 20 | −9 | −9 | 0 | 32% | −8.28 | 1 | 59% | 13% |
| S1-2 | 24 | −9 | −9 | 0 | 38% | −8.37 | 11 | 61% | 15% |
| S1-3 | 25 | −8 | −8 | 0 | 22% | −7.05 | 11 | 60% | 20% |
| S1-4 | 36 | −14 | −14 | 0 | 2% | −11.76 | 2 | 63% | 21% |
| S1-5 | 48 | −23 | −21 | 2 | 0% | −18.39 | 4 | 68% | 14% |

Up to 25 residues the published minimum is found reliably. At 36 residues it is reached in roughly 2% of runs, and at 48 the search falls not so far from ground-state. Published energies for this length come from more sophisticated algorithms (replica exchange and chain-growth Monte Carlo) with compute budgets larger than a single annealing run.

Ground-state degeneracy differs sharply between sequences. The 20-mer converges to a single structure across all 32 successful runs, while S1-2 and S1-3 yield eleven distinct ones at the same minimum energy. The two structures reported for S1-4 come from only two successful runs and carry no statistical weight.

**Extension**

No internal cavities were found in the ground-states of the 20-mer. A chain that short cannot enclose a free lattice site without losing more contacts than the enclosure gains. Cavities do appear at 48 residues, which is consistent with binding pockets being a property of larger chains.

## Key decisions

**Absolute coordinates.** A conformation is a list of `(x, y)` tuples, one per residue, with connectivity implicit in list order. Energy evaluation and the self-avoidance check — the two operations in the inner loop — are then direct on coordinates: a contact is a Manhattan distance of 1, self-avoidance is `len(set(coords)) == len(coords)`. A turn or direction encoding would make connectivity structural, but every energy evaluation would first have to unfold the chain into coordinates. Connectivity is instead guaranteed by construction of the move set and verified by assertion.

**Move set.** A mixture of 80% pivot and 20% local moves (tail move and corner flip). Pivot rotates the whole tail beyond a chosen residue; because lattice rotations map adjacent sites to adjacent sites, connectivity is preserved by construction and only self-avoidance can break. Pivot is ergodic for self-avoiding walks [2]. Local moves are not ergodic on their own, so the mixture is what makes the search both complete in principle and effective in practice.

**Energy**: Using special table with all contact pairs to make it more flexible for changes within possible extensions of HP-model.

**Lattice bipartiteness**: Square-lattice sites split into two sublattices by the parity of `x + y`, and each chain step flips the parity. Two residues can therefore be in contact only if `|i - j|` is odd, so the contact search skips even index differences entirely.

**Structure comparison** Counting distinct ground-states requires deciding when two conformations are the same. Conformations are canonicalised over the eight symmetries of the square, translation, and chain reversal, taking the lexicographic minimum over all resulting variants. Mirror images are treated as identical, since in 2D a reflection is reachable by rotation through the third dimension and the model has no chirality. Reversal is
applied to the sequence and the coordinates together, which matters for non-palindromic sequences.

**Search strategy** Independent restarts from randomized initial conformations, with each run seeds drawn from a single RNG, so that a whole experiment is reproducible from one seed while individual runs remain independently replayable for debugging.

## What I measured

**Bipartiteness effect** After restricting the contact search to odd index differences, re-running the benchmark with the same seed reproduced every hit rate and mean energy to the last digit while halving runtime (193 s goes to 130 s on the 20-mer). This is the check I would want for any optimisation claimed to preserve results.

**The optimal probability of pivot moves** I tried sweeping the pivot probability over 0.4 / 0.6 / 0.8 / 1.0 (100 runs, five sequences, shared seed), mean energy improved towards 0.8 and was worse at 1.0 (pure pivot) on all five sequences. At n = 36 residues the mixture closed the gap to the published value.

**Local moves role** I expected local moves to help by raising acceptance and a pivot on a densely packed chain usually makes the tail intersect the rest. Measuring the rejection breakdown showed the self-intersection rate essentially unchanged between regimes (59–60% with the mixture, 61–63% with pure pivot). The remaining explanation is step size: a local move changes the energy by 0–1 while a pivot changes it by several, so at low temperature the Metropolis criterion admits local moves long after it has frozen out pivots. This would be worth confirming by comparing the Metropolis rejection rates between the two regimes.

**Saturation of restarts** If a number of independent starts is increased from 20 to 100 for 20-mer seqeunce the hit rate goes from 25% to 32%. More restarts only estimate it more precisely. Improving the result therefore requires changing the individual run: move set, schedule, or algorithm.

**Reproducibility** Same seed gives the same distributions of energies across restarts.

## Limitations

**Nothing is proven exhaustively** No exhaustive enumeration was implemented, so even for the 6-mer the agreement with E = −2 is a search result matching a published value rather than a proof. Enumeration is entirely feasible at that length (284 self-avoiding walks) and is the first gap I would close.

**Degeneracy issue.** "One ground-state structure for the 20-mer" means every run reaching −9 found the same structure up to symmetry. Structures with specific conformation may simply never be sampled.

**The method does not scale beyond 36 residues.** Reaching published values at these lengths would require a different algorithm, such as basin hopping or replica exchange Monte Carlo[4].

**Ergodicity not verified** Pivot moves are ergodic for self-avoiding walks in general and the mixture inherits this because pivot is always available. My implementation has constration for pivots to rotations of the chain tail about interior residues.

**Ligand binding is an extension** The binding interaction table differs from the residue-residue scheme. Binding is evaluated on a rigid, already-folded conformation, so the model doesn't capture the way of obtaining this cavity.

**Statistical resolution.** Hit rates come from 100 runs, giving a standard error of 4–5%. Differences below roughly 10 points are not differ. Comparisons in this report therefore rest on mean energy and on consistency across several sequences.

## Discussion questions

1. Why must chain-adjacent H–H pairs be excluded from the contact count? What goes wrong — silently — if you forget the self-avoidance constraint or the chain-adjacency exclusion?

These H-H residues are covalently bonded and exist in any conformer. Therefore, it is not so important to count the same sequence when comparing different conformations. They can be excluded without any outcome (constant contribution). In case of false self-avoidance, several residues can then occupy one site, producing contacts no physical chain could form. Energy becomes artificially low, because it has a situation as a self-eating snake in the phone snake game. Therefore, my benchmark asserts that the energy found never falls below the published value.

2. What move set did you use for the stochastic search, and can it in principle reach every conformation? (Some common local move sets cannot — it is worth knowing whether yours is one of them.)

A random mixture of pivot moves (80%) and local moves (20%). Corner flips and tail moves are not ergodic, and a search built only on them would leave regions of conformation space unreachable. Pivot moves are ergodic for self-avoiding walks[2], and if pivot is available at every step, the mixture has that property. In principle, the search can therefore reach any conformation.

3. For the 20-mer, how confident are you that you found the true ground state rather than a local minimum, and what is your evidence? How does that confidence differ from the short sequence you enumerated?

The confidence is statistical and external, not exhaustive. Of 100 independent runs, 32 reached −9, and that value matches the published one[3]. 
The difference from a short sequence is a matter of scale. Self-avoiding walks on the square lattice grow roughly as 2.64^n: 284 for 6 residues, about 2,200 for 8, and about 3*10^8 for 20. At the short end, the space can be enumerated and the answer proven.  For example, at 20 residues it would be complex, because finding the ground state is NP-hard in general. My 100 runs of  near 150,000 steps each visit at most 1.5 * 10^8 conformations, well under 5% of the space even assuming no repeats. My argument is that repeated independent searches converge on the same value and that this value matches the literature.

4. Does your target sequence have a unique ground-state conformation or a degenerate one? Why does ground-state uniqueness matter if we want the model to behave like a real, foldable protein?

The 20-mer is non-degenerate in my results: all 32 runs reaching −9 gave one structure. S1-2 and S1-3 are strongly degenerate, with eleven distinct structures each at the minimum. So degeneracy is a property of the sequence. Uniqueness matters because a protein with many equally stable folds has no
defined native state. In nature, there is mostly no reliable function with multiple ground-states. Anfinsen's dogma holds that the native structure is determined by the sequence alone for the unique native state. Sequences with unique
ground-states are considered the model's
analogue of foldable proteins. Highly degenerate states can be assigned to sequences that behave
as random heteropolymers.

5. On a square lattice, the sites form two interleaved sublattices. What does that imply about which pairs of residues can ever be in contact, and how could you exploit it?

Coloring sites by the parity of `x + y` gives two sublattices, where every chain step moves to the other color. The color of residue i is therefore fixed by the parity of i. And adjacent sites always have opposite colors, two residues can be in contact only when `|i - j|` is odd. I used this rule in my cycles by iterating 'j' from 'i + 3' in steps of 2. It removes half of the candidate pairs. With a benchmark, it showed faster runtime, making this method cheaper.

6. If you ran Monte Carlo across a range of temperatures, what would you expect to see in the average energy, and what would that tell you about folding cooperativity?

Monte Carlo across a range of temperatures requires equilibrium sampling at fixed temperatures. And this is not about my annealing
trajectory. Running Metropolis Monte Carlo at a series of fixed temperatures,
discarding a burn-in at each, I would expect <E> to be near 0 at high temperature. While <E> obtains the minimum at low temperature, and some drop between them. If it is a cooperative transition that will affect the resulting 0/1 situation with two states without any intermediate in between. If it is a smooth gradient decline, it can show a continuous folding process with intermediates. To answer this, I would add a physical feature of the heat capacity that can describe this process. However, I didn't do this. It can be an extension to the current state.

7. What real feature of protein folding does this model genuinely capture, and what does it badly miss?

It misses almost everything else. The focus of HP-model only on hydrophobic interactions without other possible connections, such as hydrogen bonds, salt bridges, ionic interactions, polar interactions, S-S covalent bridges. Also, there are many amino residues with unique groups, making the possible conformation space absolutely large. But in HP-model twenty amino acid types are collapsed into two. Backbone geometry is replaced by lattice steps, side chains are absent, and solvent enters only implicitly through the H–H contact term. The final folding in reality is a balance of different interactions.

8. What approaches did you consider and reject along the way, and why? If you had another day on this, what would you try next?

If I had more time, I would add other methods that demonstrate excellent performance for such problems. For example, there is replica exchange Monte Carlo. It found global minima for lengths greater than 48 residues [4]. I would also modify my end criterion with a flexible number of steps based on the compactness of folding and the length of time without energy changes required to obtain the global minimum. I would also add heat capacity calculations and compactness measurements. One idea is to train a small machine learning model to predict future foldable ground-state energies based on compactness, H:P ratio, and other features. To do so, I would have to perform some brute-force calculations of exact ground-state conformers or find them in the literature. Then, I would train the model to predict the possible range of energies for future estimation using a stochastic method, if the range of energies is approximately true.

## Post-Deadline Updates

**Speed-Up Energy Calculations** Goes from O(n^2) to O(n). In new version each contact was found by scanning every pair of residues with hashing occupied lattice sites and checking four neighbors of each residue, that is O(n) as expected. Old version with O(n^2) was saved in 'energy.py' as 'energy_old' function and used to test the acceleration of module. 

Cost of a single energy() call for HPHP... seqeunce with different lengths (n = 10, 30, 50) has the following trend:

| n   | O(n^2), µs | O(n), µs | Acceleration |
|-----|-----------|----------|---------|
| 10  | 3.8       | 2.3      | 1.7x    |
| 30  | 32.5      | 5.5      | 5.9x    |
| 50  | 79.3      | 9.0      | 8.9x    |

Full ran with command `main.py --nstarts 100` (seed = 0) best of three: old = 121 seconds, new = 64 seconds. (1.88x acceleration). Ground-state energies are identical, as far as hit rate and mean energy on both approaches. The next bottleneck is move generations and self-avoidance criterion.  