# HP model[1,2] of protein folding on a 2D lattice with binding extension

Search for the ground state by simulated annealing, with ground-state degeneracy
analysis and ligand binding site detection.

## Model

A protein sequence is reduced to a chain of two residue types : **H** (hydrophobic) and
**P** (polar). To fold the conformer, they are placed as a self-avoiding walk on a 2D square lattice. Two H
residues on adjacent lattice sites that are not neighbors along the chain form
a contact worth −1; H–P and P–P contacts contribute nothing. The energy of a
conformation is therefore the negated count of non-bonded H–H contacts, and the
native fold is the conformation minimising it.

## Installation

Python 3.9+, standard library only. Nothing to install.

Running the tests requires pytest: `pip install -r requirements-dev.txt`, then `pytest`.

## Usage

```bash
python main.py                                        # 20-mer benchmark (from task), 20 independent runs
python main.py --sequence "H2(P2H)7H" --nstarts 20    # compressed notation (quotes required)
python main.py --nstarts 20 --binding                 # + cavities and ligand binding check analysis
python benchmark.py --max-length 25 --nstarts 100     # compare against published values
```

Important: sequences containing parentheses must be quoted, otherwise the shell tries to
expand them.

Full parameter list is avaliable with command: `python main.py --help`, `python benchmark.py --help`.

## Project structure

| File/Dir | Purpose |
|------|---------|
| `tests/` | pytest suite |
| `tools/` | timing scripts |
| `conformation.py` | lattice representation, validity checks, move set, sequence parsing |
| `energy.py` | contact table and energy evaluation |
| `annealing.py` | Metropolis criterion, annealing schedule, multi-start driver |
| `symmetry.py` | D4 symmetries, canonical forms, degeneracy counting |
| `binding.py` | cavity detection by flood fill, ligand binding energies |
| `visualization.py` | ASCII rendering of conformations |
| `reporting.py` | text reports and energy histograms |
| `benchmark.py` | benchmark sequences and comparison table |
| `main.py` | command-line entry point |

## Results

Benchmark against published ground-state energies (100 independent runs per
sequence, single seed):

| ID | Length | Published | Found | hit % | Mean E |
|----|--------|-----------|-------|-------|--------|
| HPPHPH | 6 | −2 | −2 | 100% | −2.00 |
| Unger–Moult | 20 | −9 | −9 | 36% | −8.34 |
| S1-2 | 24 | −9 | −9 | 35% | −8.33 |
| S1-3 | 25 | −8 | −8 | 23% | −7.07 |
| S1-4 | 36 | −14 | −14 | 2% | −11.55 |
| S1-5 | 48 | −23 | −22 | 0% | −18.03 |

Up to n = 25 residues the published minimum is reached commonly. Beyond 36 residues
the method reaches it only occasionally, and past 48 it no longer does. Because published
values for longer chains come from replica exchange methods Monte Carlo with far larger
compute budgets. See `WRITEUP.md` for the reasoning behind this limit.

## Departures from the canonical HP model

The ligand binding extension uses its own interaction table in which a P ligand
in a polar pocket is favourable (−1). In the canonical model P–P contacts are
neutral. Keeping that here would make a polar ligand indifferent to every pocket
and the analysis vacuous. Protein–protein energies are unchanged (−1 / 0 / 0).

## References

1. Lau, K. F. & Dill, K. A. (1989). A lattice statistical mechanics model of the
   conformational and sequence spaces of proteins. *Macromolecules* **22**,
   3986–3997. [doi:10.1021/ma00200a030](https://doi.org/10.1021/ma00200a030)
2. Madras, N. & Sokal, A. D. (1988). The pivot algorithm: a highly efficient
   Monte Carlo method for self-avoiding walks. *J. Stat. Phys.* **50**, 109–186. [doi:10.1007/BF01022990](https://doi.org/10.1007/BF01022990)
3. Unger, R. & Moult, J. (1993). Genetic algorithms for protein folding
   simulations. *J. Mol. Biol.* **231**, 75–81.
   [doi:10.1006/jmbi.1993.1258](https://doi.org/10.1006/jmbi.1993.1258)
4. Thachuk, C., Shmygelska, A. & Hoos, H. H. (2007). A replica exchange Monte
   Carlo algorithm for protein folding in the HP model. *BMC Bioinformatics*
   **8**, 342. [doi:10.1186/1471-2105-8-342](https://doi.org/10.1186/1471-2105-8-342)
5. Zhang, J., Kou, S. C. & Liu, J. S. (2007). Biopolymer structure simulation
   and optimization via fragment regrowth Monte Carlo. *J. Chem. Phys.* **126**,
   225101. [doi:10.1063/1.2736681](https://doi.org/10.1063/1.2736681)
