"""
Brute-force enumeration of SAWs.
"""

import pytest
import random
from conformation import final_validation_sequence
from bf import iter_saws, iter_symmetry_saws, exact_ground_states, exact_distinct_ground_states
from symmetry import canonical_form
from energy import energy, energy_old

#SAW space of k residues (1, 2, 3, 4, 5, 6) where n = k + 1 (source: https://oeis.org/A001411).
SAW_COUNTS = [4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100]

@pytest.mark.parametrize("k", range(1, 11))
def test_saw_count_matches(k):
    assert sum(1 for _ in iter_saws(k + 1)) == SAW_COUNTS[k - 1]

@pytest.mark.parametrize("n", range(1, 11))
def test_every_walk_is_valid_and_unique(n):
    walks = list(iter_saws(n))
    assert all(len(w) == n and w[0] == (0, 0) for w in walks)
    assert all(final_validation_sequence(w) for w in walks)
    assert len(set(walks)) == len(walks)

def test_single_residue():
    assert list(iter_saws(1)) == [((0, 0),)]

def test_empty_chain_is_rejected():
    with pytest.raises(ValueError):
        list(iter_saws(0))

@pytest.mark.parametrize("k", range(1,11))
def test_symmetry_count_matches(k):
    expected = 1 + (SAW_COUNTS[k - 1] - 4) // 8
    assert sum(1 for _ in iter_symmetry_saws(k + 1)) == expected

@pytest.mark.parametrize("n", range(1, 11))
def test_canonical_walks_are_one_per_symmetry_class(n):
    canonical = [canonical_form(w, allow_reverse=False) for w in iter_symmetry_saws(n)]
    everything = {canonical_form(w, allow_reverse=False) for w in iter_saws(n)}
    assert len(set(canonical)) == len(canonical)   # no class is produced twice
    assert set(canonical) == everything            # no class is missed

@pytest.mark.parametrize("n", range(2, 11))
def test_canonical_walks_are_valid_and_start_right(n):
    for w in iter_symmetry_saws(n):
        assert w[:2] == ((0, 0), (1, 0))
        assert final_validation_sequence(w)


def test_canonical_short_chains():
    assert list(iter_symmetry_saws(1)) == [((0, 0),)]
    assert list(iter_symmetry_saws(2)) == [((0, 0), (1, 0))]
    assert list(iter_symmetry_saws(3)) == [((0, 0), (1, 0), (2, 0)), ((0, 0), (1, 0), (1, 1))]

def random_sequences(n, count, seed):
    rng = random.Random(seed)
    return ["".join(rng.choice("HP") for _ in range(n)) for _ in range(count)]


@pytest.mark.parametrize("n", range(4, 11))
def test_fast_energy_matches_naive_on_all_walks(n):
    # O(n) hash-map energy against the O(n^2) double loop.
    for sequence in random_sequences(n, 3, seed=n):
        for w in iter_symmetry_saws(n):
            assert energy(sequence, w) == energy_old(sequence, w), (sequence, w)


@pytest.mark.parametrize("n", range(4, 11))
def test_incremental_search_matches_full_recomputation(n):
    # Incremental energy inside the DFS against energy() recomputed on every leaf.
    for sequence in random_sequences(n, 3, seed=100 + n):
        walks = list(iter_symmetry_saws(n))
        energies = [energy(sequence, w) for w in walks]
        e_min = min(energies)
        expected = {w for w, e in zip(walks, energies) if e == e_min}
        e_exact, optimal = exact_ground_states(sequence)
        assert e_exact == e_min, sequence
        assert set(optimal) == expected, sequence
        assert len(optimal) == len(expected)   # no walk reported twice.


def test_hpph_square():
    assert exact_ground_states("HPPH") == (-1.0, [((0, 0), (1, 0), (1, 1), (0, 1))])


@pytest.mark.parametrize("n", [1, 2, 3, 6, 9, 11])
def test_all_polar_chain_has_zero_energy(n):
    e, optimal = exact_ground_states("P" * n)
    assert e == 0.0
    assert len(optimal) == sum(1 for _ in iter_symmetry_saws(n))   # every walk is optimal.


@pytest.mark.parametrize("n", range(4, 11))
def test_pruning_does_not_change_the_answer(n):
    # Branch and bound must keep every tie
    for sequence in random_sequences(n, 3, seed=200 + n):
        assert exact_ground_states(sequence, prune=True) == exact_ground_states(sequence, prune=False)


@pytest.mark.slow
def test_benchmark_20mer_ground_state():
    e, optimal = exact_ground_states("HPHPPHHPHPPHPHHPPHPH")
    assert e == -9.0
    assert len(optimal) == 2

def palindromic_sequences(n, count, seed):
    rng = random.Random(seed)
    halves = ["".join(rng.choice("HP") for _ in range((n + 1) // 2)) for _ in range(count)]
    return [h + h[::-1][n % 2:] for h in halves]


@pytest.mark.parametrize("n", range(4, 11))
def test_reversal_merges_nothing_for_non_palindromes(n):
    for sequence in random_sequences(n, 5, seed=300 + n):
        if sequence == sequence[::-1]:
            continue
        _, d4 = exact_ground_states(sequence)
        _, distinct = exact_distinct_ground_states(sequence)
        assert len(distinct) == len(d4), sequence


@pytest.mark.parametrize("n", range(4, 11))
def test_reversal_merges_only_reversal_images(n):
    for sequence in palindromic_sequences(n, 3, seed=400 + n):
        _, d4 = exact_ground_states(sequence)
        _, distinct = exact_distinct_ground_states(sequence)
        assert len(distinct) <= len(d4)
        # every dropped conformation is a lattice image of some kept one, reversed.
        kept = {canonical_form(c, allow_reverse=True) for c in distinct}
        assert all(canonical_form(c, allow_reverse=True) in kept for c in d4)


def test_palindrome_ground_states_collapse_under_reversal():
    sequence = "HHPHPHPHH"   # palindromic
    assert len(exact_ground_states(sequence)[1]) == 2
    assert len(exact_distinct_ground_states(sequence)[1]) == 1


@pytest.mark.slow
def test_benchmark_20mer_is_non_degenerate_under_reversal():
    sequence = "HPHPPHHPHPPHPHHPPHPH"
    assert sequence == sequence[::-1]
    assert len(exact_ground_states(sequence)[1]) == 2
    assert len(exact_distinct_ground_states(sequence)[1]) == 1