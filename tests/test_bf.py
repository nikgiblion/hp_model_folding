"""
Brute-force enumeration of SAWs.
"""

import pytest
from conformation import final_validation_sequence
from bf import iter_saws, iter_symmetry_saws
from symmetry import canonical_form

#SAW space of k residues (1, 2, 3, 4, 5, 6) where n = k + 1 (source: https://oeis.org/A001411)
SAW_COUNTS = [4, 12, 36, 100, 284, 780]

@pytest.mark.parametrize("k", range(1, 7))
def test_saw_count_matches(k):
    assert sum(1 for _ in iter_saws(k + 1)) == SAW_COUNTS[k - 1]

@pytest.mark.parametrize("n", range(1, 9))
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

@pytest.mark.parametrize("k", range(1,7))
def test_symmetry_count_matches(k):
    expected = 1 + (SAW_COUNTS[k - 1] - 4) // 8
    assert sum(1 for _ in iter_symmetry_saws(k + 1)) == expected

@pytest.mark.parametrize("n", range(1, 7))
def test_canonical_walks_are_one_per_symmetry_class(n):
    canonical = [canonical_form(w, allow_reverse=False) for w in iter_symmetry_saws(n)]
    everything = {canonical_form(w, allow_reverse=False) for w in iter_saws(n)}
    assert len(set(canonical)) == len(canonical)   # no class is produced twice
    assert set(canonical) == everything            # no class is missed

@pytest.mark.parametrize("n", range(2, 6))
def test_canonical_walks_are_valid_and_start_right(n):
    for w in iter_symmetry_saws(n):
        assert w[:2] == ((0, 0), (1, 0))
        assert final_validation_sequence(w)


def test_canonical_short_chains():
    assert list(iter_symmetry_saws(1)) == [((0, 0),)]
    assert list(iter_symmetry_saws(2)) == [((0, 0), (1, 0))]
    assert list(iter_symmetry_saws(3)) == [((0, 0), (1, 0), (2, 0)), ((0, 0), (1, 0), (1, 1))]
