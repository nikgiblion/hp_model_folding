"""
Brute-force enumeration of SAWs.
"""

import pytest
from conformation import final_validation_sequence
from bf import iter_saws

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


