"""Sequence parsing, lattice representation and the move set."""

import random
import pytest
from conformation import (
    connectivity_check,
    corner_flip,
    expand_sequence,
    final_validation_sequence,
    linear_conformation,
    pivot_move,
    self_avoiding_walking_check,
    sequence_check,
    tail_move,
)


def test_sequence_is_trimmed_and_upcased():
    assert sequence_check(" hphp ") == "HPHP"


def test_forbidden_residue_is_rejected():
    with pytest.raises(ValueError):
        sequence_check("HPXH")


def test_empty_sequence_is_rejected():
    with pytest.raises(ValueError):
        sequence_check("   ")


def test_repeat_counts_expand():
    assert expand_sequence("H2P3") == "HHPPP"


def test_bracketed_group_expands():
    assert expand_sequence("H2(P2H)7H") == "HH" + "PPH" * 7 + "H"
    assert len(expand_sequence("H2(P2H)7H")) == 24


def test_unpaired_bracket_is_rejected():
    with pytest.raises(ValueError):
        expand_sequence("H2(P2H")


def test_linear_conformation_runs_along_the_x_axis():
    assert linear_conformation(3) == [(0, 0), (1, 0), (2, 0)]


def test_closed_square_is_connected_but_not_self_avoiding():
    # The chain returns to its own starting site.
    closed = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    assert connectivity_check(closed)
    assert not self_avoiding_walking_check(closed)
    assert not final_validation_sequence(closed)


def test_broken_chain_is_not_connected():
    assert not connectivity_check([(0, 0), (5, 5)])


def test_pivot_move_preserves_length_and_connectivity():
    rng = random.Random(0)
    start = linear_conformation(10)
    for _ in range(200):
        moved = pivot_move(start, rng)
        assert len(moved) == len(start)
        assert connectivity_check(moved)


def test_tail_move_changes_one_endpoint_only():
    rng = random.Random(1)
    start = linear_conformation(6)
    moved = tail_move(start, rng)
    differing = [i for i in range(len(start)) if start[i] != moved[i]]
    assert len(differing) == 1
    assert differing[0] in (0, len(start) - 1)
    assert connectivity_check(moved)


def test_corner_flip_returns_none_on_a_straight_segment():
    # Every interior site of a linear chain is straight, so no flip exists.
    rng = random.Random(2)
    assert corner_flip(linear_conformation(5), rng) is None