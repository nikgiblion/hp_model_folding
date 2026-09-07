"""D4 symmetries, canonical forms and degeneracy counting."""

from symmetry import SYMMETRIES, canonical_form

SQUARE = [(0, 0), (1, 0), (1, 1), (0, 1)]
L_SHAPE = [(0, 0), (1, 0), (2, 0), (2, 1)]


def test_the_group_has_eight_distinct_elements():
    # D4: four rotations, each with and without reflection.
    assert len({t(3, 7) for t in SYMMETRIES}) == 8


def test_canonical_form_is_translation_invariant():
    shifted = [(x + 5, y + 3) for x, y in SQUARE]
    assert canonical_form(SQUARE) == canonical_form(shifted)


def test_canonical_form_is_rotation_invariant():
    rotated = [(-y, x) for x, y in SQUARE]
    assert canonical_form(SQUARE) == canonical_form(rotated)


def test_canonical_form_is_reflection_invariant():
    mirrored = [(-x, y) for x, y in SQUARE]
    assert canonical_form(SQUARE) == canonical_form(mirrored)


def test_every_symmetry_maps_a_conformation_to_the_same_class():
    # The point of the canonical form: all eight images collapse into one.
    images = {canonical_form([t(x, y) for x, y in SQUARE]) for t in SYMMETRIES}
    assert len(images) == 1


def test_distinct_shapes_have_distinct_canonical_forms():
    assert canonical_form(SQUARE) != canonical_form(L_SHAPE)


def test_canonical_form_is_idempotent():
    once = canonical_form(SQUARE)
    assert canonical_form(once) == once