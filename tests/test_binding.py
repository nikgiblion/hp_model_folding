"""Cavity detection and ligand binding energies."""

from binding import check_cavities, inside_cavities
from conformation import linear_conformation

# A closed ring of eight amino residues around a single empty site.
RING = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def test_ring_encloses_one_cavity():
    assert inside_cavities(RING) == {(1, 1)}


def test_extended_chain_encloses_nothing():
    assert inside_cavities(linear_conformation(10)) == set()


def test_compact_chain_without_a_hole_encloses_nothing():
    # A 2x2 block is compact but has no interior site.
    assert inside_cavities([(0, 0), (1, 0), (1, 1), (0, 1)]) == set()


def test_hydrophobic_ligand_in_a_hydrophobic_pocket():
    # Four H residues face the cavity, each worth -1.
    result = check_cavities("HHHHHHHH", RING, ligand="H")
    assert result[0]["energy"] == -4.0


def test_hydrophobic_ligand_in_a_polar_pocket_is_neutral():
    assert check_cavities("PPPPPPPP", RING, ligand="H")[0]["energy"] == 0.0
