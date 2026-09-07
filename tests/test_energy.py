"""Energy calculation on 2D lattice."""

import pytest
from conformation import linear_conformation
from energy import ENERGY_CONTACT_PAIRS, all_contacts, energy, hh_contacts

# A 2x2 square
SQUARE = [(0, 0), (1, 0), (1, 1), (0, 1)]


def test_extended_chain_has_no_contacts():
    assert energy("HHHH", linear_conformation(4)) == 0.0


def test_square_has_one_hh_contact():
    assert energy("HHHH", SQUARE) == -1.0


def test_only_the_contacting_pair_matters():
    # Residues 0 and 3 are H, so the contact still counts.
    assert energy("HPPH", SQUARE) == -1.0


def test_polar_pair_in_contact_costs_nothing():
    # Here the H residues are chain neighbours, which never counts.
    assert energy("PHHP", SQUARE) == 0.0


def test_contact_list_on_the_square():
    assert all_contacts(SQUARE) == [(0, 3)]


def test_hh_contacts_filters_by_residue_type():
    assert hh_contacts("HHHH", SQUARE) == [(0, 3)]
    assert hh_contacts("HPPH", SQUARE) == [(0, 3)]
    assert hh_contacts("PHHP", SQUARE) == []


def test_length_mismatch_is_rejected():
    with pytest.raises(ValueError):
        energy("HHH", SQUARE)

