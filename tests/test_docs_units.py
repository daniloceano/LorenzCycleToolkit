#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The table of per-level quantities and units in ``docs/source/results.rst``
must agree with ``src/utils/output_spec.py``, and the column units of the
spec must agree with the units the term classes convert their column
values to (``OUTPUT_UNITS``).
"""

import csv
from pathlib import Path

import pint

from src.analysis.conversion_terms import ConversionTerms
from src.analysis.energy_contents import EnergyContents
from src.analysis.generation_and_dissipation_terms import GenerationDissipationTerms
from src.analysis.mass_continuity import MassContinuity
from src.utils.output_spec import (PROVISIONAL_PROFILE_SPEC, TENDENCY_UNITS,
                                   VERTICAL_PROFILE_SPEC)

RESULTS_RST = Path(__file__).resolve().parents[1] / "docs" / "source" / "results.rst"
COLUMNS = ["File", "Quantity written", "Units", "Column term", "Column units", "Factor"]


def _units_table():
    lines = RESULTS_RST.read_text(encoding="utf-8").splitlines()
    start = lines.index(".. csv-table::")
    header = lines[start + 1].strip()
    assert header.startswith(":header:")
    assert next(csv.reader([header[len(":header:"):].strip()], skipinitialspace=True)) == COLUMNS
    body = []
    for line in lines[start + 2:]:
        if line.strip().startswith(":"):
            continue
        if not line.strip():
            if body:
                break
            continue
        body.append(line.strip())
    return {row[0]: dict(zip(COLUMNS, row)) for row in csv.reader(body, skipinitialspace=True)}


def test_documented_table_matches_the_spec():
    ureg = pint.UnitRegistry()
    table = _units_table()
    assert sorted(table) == sorted(VERTICAL_PROFILE_SPEC)
    for term, row in table.items():
        spec = VERTICAL_PROFILE_SPEC[term]
        assert row["Units"] == spec["units"], term
        ureg.Quantity(1, row["Units"])
        assert row["Column term"] == spec.get("column_term", ""), term
        assert row["Column units"] == spec.get("column_units", ""), term
        assert row["Factor"] == spec.get("column_factor", ""), term


def test_unverified_profiles_are_not_in_the_table():
    assert not set(_units_table()) & set(PROVISIONAL_PROFILE_SPEC)


def test_column_units_of_the_spec_match_output_units_of_the_term_classes():
    ureg = pint.UnitRegistry()
    owner = {
        **{t: EnergyContents for t in ("Az", "Ae", "Kz", "Ke")},
        **{t: ConversionTerms for t in ("Cz", "Ca", "Ce", "Ck", "C_overturning")},
        **{t: GenerationDissipationTerms for t in ("Gz", "Ge")},
        "M": MassContinuity,
    }
    for term, cls in owner.items():
        doc = ureg.Quantity(1, VERTICAL_PROFILE_SPEC[term]["column_units"])
        out = ureg.Quantity(1, cls.OUTPUT_UNITS)
        assert doc.to(out.units).magnitude == 1.0, term


def test_tendency_units_and_scheme_are_documented():
    text = RESULTS_RST.read_text(encoding="utf-8")
    for name in TENDENCY_UNITS:
        assert name.replace(" (finite diff.)", "") in text
    assert "Tendencies of the energy reservoirs (W/m², i.e. J m⁻² s⁻¹)" in text
    assert "first-order one-sided differences at the first and last times" in text
