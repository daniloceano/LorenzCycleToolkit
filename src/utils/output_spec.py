"""
This module describes what each per-level CSV file holds: the quantity, its
units, how the column value in the main results file is obtained from it,
and how the split terms combine. It is a declarative record of the output
convention. The computation does not read it; the vertical-level figures,
the documentation tests and the unit tests do.

Notation: [x] is the zonal mean, x* = x - [x] the zonal eddy, <x> the area
mean, [x]'' = [x] - <[x]>, sigma the static stability, g the gravitational
acceleration (``metpy.constants.g``) and p the pressure in Pa.

The column value of a term is ``factor * integral(profile dp)`` with the
trapezoidal rule over the pressure levels written in the file. The factor
is not uniform across terms: the energy and conversion profiles of most
terms already include every physical factor, but ``Kz``/``Ke`` are written
without ``1/(2g)`` and ``Ck`` without ``1/g``. ``M`` has its own column
relation (see its entry).

Verification of the stable entries (Wave 2 of the post-audit corrections):
the units were derived from the equations, captured with Pint at the moment
each profile is written, and checked against the written values for the
ERA5 fixed and moving references; the factors were checked against the
main results file of the reference runs.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

from metpy.constants import g

_VERIFIED = (
    "verified: units derived from the equations and captured with Pint when the "
    "profile is written; profile-to-column factor checked against the main "
    "results file"
)
_VERIFIED_SPLIT = (
    "verified: units derived from the equations and captured with Pint when the "
    "profile is written; combination checked against the written files"
)

VERTICAL_PROFILE_SPEC = {
    "Az": {
        "quantity": "Zonal available potential energy per unit pressure",
        "expression": "<[T]''^2> / (2 sigma)",
        "units": "J/m**2/Pa",
        "display_units": "J m⁻² Pa⁻¹",
        "label": "Zonal APE per unit pressure",
        "column_term": "Az",
        "column_units": "J/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Ae": {
        "quantity": "Eddy available potential energy per unit pressure",
        "expression": "<T*^2> / (2 sigma)",
        "units": "J/m**2/Pa",
        "display_units": "J m⁻² Pa⁻¹",
        "label": "Eddy APE per unit pressure",
        "column_term": "Ae",
        "column_units": "J/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Kz": {
        "quantity": (
            "Area mean of the squared zonal-mean wind components, equal to twice "
            "the zonal kinetic energy per unit mass; written without 1/(2g)"
        ),
        "expression": "<[u]^2 + [v]^2>",
        "units": "m**2/s**2",
        "display_units": "m² s⁻²",
        "label": r"$\langle [u]^2+[v]^2 \rangle$ = 2 × zonal KE per unit mass",
        "column_term": "Kz",
        "column_units": "J/m**2",
        "column_factor": "1/(2g)",
        "role": "scaled_integrand",
        "verification": _VERIFIED,
    },
    "Ke": {
        "quantity": (
            "Area mean of the squared eddy wind components, equal to twice "
            "the eddy kinetic energy per unit mass; written without 1/(2g)"
        ),
        "expression": "<u*^2 + v*^2>",
        "units": "m**2/s**2",
        "display_units": "m² s⁻²",
        "label": r"$\langle u^{*2}+v^{*2} \rangle$ = 2 × eddy KE per unit mass",
        "column_term": "Ke",
        "column_units": "J/m**2",
        "column_factor": "1/(2g)",
        "role": "scaled_integrand",
        "verification": _VERIFIED,
    },
    "Cz": {
        "quantity": "Conversion from Az to Kz per unit pressure (positive: Az to Kz)",
        "expression": "-Cz_1 * Cz_2 = -(R/(p g)) <[omega]''[T]''>",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Cz per unit pressure",
        "column_term": "Cz",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Ca": {
        "quantity": "Conversion from Az to Ae per unit pressure (positive: Az to Ae)",
        "expression": "-(Ca_1 + Ca_2)",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Ca per unit pressure",
        "column_term": "Ca",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Ce": {
        "quantity": "Conversion from Ae to Ke per unit pressure (positive: Ae to Ke)",
        "expression": "-Ce_1 * Ce_2 = -(R/(p g)) <omega* T*>",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Ce per unit pressure",
        "column_term": "Ce",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Ck": {
        "quantity": (
            "Conversion from Ke to Kz per unit pressure multiplied by g "
            "(positive: Ke to Kz); written without 1/g"
        ),
        "expression": "Ck_1 + Ck_2 + Ck_3 + Ck_4 + Ck_5",
        "units": "W/kg",
        "display_units": "W kg⁻¹",
        "label": "g × Ck per unit pressure",
        "column_term": "Ck",
        "column_units": "W/m**2",
        "column_factor": "1/g",
        "role": "scaled_integrand",
        "verification": _VERIFIED,
    },
    "C_overturning": {
        "quantity": "Domain-mean overturning pressure work per unit pressure (diagnostic, outside the budget)",
        "expression": "-<omega> alpha / g, alpha = R <T> / p",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "C_overturning per unit pressure",
        "column_term": "C_overturning",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Gz": {
        "quantity": "Generation of Az per unit pressure",
        "expression": "<[Q]''[T]''> / (cp sigma)",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Gz per unit pressure",
        "column_term": "Gz",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "Ge": {
        "quantity": "Generation of Ae per unit pressure",
        "expression": "<Q* T*> / (cp sigma)",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Ge per unit pressure",
        "column_term": "Ge",
        "column_units": "W/m**2",
        "column_factor": "1",
        "role": "integrand",
        "verification": _VERIFIED,
    },
    "M": {
        "quantity": "Mass-continuity residual: horizontal divergence plus d<omega>/dp",
        "expression": "<D_EW> + <D_NS> + d<omega>/dp",
        "units": "1/s",
        "display_units": "s⁻¹",
        "label": "Mass-continuity residual",
        "column_term": "M",
        "column_units": "kg/m**2/s",
        "column_factor": "special",
        "column_relation": (
            "M_column = (integral((<D_EW> + <D_NS>) dp) + omega_bottom - omega_top) / g; "
            "the vertical part uses the exact boundary difference, so "
            "integral(M dp) / g does not reproduce the column value"
        ),
        "role": "special",
        "verification": (
            "verified: units derived from the equations and captured with Pint when "
            "the profile is written; column relation checked with a synthetic field"
        ),
    },
    # Split terms: not plotted, archived for analysis.
    "Ca_1": {
        "quantity": "Horizontal part of the Ca integrand",
        "expression": "<v* T* d[T]''/dphi> / (a sigma)",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Ca_1",
        "part_of": "Ca",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
    "Ca_2": {
        "quantity": "Vertical part of the Ca integrand",
        "expression": "<omega* T* d[T]''/dp> / sigma",
        "units": "W/m**2/Pa",
        "display_units": "W m⁻² Pa⁻¹",
        "label": "Ca_2",
        "part_of": "Ca",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
    "Cz_1": {
        "quantity": "Pressure factor R/(p g) of Cz (depends on pressure only)",
        "expression": "R / (p g)",
        "units": "m/Pa/K",
        "display_units": "m Pa⁻¹ K⁻¹",
        "label": "Cz_1",
        "part_of": "Cz",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
    "Cz_2": {
        "quantity": "Covariance of zonal-mean omega and temperature deviations",
        "expression": "<[omega]''[T]''>",
        "units": "K*Pa/s",
        "display_units": "K Pa s⁻¹",
        "label": "Cz_2",
        "part_of": "Cz",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
    "Ce_1": {
        "quantity": "Pressure factor R/(p g) of Ce (depends on pressure only)",
        "expression": "R / (p g)",
        "units": "m/Pa/K",
        "display_units": "m Pa⁻¹ K⁻¹",
        "label": "Ce_1",
        "part_of": "Ce",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
    "Ce_2": {
        "quantity": "Covariance of eddy omega and eddy temperature",
        "expression": "<omega* T*>",
        "units": "K*Pa/s",
        "display_units": "K Pa s⁻¹",
        "label": "Ce_2",
        "part_of": "Ce",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    },
}
for _i, _expr in enumerate(
    [
        "<u* v* cos(phi) d([u]/cos(phi))/dphi> / a",
        "<v*^2 d[v]/dphi> / a",
        "<tan(phi) u*^2 [v]> / a",
        "<omega* u* d[u]/dp>",
        "<omega* v* d[v]/dp>",
    ],
    start=1,
):
    VERTICAL_PROFILE_SPEC[f"Ck_{_i}"] = {
        "quantity": f"Subterm {_i} of the Ck integrand multiplied by g",
        "expression": _expr,
        "units": "W/kg",
        "display_units": "W kg⁻¹",
        "label": f"Ck_{_i}",
        "part_of": "Ck",
        "role": "split_term",
        "verification": _VERIFIED_SPLIT,
    }

# How the split terms combine into the written profile of their term.
SPLIT_TERM_RELATIONS = {
    "Ca": {"sign": -1, "operation": "sum", "terms": ["Ca_1", "Ca_2"]},
    "Ck": {"sign": 1, "operation": "sum", "terms": ["Ck_1", "Ck_2", "Ck_3", "Ck_4", "Ck_5"]},
    "Cz": {"sign": -1, "operation": "product", "terms": ["Cz_1", "Cz_2"]},
    "Ce": {"sign": -1, "operation": "product", "terms": ["Ce_1", "Ce_2"]},
}

# Dz and De are written only by the direct dissipation pathway (without -r,
# with friction fields in the input). Their profile convention has NOT been
# verified; nothing here should be read as a confirmed unit or formula.
PROVISIONAL_PROFILE_SPEC = {
    "Dz": {
        "status": "NOT VERIFIED",
        "note": (
            "Written by the direct dissipation pathway, which needs friction "
            "fields and is not exercised by the reference cases. Units, level "
            "and column relation are under review."
        ),
    },
    "De": {
        "status": "NOT VERIFIED",
        "note": (
            "Written by the direct dissipation pathway, which needs friction "
            "fields and is not exercised by the reference cases. Units, level "
            "and column relation are under review."
        ),
    },
}

# Columns of the main results file that are not vertical integrals of a
# per-level file.
TENDENCY_UNITS = {
    "∂Az/∂t (finite diff.)": "W/m**2",
    "∂Ae/∂t (finite diff.)": "W/m**2",
    "∂Kz/∂t (finite diff.)": "W/m**2",
    "∂Ke/∂t (finite diff.)": "W/m**2",
}
TENDENCY_SCHEME = (
    "numpy.gradient over the elapsed time in seconds: second-order centred "
    "differences at interior times (non-uniform spacing allowed) and "
    "first-order one-sided differences at the first and last times"
)

PRESSURE_AXIS_LABEL = "Pressure (hPa)"

_FACTORS = {"1": 1.0, "1/g": 1.0 / float(g.magnitude), "1/(2g)": 1.0 / (2.0 * float(g.magnitude))}


def column_factor(term):
    """
    Return the number that turns the pressure integral of a profile into its
    column value, or None when the relation is not a constant factor.

    Args:
        term (str): Name of the per-level file (e.g. ``"Kz"``).

    Returns:
        float or None: ``factor`` in ``column = factor * integral(profile dp)``
        (pressure in Pa), or None for ``M`` and for split terms.
    """
    return _FACTORS.get(VERTICAL_PROFILE_SPEC[term].get("column_factor"))


def profile_label(term):
    """
    Return the axis or colour-bar label of a per-level profile.

    Args:
        term (str): Name of the per-level file.

    Returns:
        str: The quantity and its units, e.g.
        ``"Zonal APE per unit pressure (J m⁻² Pa⁻¹)"``.

    Raises:
        KeyError: If the term has no stable entry, so that a figure never
        falls back to a generic label.
    """
    spec = VERTICAL_PROFILE_SPEC[term]
    return f"{spec['label']} ({spec['display_units']})"


def pressure_tick_label(pressure_pa):
    """
    Format a pressure in Pa as a tick label in hPa.

    Only the label is converted; the plotted coordinate stays in Pa.

    Args:
        pressure_pa (float): Pressure in Pa.

    Returns:
        str: The pressure in hPa, without trailing zeros (``70000.0`` gives ``"700"``).
    """
    return f"{float(pressure_pa) / 100.0:g}"
