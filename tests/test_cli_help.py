"""
Help text of the command line for the data sources (stage 4d).

Checks the two findings of stage 4c: ``--cdsapi`` is no longer marked as
experimental, and ``--time-resolution`` applies to every data source and states
the restriction of each route (``--arco``: multiples of 6 h only). No file is
read or written and nothing is downloaded.
"""

import lorenzcycletoolkit


def option_help(name):
    """Help text of one option, with whitespace normalised."""
    parser = lorenzcycletoolkit.create_arg_parser()
    action = next(a for a in parser._actions if name in a.option_strings)
    return " ".join(action.help.replace("%%", "%").split())


def full_help():
    return " ".join(lorenzcycletoolkit.create_arg_parser().format_help().split())


def test_cdsapi_is_not_marked_experimental():
    text = option_help("--cdsapi")
    assert "experimental" not in text.lower()
    assert "CDS API" in text


def test_time_resolution_applies_to_every_route_with_its_restrictions():
    text = option_help("--time-resolution")
    for route in ("--cdsapi", "--ncar", "--wekeo", "--arco", "--destine"):
        assert route in text
    assert "default: 3" in text
    assert "--arco only multiples of 6" in text
    assert "any whole number of hours" in text
    assert "CDS API data download" not in text  # the old text, limited to --cdsapi


def test_only_ncss_remains_experimental_in_the_data_sources():
    text = full_help().lower()
    assert text.count("experimental") == 1 and "ncss (experimental)" in text
