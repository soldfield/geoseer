
# geoseer
![Licence](https://img.shields.io/github/license/soldfield/geoseer.svg)
[![PyPI][https://img.shields.io/pypi/v/geoseer.svg)](https://pypi.python.org/pypi/geoseer)

The geoseer package provides a set of useful tools for common tasks undertaken by subsurface and surface professionals concerned with spatial data. The package talkes it's name from these individuals who may consider some element of what they work on to be a 'geo*' subject like: geology, geography, geophysics etc...

## nav_data: Navigation data calculations
This tool accepts two inputs, a list of markers with the geospatial coordinates (X & Y locations) and a list of lines
that are defined by those marker locations. This tool works for 2D seismic data or GPR surveys, where a larger survey
may be defined by a subset of well-defined points.

## plain text file scanning and loading
WORK IN PROGRESS
This tool extracts data from plain text files formatted with header fields above data tables.

## poly2fpq: Conversion to FracPaQ format
Tool for the conversion of files containing polylines to FracPaQ formats (Healy et al., 2017). Designed for use of
linestrings representing fracture traces.
