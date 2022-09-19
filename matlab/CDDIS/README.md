## Description

The CDDIS directory is used to store ephemeris files files fetched from CDDIS. 

## Usage

* Run the `fetch.py` script available in `/pyFunctions`

Note:

A user can change: `YYYY`, `DOY`, `HH`, and `IGSStationInitials` in `fetch.py` to fetch a different ephemeris file.

## Pre-requisites

* `Python 3.8.9` or greater is installed
* Register with the CDDIS earthdata login
* .netrc file is installed in the home directory with your registration details above as given below:
> machine urs.earthdata.nasa.gov login <username> password <pword>