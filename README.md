# KerrGeodesics

Version 0.5
Copyright 2017 Niels Warburton

The KerrGeodesics package provides a set of functions for calculating
properties of bound timelike geodesic orbits about a Kerr black hole.


### Getting the package

The latest development version will always be available from the project git
repository:

git clone https://github.com/BlackHolePerturbationToolkit/KerrGeodesics.git


### Requirements

Mathematica: KerrGeodesics requires a recent version of Mathematica. It is 
typically tested with only the latest available version.


### Installation

Clone the repository and place it somewhere on Mathematica's $Path.
Typical locations are inside ${HOME}/.Mathematica/Applications/ for Linux or
inside ${HOME}/Library/Mathematica/Applications/ for Mac OSX.


### Usage

The package may be loaded into Mathematica using the command:

<< KerrGeodesics`


### Examples

Examples are included in the documentation. See the
KerrGeodesics page in Documentation Center.


### Changelog

- 25 March 2019: Add propertime frequencies (thanks Maarten van de Meent). Break code into separate files to make it more manageable. Also make KerrGeoConstantsOfMotion and KerrGeoFrequencies return associations (thanks Barry Wardell). Add unit tests (thanks Sam Upton and Ollie Long).
- 22 November 2018: Merged the development branch into the master branch. Thanks to Zach, Tommy and Chuck for contributing code. This update is a major change with a lot of improvements. A key difference is the defintion of the inclination angle which is now x_inc = Cos[\theta_inc]. Bumped version number up to 0.5.
- 7 August 2018: Added generic orbit calculation (thanks for M. van de Meent for contributing code). Currently the code is separated in KerrGeoOrbit2[..] but will soon replace the earlier code
- 7 September 2017: Initial version publicly released
- 10 June 2017: Initial version created.


### Known problems

Known bugs are recorded in the project bug tracker:

https://github.com/BlackHolePerturbationToolkit/KerrGeodesics/issues


### License

This code is distributed under the University of Illinois/NCSA
Open Source License. Details can be found in the LICENSE file.


### Authors

Niels Warburton  
Maarten van de Meent  
Zach Nasipak  
Thomas Osburn  
Charles Evans  
