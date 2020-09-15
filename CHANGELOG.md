# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.3.0] - 2020-09-15

### Changed
 - KerrGeoOrbitFunction now outputs Frequencies as an Association.
 - Trajectories are returned as pure functions.
 - Trajectories are Listable.

### Fixed
 - Assorted bug fixes needed to pass unit tests
 - Fixed a problem where loading the package with Needs would generate an error message.


## [0.2.0] - 2020-05-18

### Added
 - User-friendly output format for KerrGeoOrbitFunction and KerrParallelTransportFrameFunction.


## [0.1.0] - 2020-05-18
 - Initial versioned release. Previous changes:
   - 2019-03-25: Add proper time frequencies (thanks Maarten van de Meent). Break code into separate files to make it more manageable. Also make KerrGeoConstantsOfMotion and KerrGeoFrequencies return associations (thanks Barry Wardell). Add unit tests (thanks Sam Upton and Ollie Long).
   - 2018-11-22: Merged the development branch into the master branch. Thanks to Zach, Tommy and Chuck for contributing code. This update is a major change with a lot of improvements. A key difference is the defintion of the inclination angle which is now x_inc = Cos[\theta_inc]. Bumped version number up to 0.5.
   - 2018-08-07: Added generic orbit calculation (thanks for M. van de Meent for contributing code). Currently the code is separated in KerrGeoOrbit2[..] but will soon replace the earlier code
   - 2017-09-07: Initial version publicly released.
   - 2017-06-10: Initial version created.


