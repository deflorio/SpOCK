# SpOCK
The Spacecraft Orbital Computations Kit (SpOCK) is an open source spacecraft simulator for space mission analysis, space operations and GNC systems development and testing.
SpOCK allows high fidelity simulation of spacecraft's orbital and attitude dynamics, hardware, and the computations of mission events (ground station contacts, payload data-takes and eclipses). The ambition of the author is to develop a development and operational tool which includes only few but accurate simulations features and which can be easily further developed and customized. For this reason the software is developed with a (hopefully) clear structure with elements which can be easily read and reproduced and without a graphical user interface. For detailed information about installation, configuration and use please read the user guide.

---

## Features

*Spacecraft hardware*

Model of 6 spacecraft faces, 3 solar panels, 1 Sun camera, 1 Earth camera, 3 magnetometers, 6 coarse Sun sensors, 3-axis rate sensor, 3 reaction wheels, 3 magnetorquers, 2 types of orbit control propulsion systems.

*Dynamics*

High fidelity orbit and attitude numerical propagation, mission events computation, atmospheric drag and solar radiation pressure with panels model (dependent on attitude) or reference area model, attitude and orbital maneuvers with different type of actuators. The possibility is given to build a light version of the orbit propagator by using hardcoded perturbations models (gravitational field, drag and third body). This can be useful for embedded applications.

*Orbit environment models*

Orbit environment models: gravity field (EIGEN-6S, GGM02C, GGM03S, GGM03C), atmospheric density (JB2008, NRLMSIS-2.1, Harris-Priester), Earth's magnetic field (IGRF13, WMM2020), solar radiation pressure, third body (JPL planet ephemerides, low precision analytical formulas), gravity gradient.

*Time and Reference Systems*

Reach library for time, coordinate and parameterization systems transformations.

*Customizability*

Easy implementation of new spacecraft hardware, ready-to-use insertion of attitude and orbit control modules, set up for hardware-in-the-loop simulations.

---

## System Requirements

Linux Ubuntu 19 or higher

GNU make (sudo apt-get install build-essential)

gcc >= 10 supporting -std=c++20 (sudo apt-get install build-essential)

Libraries usually not included by default in Linux Ubuntu distribution: libboost-all-dev libxerces-c-dev xsdcxx gfortran libglut3.12 libglut-dev freeglut3-dev libxi-dev libxmu-dev mesa-utils libsdl2-2.0-0 libsdl2-dev libsoil* doxygen graphviz. These libraries (~1000 MB of addtional disk space) will be installed automatically during the installation process.

---

## Installation

*In root directory:*

make install_libs (only first compilation or after OS upgrade)

make install_atmo (only first compilation or after OS upgrade)

make install_mag (only first compilation or after OS upgrade)

make orbit

make sgp4

make attitude

make events

make tools

All the executable will be generated in the 'bin' folder. Command 'make clean' deletes all executables in 'bin' folder and folder 'obj'.

Download most recent de###.bsp file from SPICE library kernels ftp (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/) and put it into folder data/cspice (file de###.bsp is too big for the repository)

*Remarks:*

Installation and run were tested successfully on different real and virtual machines (VMware on Windows 10) running Linux Ubuntu 19.04 or higher. It is likely but not guaranteed that the installation on Ubuntu < 19 or other Linux distributions could work

Coomands make install_atmo and install_mag are used to compile external Fortran libraries and could give some Fortran compilation warnings.

---

## Configuration and Use

Refer to SPOCK_UserGuide_v1.7.pdf for the correct setup of simulation parameters and inputs/outputs paths in input files simparam_sample.xml and eventsparam_sample.xml, and the correct arguments of the executables.

---

## Run

*In linux terminal:*

Orbit numerical propagation: &nbsp;&nbsp; ./bin/OrbitPropagator &nbsp;&nbsp;&nbsp;&nbsp; input/SimulationParameters/simparam_sample.xml

SGP4 propagation: &nbsp;&nbsp; ./bin/SGP4Propagator &nbsp;&nbsp;&nbsp;&nbsp; input/SimulationParameters/simparam_sample.xml

Attitude numerical propagation: &nbsp;&nbsp; ./bin/AttitudePropagator &nbsp;&nbsp;&nbsp;&nbsp; input/SimulationParameters/simparam_sample.xml

Events computation: &nbsp;&nbsp; ./bin/EventsComputation &nbsp;&nbsp;&nbsp;&nbsp; input/SimulationParameters/eventsparam_sample.xml

Orbit evaluation: &nbsp;&nbsp; ./bin/EphEval &nbsp;&nbsp;&nbsp;&nbsp; /your_absolute_path/YourEphemeris.csv

Orbit interpolation: &nbsp;&nbsp; ./bin/EphIntpl &nbsp;&nbsp;&nbsp;&nbsp; /your_absolute_path/YourEphemeris.csv &nbsp;&nbsp;&nbsp;&nbsp; IntplType &nbsp;&nbsp;&nbsp;&nbsp; IntplPoints &nbsp;&nbsp;&nbsp;&nbsp; IntplStep

Conversion to sp3 format: &nbsp;&nbsp; ./bin/Eph2sp3 &nbsp;&nbsp;&nbsp;&nbsp; /your_absolute_path/YourEphemeris.csv &nbsp;&nbsp;&nbsp;&nbsp; PosVel &nbsp;&nbsp;&nbsp;&nbsp; SatName (optional) &nbsp;&nbsp;&nbsp;&nbsp; DataUsed (optional) &nbsp;&nbsp;&nbsp;&nbsp; CoordSys (optional) &nbsp;&nbsp;&nbsp;&nbsp; OrbType (optional) &nbsp;&nbsp;&nbsp;&nbsp; AgencyName (optional)

---

## Author

Sergio De Florio

---

## License
SpOCK is free software: you can redistribute the original source code and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3. The use of the third parties libraries used by SpOCK (SPICE, NRLMSIS-2.1, etc.) is instead subjected to the specific licenses of these libraries.

---

## Disclaimer
This software is the product of a free-time-project carried on by a single person and is provided 'as is' without warranty of any kind. There can be no warranty that:

* The software will meet your requirements
* The software will be uninterrupted, timely, secure or error-free
* The results that may be obtained from the use of the software will be effective, accurate or reliable
* The quality of the software will meet your expectations
* Any errors in the software will be corrected

The author will try to maintain it and fix any found or reported bug in a resonable time but he cannot guarantee it.

Software and its documentation:

* Could include technical or other mistakes, inaccuracies or typographical errors
* May be out of date, and the author makes no commitment to update such materials
