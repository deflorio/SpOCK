# SpOCK
The Spacecraft Orbital Computations Kit (SpOCK) is an open source tool for spacecraft mission analysis and simulation. SpOCK allows the simulation of spacecrafts' hardware, orbital and attitude dynamics, and the computations of mission events (ground station contacts, payload data-takes and eclipses). The ambition of the author is to develop an operational tool (mission analysis, flight dynamics and operations) which includes only few but accurate simulations features and which can be easily further developed and customized. For this reason the software is developed with a (hopefully) clear structure with elements which can be easily read and reproduced and without a graphical user interface. For detailed information about installation, configuration and use please read the user guide.

---

## Features

*Spacecraft hardware*

Model of 6 spacecraft faces, 3 solar panels, 1 Sun camera, 1 Earth camera, 3 magnetometers, 6 coarse Sun sensors, 3-axis rate sensor, 3 reaction wheels, 3 magnetorquers, 2 types of orbit control propulsion systems.

*Dynamics*

Precise orbit and attitude propagation, mission events computation, atmospheric drag and solar radiation pressure with panels model (dependent on attitude) or reference area model, attitude and orbital maneuvers with different type of actuators. The simplified general perturbations (SGP4) propagator used with two-line mean element (TLE) sets is also implemented.

*Orbit environment models*

Gravity field (all GGM0\*, EIGEN\* and models with .gfc files format of this type), atmospheric density (JB2008, NRLMSISE-00), Earth's magnetic field (IGRF13, WMM2020), solar radiation pressure, third body (JPL planet ephemerides), gravity gradient.

*Customizability*

Easy implementation of new spacecraft hardware, ready-to-use insertion of attitude and orbit control modules, set up for hardware-in-the-loop simulations.

---

## System Requirements

Linux Ubuntu 19 or higher

GNU make (sudo apt-get install build-essential)

gcc >= 7 supporting -std=c++17 (sudo apt-get install build-essential)

Libraries usually not included by default in Linux Ubuntu distribution: libboost-all-dev, libxerces-c-dev, xsdcxx, gfortran, freeglut3, freeglut3-dev, mesa-utils, libsdl2\*, libsoil\*, doxygen, graphviz. These libraries (1125 MB of addtional disk space) will be installed automatically during the installation process.

---

## Installation

*In root directory:*

make install_libs (first compilation or after OS upgrade)

make install_atmo (first compilation or after OS upgrade)

make install_mag (first compilation or after OS upgrade)

make orbit

make sgp4

make attitude

make events

All the executable will be generated in the 'bin' folder. Command 'make clean' will delete all executables in 'bin' folder and folder 'obj'.

Download file de440.bsp from SPICE library kernels ftp (ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/) and put it into folder data/cspice (file de440.bsp is too big for the repository)

*Remarks:*

Installation and run were tested successfully on different real and virtual machines (VMware on Windows 10) running Linux Ubuntu 19.04 or higher. It is likely but not guaranteed that the installation on Ubuntu < 19 or other Linux distributions could work

Coomands make install_atmo and install_mag are used to compile external Fortran libraries and could give some Fortran compilation warnings.

---

## Configuration

Put correct inputs/outputs paths in input files simparam_sample.xml and eventsparam_sample.xml as explained in detail in SPOCK_UserGuide.pdf.

---

## Run

*In linux terminal:*

./bin/OrbitPropagator input/SimulationParameters/simparam_sample.xml

./bin/SGP4Propagator input/SimulationParameters/simparam_sample.xml

./bin/AttitudePropagator input/SimulationParameters/simparam_sample.xml

./bin/EventsComputation input/SimulationParameters/eventsparam_sample.xml

---

## Author

Sergio De Florio

---

## License
SpOCK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3.

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
