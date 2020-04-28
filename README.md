# SpOCK
The Spacecraft Orbital Computations Kit (SpOCK) is an open source tool for spacecraft mission analysis and simulation. SpOCK allows the simulation of spacecrafts' hardware, orbital and attitude dynamics, and the computations of mission events (ground station contacts, payload data-takes and eclipses). The ambition of the author is to develop an operational tool (mission analysis, flight dynamics and operations) which includes only few but accurate simulations features and which can be easily further developed and customized. For this reason the software is developed with a (hopefully) clear structure with elements which can be easily read and reproduced and without a graphical user interface. For detailed information about installation, configuration and use please read the user guide.

---

## System Requirements

Linux Ubuntu 18 or higher

GNU make (sudo apt-get install build-essential)

gcc $>=$ 7 supporting -std=c++17 (sudo apt-get install g++)

Libraries usually not included by default in Linux Ubuntu distribution: libboost-all-dev, libxerces-c-dev, xsdcxx, gfortran, freeglut3, freeglut3-dev, mesa-utils, libsdl2*, libsoil*, doxygen, graphviz. These libraries will be installed automatically during the installation process.

For the installation of these libraries, 1125 MB of addtional disk space will be used.

---

## Installation

In root directory:

make install_libs (if it is the first compilation)

make install_atmo

make install_mag

make orbit

make attitude

make events

All the executable will be generated in the 'bin' folder. Command 'make clean' will delete all executables in 'bin' folder and folder 'obj'.

Download file de438.bsp from SPICE library kernels ftp (ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/) and put it into folder data/cspice (file de438.bsp is too big for the repository)

---

## Configuration

Put correct inputs/outputs paths in input files simparam_sample.xml and eventsparam_sample.xml

---

## Run

In linux terminal:

./bin/OrbitPropagator input/SimulationParameters/simparam_sample.xml

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
