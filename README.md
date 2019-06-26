# OpenLSTO
**Open Source Code for Level Set Topology Optimization**

OpenLSTO is an open-source software for level set based structural topology optimization, written in C++, developed and Maintained by researchers at UC San Diego and Cardiff University.

The topology optimization method deals with the problem of determining the optimal distribution inside a design domain in order to obtain the best structural performance. Level Set Method was originally developed as a mathematical tool for tracking the motion of interfaces. Its natural handling of topological changes coupled with a clear and smooth interface representation led to its use for structural topology optimization.

This first full release implements the M2DO lab’s level set topology optimization method to solve the problems of minimizing compliance (maximizing stiffness) under a volume constraint and minimizing stress under a volume constraint in 2D and 3D design domains. A hole nucleation algorithm for the compliance problem is also included. While OpenLSTO has been designed with ease of installation and use in mind, some third-party packages have been used in order for the code to be able to handle more complex problems. A light version of the code without third-party packages, and with less capabilities, is available [here](https://github.com/M2DOLab/OpenLSTO-lite)

The level-set module of this software, M2DO_LSM, was adapted from code written by [Lester Hedges](https://github.com/lohedges/slsm), which can be found [here](https://github.com/lohedges/slsm).

## Installation
-  Download

To clone the repository to the local machine,(git client) 
```git clone https://github.com/M2DOLab/OpenLSTO.git```

The Zip file is also available with the following static link:
[Download link](http://m2do.ucsd.edu/static/zip/OpenLSTO-v1.0.zip)

- CLI environment

In general, all OpenLSTO execution occurs via command line arguments within a terminal. For Unix/Linux or Mac OS X users, the native terminal applications are needed. 

- Execution

Go to the relevant folder under projects to run a specific module, e.g to compliance minimization:
> cd projects/compliance

Compile & Build: 
> make main

Execution
> ./bin/a.out

Clean existing compilation files:
> make clean

clean existing object files:
> make clean_obj

clean existing output files:
> make clean_results

- External dependency

Users of OpenLSTO need a data visualization tool to post-process solution files. The software currently supports .vtk output format natively read by [*ParaView*](https://www.paraview.org/).

## Licensing
OpenLSTO is available for download under the [Apache V. 2.0 license](http://www.apache.org/licenses/LICENSE-2.0). Please refer to the License page for terms and conditions.

## Contributors

[Prof. Hyunsun Alicia Kim](http://m2do.ucsd.edu) (UCSD, Cardiff U.)

Dr. Sandilya Kambampati (UCSD)

[Dr. Lester Hedges](http://lesterhedges.net)

Dr. Zongliang Du (UCSD)

Dr. Renato Picelli (Cardiff U.)

Dr. Scott Townsend (Cardiff U.)

Dr. Xiao-Yi Zhou (Cardiff U.)

Dr. Hayoung Chung (UCSD)

Dr. Eliana Bortot (UCSD)

Ms. Carolina Jauregui (UCSD)

Mr. Andreas Neofytou (Cardiff U.)

Mr. Douglas de Aquino Castro (UCSD)
