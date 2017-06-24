# OpenFOAMTrainingCombustion
**OpenFOAM Training: Combustion with OpenSMOKE++ (3-5 July 2017, Brussels)**

This repository contains the material (source code, kinetic mechanisms, documentation) adopted during the OpenFOAM training session on combustion, scheduled for **Tuesday, 4 July: Use of external libraries for chemistry - Getting started with OpenSMOKE++**

It is important that you read carefully the notes and instructions reported below, before the beginning of the training.  
If you already successfully configured your environment according to the instructions reported in the [OpenFOAMTrainingCombustionTest](https://github.com/acuoci/OpenFOAMTrainingCombustionTest) repository, you can skip what is reported in the following.

# Prerequisites
* OpenFOAM-dev or OpenFOAM-4.x (https://openfoam.org/)
* Boost C++ (http://www.boost.org/)

# Preparation of your environment
The OpenSMOKE++ library relies some third-party libraries for some important tasks:
* Eigen++ (for linear algebra operations)
* RapidXML (for XML file management)
* Boost C++ (for string manipulation, file system operations, etc.)

The Eigen++ and RapidXML libraries do not require precompilation, since they are entirely based on header files and they are provided together with the OpenSMOKE++ libraries. The Boost C++ libraries are not provided with OpenSMOKE++. Usually they are already available in your Linux distribution. If not, you need to install them.  

If the distribution you are using is Ubuntu or SuSE, you can follow procedure 1 reported below.  

If this is not the case, compilation from source code is required (see procedure 2).

## 1. Installation of Boost C++ (Ubuntu and SuSE)
For Ubuntu and SuSE distribution you can use one of the following commands:

**Ubuntu (versions 14.04 or above)**  
`sudo apt-get install libboost-all-dev`

**SuSE (OpenSuSE/SLES v12 or above, or Tumbleweed)**   
`sudo zypper install boost-devel`

## 2. Compilation and installation of Boost C++ from source code
If you cannot install the Boost C++ libraries using the commands reported above, you need to download the source code and compile it.

1. Download the source code from:  
https://dl.bintray.com/boostorg/release/1.64.0/source/

2. Follow the compilation and installation instructions reported at:  
http://www.boost.org/doc/libs/1_64_0/more/getting_started/unix-variants.html

After compilation and installation, set the environment variable pointing at the location where the compiled Boost C++ libraries have been installed (i.e. the folder containing the `include` and `lib` subfolders):

**bash or ksh**  
`export BOOST4OPENSMOKEPP=/path/to/boost`

**tcsh or csh**  
`setenv BOOST4OPENSMOKEPP /path/to/boost`

# Testing the environment
In order to check if the environment is correctly configured, we compile and run a simple application based on OpenSMOKE++ and OpenFOAM.
1. Go to the `Training/TestEnvironment` folder
2. Compile the source code by typing: `wmake`
3. If compilation succeeded, run the solver by typing: `environmentTest`  

If everything was done properly, you should have the following output on the screen:
```
...
Selecting ODE solver seulex
0.000000e+00
1.000000e-06
2.000000e-06
3.000000e-06
4.000000e-06
5.000000e-06
6.000000e-06
7.000000e-06
8.000000e-06
9.000000e-06
```

If you have any questions, please email me at alberto.cuoci(at)gmail.com (Subject: *OpenFOAMTrainingCombustionTest Issues*)
