# Q-GCM User Guide

This document outlines everything you need to know as a user of Q-GCM.
If you have problems or questions which are not covered here, please address them to the Q-GCM [mailing list](https://groups.google.com/forum/#!forum/q-gcm-users).

## Downloading Q-GCM

### Installing git

Q-GCM uses git for version control, with an official central repository hosted at [GitHub](https://github.com/BreakawayLabs/q-gcm).
The recommended way of accessing the Q-GCM source is to clone the official repository using git.
Git can be installed on Debian/Ubuntu systems with the following command.

    sudo apt-get install git

Installation instructions for other platforms can be found on the [git website](http://git-scm.com/)

### Downloading Stable Releases

#### Downloading as Zip/Tar File

### Downloading Development Branch

To get the latest development version of Q-GCM, run the following command:

    git clone https://github.com/BreakawayLabs/q-gcm.git

This will download the development branch into a local repository in the directory `q-gcm`.
While all efforts are made to ensure the development branch is stable, it is possible that there will be bugs ranging from the minor to the catestrophic.
As such it is not recommended to user the development branch for research unless there is a specfic new feature which is needed.

If you find a problem please check the [task tracker](http://qgcm.breakawaylabs.com.au/youtrack) to see if it is already known about.
If you cannot find an existing issue which covers your problem, please create a new issue so that the problem can be fixed ASAP.
See the [Reporting Bugs](#bugs) section below.

### Accessing Annexed Files

## Dependencies

## Compiling Q-GCM

Q-GCM is written in [Fortran 90](http://en.wikipedia.org/wiki/Fortran#Fortran_90) and distributed in source only form. As such, users must compile the model before running it.

By default Q-GCM uses the [`gfortran`](http://gcc.gnu.org/wiki/GFortran) compiler. Version 4.6.3 has been successfully tested, however any recent version should work fine.

To compile the model, run the following commands:

    cd src
    make

This will create the executable file `src/q-gcm`.

## Running Predefined Example Experiments

Q-GCM comes with a number of predefined example experiments.
These can be sued as a starting point for customising your own experiments.
It is advised to run a number of these examples to ensure the model is running correctly on your system.

## Experiment Directory Structure

## Running Custom Experiments

## Configuring Q-GCM

## Diagnostics

## Advanced Compiler Options

## <a id="bugs"></a>Reporting Bugs