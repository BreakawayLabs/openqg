# OpenQG Quickstart Guide

The following instructions are designed to get you up and running with a simple example experiment. For more details on running the model please consult the user guide.

## Dependencies

OpenQG is targeted at Linux based, desktop class machines. The following packages are required to run on [Ubuntu 12.04 LTS](http://releases.ubuntu.com/precise/). 

    sudo apt-get install git make gfortran python libfftw3-dev libnetcdf-dev liblapack-dev netcdf-bin

Users of other systems should ensure equivalent packages are installed.

## Download

The source code for OpenQG is hosted on [github](https://github.com/BreakawayLabs/openqg). To download the latest development version, run:

    git clone git@github.com:BreakawayLabs/.git

For more details on how to download different releases, or to download as a zip/tar file, check the [downloads](http://openqg.breakawaylabs.com.au/web/downloads) page.

## Compile

To compile the model with OpenMP enabled (recommended), run

    make -C src

To compile the model in single processor mode, with OpenMP disabled, run

    make -C src single

## Run

The model should always be run via the `run_model.py` script, which ensures all input and output files are managed correctly. To get started running a predefined North Atlantic double gyre model for 10 days, run

    python run_model.py -x src/openqg -o output -e dg_fast

If you have compiled in single processor mode, run

    python run_model.py -x src/openqg-single -o output -e dg_fast

This will run the model and place the output in `/output`.

## What Next?

To find out more about configuring your own experiments and analysing results, please consult the user guide.
