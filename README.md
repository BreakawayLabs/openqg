# Q-GCM - Quasigeostropic Coupled Model

Q-GCM is an idealised ocean-atmosphere model suitable for simulations of mid-latitude physical processes.
Q-GCM is suitable for use by both researchers and students.
It can be run on a standard desktop machine in either single- or multi-threaded mode.
Parallelisation is provided by OpenMP.

## Super quick-start

To get up and running, try the following commands:

    make -C src
    python run_model.py -x src/q-gcm -o output -e dg_fast

Congratuations, you just ran a North Atlantic double gyre experiment for 10 days using an OpenMP enabled version of Q-GCM.
The output can be found in `output/`.
For a slightly gentler introduction see the quickstart guide:

[http://qgcm.breakawaylabs.com.au/web/docs/project/quickstart](http://qgcm.breakawaylabs.com.au/web/docs/project/quickstart)

## Dependencies

Q-GCM links against LAPACK, NetCDF and FFTW3 and uses `gfortran` as the compiler.

## Documentation

Full documentation can be found at:

[http://qgcm.breakawaylabs.com.au/web/docs](http://qgcm.breakawaylabs.com.au/web/docs)

## Licence

Q-GCM is free software licenced under the standard GPLv3 licence.
