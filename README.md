# Quasi-Love Wave Detections with SAC
- Codes for detecting Quasi-Love waves, featuring a Bash/SAC/GMT option and a Python/ObsPy/SAC/PyGMT option.
- Written by Tom Merry
- Based on work by: Eakin (2021), Quasi-Love wave scattering reveals tectonic history of Australia and its margins reflected by mantle anisotropy, _Communications Earth & Environment_, 2, 210 https://doi.org/10.1038/s43247-021-00276-7 , original code found at https://github.com/SeismoCaro/Quasi-Love-Wave-Detections
- Used for: Merry & Eakin (2024), Insights on the African Upper Mantle From Quasi-Love Wave Scattering, _Geochemistry, Geophysics, Geosystems_, 25, https://doi.org/10.1029/2023gc011385 
- Code reproduces Figure ???


**Overview of Files**

- *QuasiLove*
  - This folder contains codes to be imported in for use in detecting and locating Quasi-Love waves
  - quasilove_fns.sh: bash functions that run the QL analysis. The main function is qlmain

More info coming soon

**Steps to run**

Coming soon

**Requirements**

**SAC** should be installed and callable by the command ```sac``` (but the path could be defined in the quasilove_fns.sh file). This is written for the standard (Linux) distribution of SAC (note that this can also be installed on a Mac!), which means we expect the SAC binary files to be little-endian. Not sure what minimum version of SAC is needed but any reasonably recent one should do. (Maths in sac macros here is written as, e.g., ```a * b``` rather than ```mul a b```; this was because ```mul``` and ```div``` were failing for me on one machine. Hopefully this way is more portable?)

**GMT** version 6+. NB: some versions of ghostscript fail to preserve semi-transparencies from GMT, which would be very annoying for the example figure produced.

**Python**: Only built-in modules (math, sys) are used if using the bash version of this code. I've only tried it with Python3 but perhaps Python2 works too.
