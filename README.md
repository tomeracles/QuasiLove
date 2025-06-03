# Quasi-Love Wave Detections with SAC
- Codes for detecting Quasi-Love waves. The main working method uses Bash/SAC/GMT (with a tiny bit of Python). There is also a Python/ObsPy/PyGMT option which may be working to some extent but is effectively an abandoned work in progress, so may or may not be of interest to anyone.
- Written by Tom Merry while at the Research School of Earth Sciences, the Australian National University (ANU), in 2023.
- Based on work by: Eakin (2021), Quasi-Love wave scattering reveals tectonic history of Australia and its margins reflected by mantle anisotropy, _Communications Earth & Environment_, 2, 210 https://doi.org/10.1038/s43247-021-00276-7 , original (Matlab) code found at https://github.com/SeismoCaro/Quasi-Love-Wave-Detections
- Used for: Merry & Eakin (2024), Insights on the African Upper Mantle From Quasi-Love Wave Scattering, _Geochemistry, Geophysics, Geosystems_, 25, https://doi.org/10.1029/2023gc011385 
- This is somewhat of a minor update to the Eakin (2021) method in that it uses the Hilbert transform to compare the vertical and radial motions, plots particle motion and quantifies ellipticity (during the process of TM translating the code to Python, and then, perhaps ill-advisedly, to SAC). See Merry & Eakin (2024) for details.
- See below for instructions to reproduce Figure 3 of Merry & Eakin (2024).

***

**Overview of directories**

- *QuasiLove*
  - This folder contains codes to be imported in for use in detecting and locating Quasi-Love waves.
  - quasilove_fns.sh: bash functions that run the QL analysis. The main function is qlmain. Most of the maths is done directly in SAC, which is convenient and quick, but the language of the SAC macros is a little obscure and not terribly well documented. 
  - QL_detection.py: python functions (using obspy, pygmt) to process and plot the QL analysis. This is not finished I think —- it looks like it started to introduce the Hilbert transform but the derivative is still in there. I haven't spent the time to check how close this is to the final method (I have to imagine it would not be that hard to make them identical).

- *scripts*
  - *run_example_ql.sh*: This runs the ql analysis defined in QuasiLove/quasilove_fns.sh on the example data located in exampledata; it should recreate Fig. 3 of Merry & Eakin, 2024 (minus a couple of annotations). You should find the results and figure in exampledata/results.
  - *bulk_ql_analysis.sh*: This is a simple script (that will however definitely only work on a Mac) that goes through doing the same analysis as above on a whole set of data (see comments at the top of the script for more details). Each time it does the analysis, it throws up the figure to the screen and asks for a rating, which is saved in a logfile.
  - *bulk_ql_noview.sh*: Same as above but doesn't stop to let you view and rate the results – this is a much better use of time, i.e. let the analysis run through, then subsequently view the results with...
  - *view_ql_results*.sh. As with ```bulk_ql_analysis.sh```, this is written for a Mac — it does the viewing and rating bit of that script.
  - *download_waveforms.py* Script using obspy to download the waveforms that I used in SAC format. Includes some hard-coded info about the event catalogue, networks and data clients. Not well tested. See the comments and/or run ```./scripts/download_waveforms.py --help``` for usage info.

- *exampledata*
  - Some example data from station II.SUR that is used in Fig. 3 of Merry & Eakin, 2024.

***

**Notes for running the example scripts**

These are designed to be run from this parent directory, e.g., run 

```./scripts/run_example_ql.sh exampledata/waveforms/20150512211258_II.SUR.LHZ```

Note that you can choose any of the Z, R or T channels as the input for this script, it doesn't make a difference.

***

**Requirements**

-  **SAC** should be installed and callable by the command ```sac``` (but the path could be defined in the quasilove_fns.sh file if not). This is written for the standard (Linux) distribution of SAC (note that this can also be installed on a Mac!), which means we expect the SAC binary files to be little-endian. Not sure what minimum version of SAC is needed but any reasonably recent one should do. (Maths in sac macros here is written as, e.g., ```a * b``` rather than ```mul a b```; this was because ```mul``` and ```div``` were failing for me on one machine. Hopefully this way is more portable?)

  For more information on SAC macros, which are used for a lot of the calculations in these codes, see _Seismic Analysis Code: A Primer and User's Guide_ by Helffrich, Wookey and Bastow (2013).

-  **GMT** version 6+. NB: some versions of ghostscript fail to preserve semi-transparencies from GMT, which would be very annoying for the example figure produced.

-  **Python**: Only built-in modules (math, sys) are used if using the bash version of this code. I've only tried it with Python3 but perhaps Python2 works too.

-  For the download script, *obspy* is needed (as well as certain things obspy depends on, such as numpy). Any reasonable recent version of obspy should do.

For the (work-in-progress) python-based QL analysis, obspy and pygmt are needed. See https://www.pygmt.org/latest/ for PyGMT --- this is changing all the time, and I was probably using approximately v0.9 when writing, but if I were to pick this up again I'd use a more recent version...
