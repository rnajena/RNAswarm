# iav-splash_src
Source code related to the IAV SPLASH Project, a collaboration with Freiburg.
The original code was provided by Celia. Untill now we implemented a prototype for elementwise variance calculation across NumPy arrays.

## `variance_by_element.py` 
### Take-home messages
- Due to the way we implemented the 'degeneration' of arrays (by adding an `int` between 0 and 9) the variance distribution is flat, as seen on the heatmap.
- Implementing different 'degeneration' functions might change the variance distribution.

### known issues
- The histogram plotting is not yet function

## `support`
This repository has some scripts from the ['original' splash analysis](https://github.com/CSB5/splash) (which is not being maintained and uses Python2) script from [Nagarajan's laboratory](http://csb5.github.io/). It also includes new source code (probably from Celia). They removed two scripts from the original repo: `filter_chimeras.py` and `pickJunctionReads.awk`.

## `main.py`
As the name states, this is the driver script of the pipeline, it calls all the functions and the `run.conf` file is used to configure the run.

## `run.conf`
Configuration file for each run.




