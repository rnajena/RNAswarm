# Legacy codebase
This set of scripts were provided by Celia Jakob, and includes code from Redmond Smith and [code from the original SPLASH publication](https://github.com/CSB5/splash). Code from the original paper is under an MIT license, while the rest was provided without clear licensing.

## Code structure

### [job.sh](job.sh)
This is a SLURM job file that the user calls when using the script in a SLURM-managed system. It basically activates a conda environment that has all the dependencies and then runs [main.py](main.py) with [run.conf](run.conf) as an argument to do the actual work. The command is the following:

`python main.py run.conf`

### [run.conf](run.conf)
A basic configuration file that is used by [main.py](main.py) to specify parameters of the run.

### [main.py](main.py)
This is the main script of the codebase, it imports a bunch of packages, specifies a python `class` called `Configuration`, a `main` function and a bunch a functions that seem to be called by `main`.

#### Imported packages
- os
- datetime
- configparser
- argparse
- [pytoml](https://pypi.org/project/pytoml/): The project is no longer being actively maintained, the developers sugest using [toml](https://github.com/uiri/toml) instead.

#### Imported modules
- [from support: log](support/log.py)
- [from alignment: Aligner](alignment.py)
- [from support.helper_v1CJ: timeStamp](support/helper_v1CJ.py)
- [from chim_handler import ChimerHandler](chim_handler.py)

### [alignment.py](alignment.py)


### [chim_handler.py](chim_handler.py)


### [support/analyse_interactions.py](support/analyse_interactions.py)


### [support/chimeras2_rs.py](support/chimeras2_rs.py)


### [support/find_chimeras_rs.py](support/find_chimeras_rs.py)


### [support/helper_v1CJ.py](support/helper_v1CJ.py)


### [support/log.py](support/log.py)


### [support/RRseq_plottingCJ.py](support/RRseq_plottingCJ.py)


### [support/SPLASH_singCJ.py](support/SPLASH_singCJ.py)

