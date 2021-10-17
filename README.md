# SLOSH
SLOSH is a pipeline designed to analyze SPLASH (**S**equencing of **P**soralen crosslinked, **L**igated, **A**nd **S**elected **H**ybrids) datasets.

According to [Wikipedia](https://en.wikipedia.org/wiki/Slosh_dynamics): "In fluid dynamics, slosh refers to the movement of liquid inside another object (which is, typically, also undergoing motion)."

## Dependencies (dev state)
### iav-splash conda env

### conda envs

#### Channel priority

```Bash
channels:
  - conda-forge
  - bioconda
  - r
  - defaults
```

#### iav_splash

This environment has all dependencies needed to run slosh, but doesn't include development dependencies (black, ipython, etc)

```Bash
conda create --name iav_splash bioconductor-deseq2 r-base numpy seaborn scipy pandas circos
```

#### iav_splash-dev

This environment has all dependencies needed to run slosh, including development dependencies (ipykernel)

```Bash
conda create --name iav_splash-dev bioconductor-deseq2 r-base numpy seaborn scipy pandas circos ipykernel black nbconvert
```

#### iav_splash-r

This environment only has dependencies related to R code

```Bash
conda create --name iav_splash-r bioconductor-deseq2 r-base
```
