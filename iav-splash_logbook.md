# IAV SPLASH Project Logbook
## 2021.04.12 - Exploring fontionality of vRNAsite
### Context
`vRNAsite` is a software developed by an alumni of the group, Daniel Desirò. We are planning on using it to 'annotate' RNA structures in our samples and use those annotations to count SPLASH reads.

## 2021.04.16 - Finish refactoring and regex implementation
Tags: #regex
## Goal
After playing a little with Python's re library, I will finish implementing a regular expression (regex) to get the segment pairs from the np.array's filenames, and refactoring the code according to Kevin's suggestions.
## Useful resources
- [Regular Expression HOWTO](https://docs.python.org/3/library/re.html)
- [`re` module documentation](https://docs.python.org/3/howto/regex.html#regex-howto)
- 
## Results
The refactoring of the code is finished, but the `check_heatmaps` is not yet working and has to be fixed.

## 2021.04.20 - Finish refactoring and regex implementation
### Quick log
The problem was not `check_heatmaps` function, it was just the fact that I was using the wrong filepaths in the file.

## 2021.04.21 - Implement functions to check magnitude of variances
### Goal
- Work on [issue #5](https://github.com/gabriellovate/iav-splash_src/issues/5)

## 2021.04.27 - Test vRNAsite and work on issue #5
### Goal
- Run vRNAsite for genomes provided by Celia
- Work on [issue #5](https://github.com/gabriellovate/iav-splash_src/issues/5)

## 2021.04.29 - Generate visualisations for presentation
### Goals
- Run vRNAsite for genomes provided by Celia
    - Command used: 
```bash
python3 vRNAsite.py -pfx \
/home/ru27wav/projects/gl_iav-splash_freiburg/results/202104/20210429/vRNAsite/plots \
-fsa \
/home/ru27wav/projects/gl_iav-splash_freiburg/results/202104/20210429/vRNAsite/SC35M.fasta \
-ovr -rev -cmp -rvp -mat
```
- Generate heatmaps for some interactions