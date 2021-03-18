# iav-splash_src
Source code related to the IAV SPLASH Project, a collaboration with Freiburg.
The original code was provided by Celia. Untill now we implemented a prototype for elementwise variance calculation across NumPy arrays.

## `variance_by_element.py` 
### Take-home messages
- Due to the way we implemented the 'degeneration' of arrays (by adding an `int` between 0 and 9) the variance distribution is flat, as seen on the heatmap.
- Implementing different 'degeneration' functions might change the variance distribution.

### known issues
- The histogram plotting is not yet function