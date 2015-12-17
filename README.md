# Heterogenous firing response of mice layer V pyramidal neurons in the fluctuation-driven regime

- Pre-analyzed electrophysiological data (the intracellular recordings
  in current clamp, ~10 GBytes over all cells, have been reduced to a
  set of episodes containing 1) the desired fluctuations properties,
  2) the monitored quantities and 3) the output firing rate, see the
  manuscript for details)

- python code for the data analysis and theoretical modeling

- manuscript in org mode

## Requirements

- doit : the automation tool in python, see http://pydoit.org/ 

### For data analysis and theoretical modeling

- scientific distribution of python (e.g. anaconda, visit
  https://www.continuum.io/downloads)
  needed packages :
  - numpy
  - scipy
  - numba
  - matplotlib

### To generate the pdf manuscript from the =Org-Mode= source file and the svg figures

  - pdflatex
  - emacs > 24.5.0, get it from https://www.gnu.org/software/emacs/#Obtaining
  
## Perform numerical simulations

Go to the 'numerical simulations/' folder and run 'doit'

see the individual commands executed to 

## Analyze data

Go to the 'analysis/' folder and run 'doit'

