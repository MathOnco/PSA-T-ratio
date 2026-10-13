# PSA-T-ratio
Contains data and code to reproduce biomarker analysis study. 

* Data is included in the data folder.
  1) Table tab contains patient-specific characteristics. Initial PSA values with + indicate a difference with previously published data (see data from each patient for full data set). Initial testosterone values with * indicate values from outside labs taken before onset of trial.
  2) Model1 tab contains patient-specific best parameter fits and parameter ranges to exponential growth and decay functions described in the Methods section.
  3) Model2 tab contains baseline PSA, T, and patient-specific best parameter fits for ODE model.
  4) Each patient has a tab with their data. Days is given from start of first ADT. Cycle 1 starts when first ADT starts and ends when treatment break ends and next cycle of treatment starts. PSA is given in ng/mL and testosterone is given in ng/dL. Notes indicate any extra notes that were given and may contain any time PSA or testosterone were below/above a threshold. ADT is the intended dose of androgen deprivation therapy (a 12 week dose is intended to go for 12 weeks even if the effects last longer/shorter). Abi is on=1 or off=0. ARSI is on=1 or off=0.

* Jupyter notebooks can be used to reproduce the results and figures. Please use the requirements.txt file to load relevant packages.

  1) data_analysis.ipynb contains statistical analysis, multivariate analysis, and other plots.

  2) modelRun.ipynb will execute the fitting algorithm for the mathematical model over the cohort and plot the results.


* Clinical protocol is given in 19367.protocolv4.3.2022-08-23.pdf.
