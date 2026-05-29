# impower133

Reproduce the IMpower133 Clinical Trial Results

## Overview

The impower133 package provides functions to reproduce the results of the
IMpower133 clinical trial (Horn et al., 2018, NEJM), a phase III randomized
controlled trial evaluating atezolizumab plus chemotherapy versus placebo plus
chemotherapy in extensive-stage small-cell lung cancer.

This package was developed as a course project for the Biostatistics course
at East China Normal University.

## Installation

install.packages("impower133")

Or install the development version from GitHub:

remotes::install_github("fanfande131/impower133")

## Functions

rtriang() - Triangular distribution random number generator

simulate_impower133() - Simulate complete ITT population with baseline, efficacy, and safety data

reconstruct_survival() - Reconstruct individual patient data from published K-M curves

analyze_survival() - Analyze survival data and print key statistics

plot_survival() - Plot K-M survival curves for OS or PFS

make_table1() - Generate Table 1: Baseline characteristics

make_table2() - Generate Table 2: Response rate and disease progression

make_table3() - Generate Table 3: Adverse events

plot_figure2c() - Generate Figure 2C: Subgroup forest plot

## Usage

See `inst/examples/test_impower133.R` for a complete example. 
Simulated data, output figures (tables and plots), and the reconstructed IPD 
data are also available in `inst/examples/`.

## Data Sources

The individual patient data was reconstructed from published Kaplan-Meier
curves using the IPDfromKM web application:

https://biostatistics.mdanderson.org/shinyapps/IPDfromKM/

Method reference:
Liu N, Zhou Y, Lee JJ. IPDfromKM: reconstruct individual patient data
from published Kaplan-Meier survival curves. BMC Med Res Methodol.
2021;21(1):111.

## Trial Reference

Horn L, Mansfield AS, Szczesna A, et al. First-Line Atezolizumab plus
Chemotherapy in Extensive-Stage Small-Cell Lung Cancer. N Engl J Med.
2018;379:2220-2229. doi:10.1056/NEJMoa1809064

## License

MIT
