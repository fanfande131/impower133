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

remotes::install_github("your-username/impower133")

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

library(impower133)

# Simulate baseline data with OS survival outcomes
df <- simulate_impower133(
  seed = 21,
  os_path = "path/to/reconstructed_os.csv",
  method = "random",
  output_path = "simulated_data.csv"
)

# Generate tables
t1 <- make_table1(df)
t2 <- make_table2(df)
t3 <- make_table3(df)

# Save tables as HTML
gt::gtsave(t1, "table1.html")
gt::gtsave(t2, "table2.html")
gt::gtsave(t3, "table3.html")

# Plot OS K-M curve
plot_survival(
  "path/to/reconstructed_os.csv",
  type = "OS",
  output_path = "figure2a.png"
)

# Plot subgroup forest plot
plot_figure2c(df, output_path = "figure2c.png")

## Data Sources

The individual patient data (IPD) used in this package was reconstructed from
the published Kaplan-Meier curves of the IMpower133 trial. The reconstruction
process was performed using the IPDfromKM web application, available at:

https://biostatistics.mdanderson.org/shinyapps/IPDfromKM/

This tool implements the method described by Liu et al. (2021) to extract
coordinates from Kaplan-Meier plots and estimate the underlying patient-level
survival data.

## Dependencies

R (>= 4.0.0), survival, survminer, gt, forestplot, ggplot2, IPDfromKM

## License

MIT

## References

Horn L, Mansfield AS, Szczesna A, et al. First-Line Atezolizumab plus
Chemotherapy in Extensive-Stage Small-Cell Lung Cancer. N Engl J Med.
2018;379:2220-2229. doi:10.1056/NEJMoa1809064
