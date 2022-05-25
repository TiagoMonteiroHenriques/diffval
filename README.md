
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diffval

<!-- badges: start -->
<!-- badges: end -->

Find, visualize and explore patterns of differential taxa in a
phytosociological table (vegetation data), using the Differential Value
(DiffVal) index. The patterns are searched through mathematical
optimization algorithms, maximizing the sum of the DiffVal index for all
the taxa in the table, i.e. the Total Differential Value (TotDiffVal or
TDV). Ultimately, the TDV optimization aims at obtaining classifications
of vegetation data based on differential taxa, i.e. taxa that are
exclusive of a certain group (or exclusive of some of the groups), as in
the traditional phytosociological/geobotanical approach.

## Acknowledgments

T.M.H. was funded by the European Social Fund (POCH and NORTE 2020) and
by National Funds (MCTES), through a FCT – Fundação para a Ciência a
Tecnologia postdoctoral fellowship (SFRH/BPD/115057/2016), as well as by
National Funds, through the same foundation, under the project
UIDB/04033/2020.

## Installation

You can install the package from GitLab.

``` r
devtools::install_git("https://gitlab.com/point-veg/diffval")
```
