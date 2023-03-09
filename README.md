
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diffval

<!-- badges: start -->
<!-- badges: end -->

Find, visualize and explore patterns of differential taxa in vegetation
data (namely in a phytosociological table), using the Differential Value
(DiffVal). DiffVal captures the degree of exclusiveness of a taxon to
each of the different clusters of relevés of a sorted phytosociological
table. The patterns of differential taxa are searched through
mathematical optimization algorithms, resorting the table while
maximizing the sum of the DiffVal for all the taxa in the table, i.e.,
maximizing the Total Differential Value (TotDiffVal or TDV). Ultimately,
TDV optimization aims at obtaining classifications of vegetation data
based on differential taxa, as in the traditional phytosociological (or
geobotanical) approach. The Gurobi optimizer, as well as the R package
‘gurobi’, can be installed from
<https://www.gurobi.com/products/gurobi-optimizer/>. The vignette
“Gurobi Installation Guide”, from package ‘prioritizr’, is useful and
can be found here:
<https://prioritizr.net/articles/gurobi_installation_guide.html>.

## Acknowledgements

TMH was funded by the European Social Fund (POCH and NORTE 2020) and by
National Funds (MCTES), through a FCT – Fundação para a Ciência e a
Tecnologia (Portuguese Foundation for Science and Technology)
postdoctoral fellowship (SFRH/BPD/115057/2016), as well as by National
Funds, through the same foundation, under the project UIDB/04033/2020
(CITAB - Centre for the Research and Technology of Agro-Environmental
and Biological Sciences). JOC was financially supported by FCT through
the projects UIDB/MAT/00297/2020, UIDP/MAT/00297/2020 (Centro de
Matemática e Aplicações).

## Related articles

Portela-Pereira E, Monteiro-Henriques T, Casas C, Forner N,
Garcia-Cabral I, Fonseca JP & Neto C 2021. Teixedos no noroeste da
Península Ibérica. Finisterra LVI(117): 127–150. DOI:
10.18055/FINIS18102
<https://revistas.rcaap.pt/finisterra/article/view/18102>

## Installation

You can install the development version of the package from GitLab.

``` r
remotes::install_gitlab("point-veg/diffval")
```
