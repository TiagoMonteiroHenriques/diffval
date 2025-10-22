
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
maximizing the Total Differential Value (TotDiffVal/TDV). Ultimately,
TDV-optimization aims at obtaining classifications of vegetation data
based on differential taxa, as in the traditional phytosociological (or
geobotanical) approach (Monteiro-Henriques 2025,
<https://doi.org/10.3897/VCS.140466>). The Gurobi optimizer, as well as
the R package ‘gurobi’, can be installed from
<https://www.gurobi.com/products/gurobi-optimizer/>. The vignette
“Gurobi Installation Guide”, from package ‘prioritizr’, is useful and
can be found here:
<https://prioritizr.net/articles/gurobi_installation_guide.html>.

## How to use

The package provides several optimization approaches for maximizing TDV.
Finding a partition into $k$ subsets that maximizes TDV is NP-hard for
any fixed $k$. For the case $k=2$, integer linear programming
formulations are implemented in `optim_tdv_gurobi_k_2()`, allowing
optimal solutions to be obtained for moderately sized data sets.

For larger data sets and for $k>2$, (meta)heuristic algorithms are
provided: a *greedy randomized adaptive search procedure* (GRASP), a
*simulated annealing* approach, and a *hill climbing* approach.
Good-quality partitions can often be achieved by running the simulated
annealing algorithm initialized with solutions obtained from GRASP (this
is implemented in the function `optim_tdv_simul_anne()`). The function
`partition_tdv_grasp()` implements GRASP on its own.

For the hill climbing approach, it is both possible and advisable to
combine a first stage of *stochastic hill climbing* followed by a second
stage of *greedy hill climbing* (see Monteiro-Henriques 2025). This is
implemented in the function `optim_tdv_hill_climb()`. Multiple runs of
this approach explore the solution space broadly, allowing not only
global optimization but also the recording of local optima, which can be
relevant for understanding the data set’s structure (see
Monteiro-Henriques 2025).

For very large instances, the above (meta)heuristics may become
computationally expensive. In such cases, an alternative is to use the
function `partition_tdv_grdtp()`, which implements a simplified version
of the greedy algorithm. This function performs faster than
`optim_tdv_simul_anne()` or `optim_tdv_hill_climb()`, albeit at the
expense of the quality of the solutions produced.

## Acknowledgements

TMH was partially funded by the European Social Fund (POCH and NORTE
2020) and by National Funds (MCTES), through a FCT – Fundação para a
Ciência e a Tecnologia (Portuguese Foundation for Science and
Technology) postdoctoral fellowship (SFRH/BPD/115057/2016), as well as
by National Funds, through the same foundation, under the project
UIDB/04033/2020 (CITAB - Centre for the Research and Technology of
Agro-Environmental and Biological Sciences). JOC was financially
supported by FCT through the projects UID/297/2025 and UID/PRR/297/2025
(Center for Mathematics and Applications – NOVA Math).

## Related articles

Monteiro-Henriques T 2025. TDV-optimization: A novel numerical method
for phytosociological tabulation. Vegetation Classification and Survey
6: 99-127. DOI: <https://doi.org/10.3897/VCS.140466>

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
