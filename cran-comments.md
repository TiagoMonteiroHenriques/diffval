Dear CRAN maintainers,

Thank you for your message and comments.
I now resubmit addressing all raised issues.

Sincerely,

Tiago Monteiro-Henriques

## Resubmission

This is a resubmission where:

* I included an explanation of the acronym TDV in the description text.

* The CRAN URL that was not in canonical form was replaced in the README files.

* I removed \dontrun{} from the example in the `explore_tabulation()` function,
  and added an if clause to guarantee that R is being used interactively. In
  function `optim_tdv_gurobi_k_2()` example, I kept the \dontrun{}, as that
  function calls Gurobi, which is an additional software that must be installed
  beforehand.

* Considering the NOTE about package 'gurobi' (suggested but not available for
  checking):
  Package 'gurobi' provides an interface to the Gurobi optimization software.
  DESCRIPTION file contains instructions for the installation of Gurobi,
  as well as for the installation of the 'gurobi' package.
  Following the practice found in package 'prioritizr', I had now included
  such instructions also in the documentation of the function
  `optim_tdv_gurobi_k_2()`, which is the only function from 'diffval' package
  that relies on Gurobi optimization software (see `?optim_tdv_gurobi_k_2`). The
  use of this software is relevant for the users, as it is the only way to
  guarantee the optimality of a solution for the complex problem being
  addressed.
  
## Test environments

* Local OS X 12.6 install, R 4.2.2
* R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R-hub Fedora Linux, R-devel, clang, gfortran
* R-hub Windows Server 2022, R-devel, 64 bit
* Win-builder x86_64-w64-mingw32 (64-bit), R Under development (unstable)
  (2022-10-11 r83083 ucrt)

## R CMD check results

R-hub checks were done setting `_R_CHECK_FORCE_SUGGESTS_ to false`, as package
'gurobi' is not in mainstream repositories.

There were no ERRORs or WARNINGs in all test environments.

### Local OS X

There were no NOTEs:

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

### R-hub Ubuntu Linux

There were 3 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’

New submission

Suggests or Enhances not in mainstream repositories:
  gurobi
````
This is a first submission.
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gurobi’
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking examples ... [8s/29s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.404  0.032  12.909
optim_tdv_hill_climb 2.673  0.023   9.920
```
Optimization can be time consuming. This example is small in size.
This is not happening in Windows or OS X environments.

### R-hub Fedora Linux

There were 4 NOTEs:
```
* checking CRAN incoming feasibility ... [6s/32s] NOTE
Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’

New submission

Suggests or Enhances not in mainstream repositories:
  gurobi
```
This is a first submission.
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gurobi’
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking examples ... [8s/27s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.545  0.017  12.572
optim_tdv_hill_climb 2.631  0.011   9.552
```
Optimization can be time consuming. This example is small in size.
This is not happening in Windows or OS X environments.

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable
```
Related to test environment.

### R-hub Windows Server 2022

There were 4 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Tiago Monteiro-Henriques <tmh.dev@icloud.com>'

New submission

Suggests or Enhances not in mainstream repositories:
  gurobi
```
This is a first submission.
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'gurobi'
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking HTML version of manual ... NOTE
Skipping checking math rendering: package 'V8' unavailable
```
Related to test environment.

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
Related to test environment.
As stated in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this
seems a bug/crash in MiKTeX.

### Win-builder x86_64-w64-mingw32 (64-bit)

There were 2 NOTEs:

```
* checking CRAN incoming feasibility ... [14s] NOTE
Maintainer: 'Tiago Monteiro-Henriques <tmh.dev@icloud.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  DiffVal (14:25)
  Gurobi (20:26)
  TDV (15:55)
  geobotanical (17:48)
  phytosociological (13:37)

Suggests or Enhances not in mainstream repositories:
  gurobi

Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation_guide.html
    From: README.md
    Status: 200
    Message: OK
    CRAN URL not in canonical form
  The canonical URL of the CRAN page for a package is 
    https://CRAN.R-project.org/package=pkgname
```
This is a first submission.
Spell check was performed. These terms are correct.
DESCRIPTION has full instructions for 'gurobi' package installation.
The URL is a vignette of package 'prioritizr' with instructions for 'gurobi'
package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'gurobi'
```
DESCRIPTION has full instructions for 'gurobi' package installation.

## Downstream dependencies

There are currently no downstream dependencies for this package.
