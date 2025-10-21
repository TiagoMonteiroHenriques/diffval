Dear CRAN maintainers,

I submit diffval version 1.2.0, a minor release to the current CRAN version.

Sincerely,

Tiago Monteiro-Henriques

## General comment

* Considering the NOTE about package 'gurobi' (suggested but not available for
  checking):
  Package 'gurobi' provides an interface to the Gurobi optimization software.
  DESCRIPTION file contains instructions for the installation of Gurobi
  software, and for the installation of the 'gurobi' package. These instructions
  have also been included in the documentation of the function
  `optim_tdv_gurobi_k_2()`, which is the only function from 'diffval' package
  relying on the Gurobi optimization software.
  The use of this software is relevant for the users, as it is the only way to
  guarantee the optimality of a solution for the complex problem being
  addressed.
  
## Test environments

* Local macOS 15.6.1 (24G90) install, R 4.5.0
* macOS builder (r-release-macosx-arm64|4.2.1|macosx|macOS 11.5.2 (20G95)|Mac
  mini|Apple M1||en_US.UTF-8)
* R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R-hub Fedora Linux, R-devel, clang, gfortran
* R-hub Windows Server 2022, R-devel, 64 bit
* Win-builder x86_64-w64-mingw32 (64-bit), R Under development (unstable)
  (2023-03-08 r83956 ucrt)

## R CMD check results

R-hub checks were done setting `_R_CHECK_FORCE_SUGGESTS_ to false`, as package
'gurobi' is not in mainstream repositories.

There were no ERRORs or WARNINGs in all test environments. I list all the NOTEs
that were returned in each of the tested environments:

### Local OS X

There were no NOTEs:

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

### macOS builder 

There was 1 NOTE:

```
 checking CRAN incoming feasibility ... [7s/26s] NOTE
  Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’
  
  Suggests or Enhances not in mainstream repositories:
    gurobi
```
DESCRIPTION has full instructions for 'gurobi' package installation.

### R-hub Ubuntu Linux

There were 3 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’

Suggests or Enhances not in mainstream repositories:
  gurobi
````
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gurobi’
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking examples ... [8s/17s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.234  0.044   7.558
optim_tdv_hill_climb 2.850  0.012   5.997
```
Optimization can be time consuming. This example is small in size.
This is not happening in Windows or OS X environments.

### R-hub Fedora Linux

There were 4 NOTEs:

```
* checking CRAN incoming feasibility ... [8s/29s] NOTE
Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’

Suggests or Enhances not in mainstream repositories:
  gurobi
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘gurobi’
```
DESCRIPTION has full instructions for 'gurobi' package installation.

```
* checking examples ... [8s/17s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.248  0.008   7.543
optim_tdv_hill_climb 2.878  0.023   6.003
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

Suggests or Enhances not in mainstream repositories:
  gurobi
```
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

There was 2 NOTEs:

```
* checking CRAN incoming feasibility ... [13s] NOTE
Maintainer: 'Tiago Monteiro-Henriques <tmh.dev@icloud.com>'

Possibly misspelled words in DESCRIPTION:
  Henriques (18:37)
  Monteiro (18:28)

Suggests or Enhances not in mainstream repositories:
  gurobi
```
Words are spelled correctly.
DESCRIPTION has full instructions for 'gurobi' package installation.

## Downstream dependencies

There are currently no downstream dependencies for this package.
