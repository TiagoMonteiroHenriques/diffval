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
* checking examples ... [7s/28s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.383  0.028  12.895
optim_tdv_hill_climb 2.583  0.019   9.804
```
Optimization can be time consuming. This example is small in size.
This is not happening in Windows or OS X environments.

### R-hub Fedora Linux

There were 4 NOTEs:
```
* checking CRAN incoming feasibility ... [6s/31s] NOTE
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
* checking examples ... [8s/29s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
optim_tdv_simul_anne 3.410  0.043  14.184
optim_tdv_hill_climb 2.584  0.005   9.527
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

There was 1 NOTE:

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
```
This is a first submission.
Spell check was performed. These terms are correct.
DESCRIPTION has full instructions for 'gurobi' package installation.

## Downstream dependencies

There are currently no downstream dependencies for this package.
