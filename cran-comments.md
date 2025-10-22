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

* **Local macOS** 15.6.1 (24G90) install, R 4.5.0
* **macOS builder**, r-devel-macosx-arm64|4.5.1|macosx|macOS 13.3.1 (22E261)|Mac
  mini|Apple M1||en_US.UTF-8|macOS 11.3|clang-1403.0.22.14.1|GNU Fortran (GCC)
  14.2.0
* **Win-builder**, x86_64-w64-mingw32, R Under development (unstable)
  (2025-10-21 r88957 ucrt)
* **R-hub** (GitHub) VM **'linux'**, ubuntu-latest
* **R-hub** (GitHub) VM **'m1-san'**, macos-15, ASAN + UBSAN on macOS
* **R-hub** (GitHub) VM **'macos'**, macos-13
* **R-hub** (GitHub) VM **'macos-arm64'**, macos-latest
* **R-hub** (GitHub) VM **'windows'**, windows-latest
* **R-hub** (GitHub) Container **'nosuggests'**, R Under development (unstable) (2025-10-20 r88955) on Fedora Linux 38
* **R-hub** (GitHub) Container **'ubuntu-release'**, R version 4.5.1 (2025-06-13) on Ubuntu 24.04.3 LTS

## R CMD check results

R-hub (GitHub) checks were done setting `extra-packages:` to `gurobi=?ignore`, as
package 'gurobi' is not in mainstream repositories.

There were no ERRORs or WARNINGs in all test environments. I list all the NOTEs
and INFOs that were returned in each of the tested environments:

### Local macOS

There was 1 NOTE:

```
checking CRAN incoming feasibility ... [7s/26s] NOTE
  Maintainer: ‘Tiago Monteiro-Henriques <tmh.dev@icloud.com>’
  
  Suggests or Enhances not in mainstream repositories:
    gurobi

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

DESCRIPTION has full instructions for 'gurobi' package installation.

### macOS builder

There was 1 INFO:

```
* checking package dependencies ... INFO
Package suggested but not available for checking: ‘gurobi’

Status: OK
```

### Win-builder

There was 1 NOTE:

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

### All seven R-hub (GitHub) checks

There was 1 INFO:

```
* checking package dependencies ... INFO
Package suggested but not available for checking: ‘gurobi’

Status: OK
```
## Downstream dependencies

There are currently no downstream dependencies for this package.
