
## Code Repository for Mixture Cure Survival Analysis

This repository includes example codes for generating and analyzing
simulation data as described in the manuscript titled “A weighted
generalized estimating equation approach to mixture cure survival with
informative cluster size” by Weixi Zhu, Jonathan W. Yu, Dipankar
Bandyopadhyay, Sy Han Chiou, and Sangwook Kang.

The proposed methods can be conveniently implemented using existing
packages such as `MASS`, `survival`, `aftgee`, `geepack`, `SQUAREM`, and
`daarem`. The `MASS` and `survival` packages are typically included with
standard R installations. The remaining packages can be installed using
`install.packages()`. Once installed, these packages can be loaded using
the following code snippets, which include the version numbers.

``` r
> pkgs <- c("MASS", "aftgee", "survival", "geepack", "SQUAREM", "daarem")
> invisible(sapply(pkgs, require, character.only = TRUE))
> sapply(pkgs, packageVersion, simplify = FALSE)
```

    $MASS
    [1] '7.3.65'
    
    $aftgee
    [1] '1.2.1'
    
    $survival
    [1] '3.8.6'
    
    $geepack
    [1] '1.3.13'
    
    $SQUAREM
    [1] '2026.1'
    
    $daarem
    [1] '0.7'

This repository is organized into two folders: `example` and `code`. The
`example` folder contains sample codes for running the simulations, and
can be accessed
[here](https://htmlpreview.github.io/?https://github.com/bandyopd/CureAFT-GEE-ICS/main/example/run.html)
for easy reference. The code folder includes our implementations of the
proposed method.
