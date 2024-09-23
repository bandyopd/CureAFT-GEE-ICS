
## Code Repository for Mixture Cure Survival Analysis

This repository includes example codes for generating and analyzing
simulation data as described in the manuscript titled “A weighted
generalized estimating equation approach to mixture cure survival with
informative cluster size” by Weixi Zhu, Jonathan W. Yu, Dipankar
Bandyopadhyay, Sy Han Chiou, and Sangwook Kang.

The proposed methods can be conveniently implemented using existing
packages such as `survival`, `aftgee`, `geepack`, `SQUAREM`, and
`daarem`. Except for `survival`, which is a built-in package, the other
packages can be installed using `install.packages()`. Once installed,
these packages can be loaded using the following code snippets, which
include the version numbers.

``` r
> pkgs <- c("aftgee", "survival", "geepack", "SQUAREM", "daarem")
> invisible(sapply(pkgs, require, character.only = TRUE))
> sapply(pkgs, packageVersion, simplify = FALSE)
```

    $aftgee
    [1] '1.2.0'

    $survival
    [1] '3.7.0'

    $geepack
    [1] '1.3.10'

    $SQUAREM
    [1] '2021.1'

    $daarem
    [1] '0.7'

This repository is organized into two folders: `example` and `code`. The
`example` folder contains sample codes for running the simulations, and
can be accessed
[here](https://htmlpreview.github.io/?https://github.com/bandyopd/CureAFT-GEE-ICS/example/run.html)
for easy reference. The code folder includes our implementations of the
proposed method.
