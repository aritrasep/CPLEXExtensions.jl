# CPLEXExtensions.jl #

**Build Status:** 
[![Build Status](https://travis-ci.org/aritrasep/CPLEXExtensions.jl.svg?branch=master)](https://travis-ci.org/aritrasep/CPLEXExtensions.jl)

**DOI:** 
[![DOI](https://zenodo.org/badge/100563601.svg)](https://zenodo.org/badge/latestdoi/100563601)

This repository extends [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) for solving both single and multiobjective ([Modof.jl](https://github.com/aritrasep/Modof.jl)) optimization problems.

## Dependencies: ##

1. [Julia v0.6.0](https://julialang.org/downloads/)
1. [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)

## Installation ##

Once, Julia v0.6.0 and [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) has been properly installed, the following instructions in a **Julia** terminal will install **CPLEXExtensions.jl** on the local machine:

```julia
Pkg.clone("https://github.com/aritrasep/CPLEXExtensions.jl")
Pkg.build("CPLEXExtensions")
```

In case `Pkg.build("CPLEXExtensions")` gives you an error on Linux, you may need to install the GMP library headers. For example, on Ubuntu/Debian and similar, give the following command from a terminal:

```
$ sudo apt-get install libgmp-dev
```

After that, restart the installation of the package with:

```
Pkg.build("CPLEXExtensions")
```

## Contributions ##

This package is written and maintained by [Aritra Pal](https://github.com/aritrasep). Please fork and send a pull request or create a [GitHub issue](https://github.com/aritrasep/CPLEXExtensions.jl/issues) for bug reports or feature requests.
