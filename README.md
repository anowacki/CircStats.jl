# CircStats

[![Build Status](https://travis-ci.org/anowacki/CircStats.jl.svg?branch=master)](https://travis-ci.org/anowacki/CircStats.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ikal62afnbwl4q9d?svg=true)](https://ci.appveyor.com/project/AndyNowacki/circstats-jl)
[![Coverage Status](https://coveralls.io/repos/github/anowacki/CircStats.jl/badge.svg?branch=master)](https://coveralls.io/github/anowacki/CircStats.jl?branch=master)

[Julia](https://julialang.org) functions for the statistical analysis of
circular data.

At present, the repo contains only a minimal amount of functionality that I
routinely use in Julia.  Pull requests for new functionality are welcome,
though I hope that long-term, distributions are contributed to
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl), parametric
estimation goes to [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)'s
`fit` function, and hypothesis testing goes in
[HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl).  I will
contribute the code here when I can and replace the functionality in `CircStats`
with calls to methods elsewhere (as in the case of `von_mises_cdf`).

## Install

```julia
julia> import Pkg; Pkg.add(url="https://github.com/anowacki/CircStats.jl")
```
