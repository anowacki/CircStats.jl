# CircStats

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

```bash
git clone https://github.com/anowacki/CircStats.jl
```