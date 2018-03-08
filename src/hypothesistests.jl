"""
    watson_U2n(θ, cdf, α=0.05, degrees=false; axial=false) -> test, U², U²crit

Perform Watson's U²n test to determine if a sample of angles `θ`
fit the supplied cdf of a distribution at the `α` level of
significance.  *`test` is `true` if the sample differs significantly*.
Return also the value of the test statistic `U²` and the critical value
`U²crit`.

    watson_U2n(θ, α=0.05, degrees=false; axial=false) -> test, U², U²crit

Perform the test against the maximum-likelihood von Mises distribution
for the sample data.  This is fit using `fit_vonMises()`.
"""
function watson_U2n(θ, cdf::Function, α=0.05, degrees=false; axial=false)
    axial && (θ = 2θ)
    degrees && (θ = deg2rad.(θ))
    θ = sort(mod.(θ, 2π))
    n = length(θ)
    V = cdf.(θ) - cdf(0)
    V̄ = sum(V)/n
    U² = sum(V.^2) - sum((2*(1:n) - 1).*V/n) + n*(1/3 - (V̄ - 1/2)^2)
    U²crit = watson_U2_crit(n, α)
    U² > U²crit, U², U²crit
end

function watson_U2n(θ, α::Number=0.05, degrees=false; kwargs...)
    μ, κ = fit_vonMises(θ, degrees; kwargs...)
    watson_U2n(θ, x->von_mises_cdf(x, μ, κ, degrees), α, degrees; kwargs...)
end

"Table of critical values of Watson's U² test"
const WATSON_U2_TABLE = [0.143 0.000 0.161 0.164 0.165
                         0.145 0.173 0.194 0.213 0.224
                         0.146 0.176 0.202 0.233 0.252
                         0.148 0.177 0.205 0.238 0.262
                         0.149 0.179 0.208 0.243 0.269
                         0.149 0.180 0.210 0.247 0.274
                         0.150 0.181 0.211 0.250 0.278
                         0.150 0.182 0.212 0.252 0.281
                         0.150 0.182 0.213 0.254 0.283
                         0.150 0.183 0.215 0.256 0.287
                         0.151 0.184 0.216 0.258 0.290
                         0.151 0.184 0.216 0.259 0.291
                         0.151 0.184 0.217 0.259 0.292
                         0.151 0.185 0.217 0.261 0.293
                         0.152 0.185 0.219 0.263 0.296
                         0.152 0.186 0.219 0.264 0.298
                         0.152 0.186 0.220 0.265 0.299
                         0.152 0.186 0.221 0.266 0.301
                         0.152 0.187 0.221 0.267 0.302]
"Values of α (columns) for WATSON_U2_TABLE"
const WATSON_U2_ALPHAS = (0.100, 0.050, 0.025, 0.010, 0.005)
"Values of n (rows) for WATSON_U2_table"
const WATSON_U2_N = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100, 200)

"""
    watson_U2_crit(n, α) -> U²crit

Return the critical value at the α level for a sample of `n` angles.
"""
function watson_U2_crit(n, α)
    2 <= n || throw(ArgumentError("`n` (supplied $n) must be 2 or more"))
    any(isapprox.(α, WATSON_U2_ALPHAS, atol=0.0001)) == 1 ||
        throw(ArgumentError("Significance level for Watson's U² test must " *
                            "be one of $(WATSON_U2_ALPHAS); asked for $α"))
    indn = indmin(abs.(n .- WATSON_U2_N))
    indα = indmin(abs.(α .- WATSON_U2_ALPHAS))
    WATSON_U2_TABLE[indn,indα]
end

"""
    test_two_populations_differ(θ₁, θ₂, α, degrees=false; axial=false) -> ::Bool

Perform the Watson-Williams test to determine if two independent
circular measurements differ significantly from each other.

Assumes that samples are drawn from a von Mises distribution with a
concentration parameter κ that is the same, and κ > 2.
"""
function test_two_populations_differ(a, b, α, degrees=false; axial=false)
end
