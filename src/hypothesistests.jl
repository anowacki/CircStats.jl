"""
    V_test(θ, θ₀, degrees=false) -> is_significant, V, Vcrit
    
Perform the modified Rayleigh or V test to determine if a sample of angles
`θ` differ signficantly from random towards the angle `θ₀` at the `α` level.

Return `true` if they do differ from random significantly, the test
statistic `V` and the critical value `Vcrit`.
"""
function V_test(θ, θ₀, α=0.05, degrees=false)
    degrees && (θ = deg2rad.(θ); θ₀ = deg2rad(θ₀))
    n = length(θ)
    r = cresultant(θ)
    θ̄ = cmean(θ)
    ϑ = r*cos(θ̄ .- θ₀)
    V = sqrt(2n)*ϑ
    Vcrit = V_critical_value(n, α)
    V > Vcrit, V, Vcrit
end

"""Table of critical values of the V (modified Rayleigh) test.  Table 34 of Kanji (2006)"""
V_TEST_TABLE = [1.3051 1.6524 2.2505 2.4459 2.7938 3.0825
                1.3009 1.6509 2.2640 2.4695 2.8502 3.2114
                1.2980 1.6499 2.2734 2.4858 2.8886 3.2970
                1.2958 1.6492 2.2803 2.4978 2.9164 3.3578
                1.2942 1.6484 2.2856 2.5070 2.9375 3.4034
                1.2929 1.6482 2.2899 2.5143 2.9540 3.4387
                1.2918 1.6479 2.2933 2.5201 2.9672 3.4669
                1.2909 1.6476 2.2961 2.5250 2.9782 3.4899
                1.2902 1.6474 2.2985 2.5290 2.9873 3.5091
                1.2895 1.6472 2.3006 2.5325 2.9950 3.5253
                1.2890 1.6470 2.3023 2.5355 3.0017 3.5392
                1.2885 1.6469 2.3039 2.5381 3.0075 3.5512
                1.2881 1.6467 2.3052 2.5404 3.0126 3.5617
                1.2877 1.6466 2.3064 2.5424 3.0171 3.5710
                1.2874 1.6465 2.3075 2.5442 3.0211 3.5792
                1.2871 1.6464 2.3085 2.5458 3.0247 3.5866
                1.2868 1.6464 2.3093 2.5473 3.0279 3.5932
                1.2866 1.6463 2.3101 2.5486 3.0308 3.5992
                1.2864 1.6462 2.3108 2.5498 3.0335 3.6047
                1.2862 1.6462 2.3115 2.5509 3.0359 3.6096
                1.2860 1.6461 2.3121 2.5519 3.0382 3.6142
                1.2858 1.6461 2.3127 2.5529 3.0402 3.6184
                1.2856 1.6460 2.3132 2.5538 3.0421 3.6223
                1.2855 1.6460 2.3136 2.5546 3.0439 3.6258
                1.2853 1.6459 2.3141 2.5553 3.0455 3.6292
                1.2852 1.6459 2.3145 2.5560 3.0471 3.6323
                1.2843 1.6456 2.3175 2.5610 3.0580 3.6545
                1.2837 1.6455 2.3193 2.5640 3.0646 3.6677
                1.2834 1.6454 2.3205 2.5660 3.0689 3.6764
                1.2831 1.6453 2.3213 2.5674 3.0720 3.6826
                1.2826 1.6452 2.3228 2.5699 3.0775 3.6936
                1.2818 1.6449 2.3256 2.5747 3.0877 3.7140
                1.2817 1.6449 2.3260 2.5752 3.0890 3.7165]
"Values of α (columns) for V_TEST_TABLE"
V_TEST_ALPHAS = (0.1000, 0.0500, 0.0100, 0.0050, 0.0010, 0.0001)
"Values of n (rows) for V_TEST_TABLE"
V_TEST_NS = (5:30..., 40, 50, 60, 70, 100, 500, 1000)

"""
    V_critical_value(n, α) -> Vcrit

Return the critical value for the V test `Vcrit` for a sample size of `n` at the
`α` confidence level.
"""
function V_critical_value(n, α)
    n >= V_TEST_NS[1] || throw(ArgumentError("V test tables only valid for `n >= 5` (have $n)"))
    indn = indmin(abs.(n .- V_TEST_NS))
    indα = indmin(abs.(α .- V_TEST_ALPHAS))
    V_TEST_TABLE[indn,indα]
end

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
    U²crit = watson_U2n_crit(n, α)
    U² > U²crit, U², U²crit
end

function watson_U2n(θ, α::Number=0.05, degrees=false; kwargs...)
    μ, κ = fit_vonMises(θ, degrees; kwargs...)
    watson_U2n(θ, x->von_mises_cdf(x, μ, κ, degrees), α, degrees; kwargs...)
end

"Table of critical values of Watson's U²n test.  Table 35 of Kanji (2006)."
const WATSON_U2N_TABLE = [0.143 0.000 0.161 0.164 0.165
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
"Values of α (columns) for WATSON_U2N_TABLE"
const WATSON_U2N_ALPHAS = (0.100, 0.050, 0.025, 0.010, 0.005)
"Values of n (rows) for WATSON_U2N_table"
const WATSON_U2N_N = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100, 200)

"""
    watson_U2n_crit(n, α) -> U²crit

Return the critical value at the α level for a sample of `n` angles.
"""
function watson_U2n_crit(n, α)
    2 <= n || throw(ArgumentError("`n` (supplied $n) must be 2 or more"))
    any(isapprox.(α, WATSON_U2N_ALPHAS, atol=0.0001)) == 1 ||
        throw(ArgumentError("Significance level for Watson's U² test must " *
                            "be one of $(WATSON_U2_ALPHAS); asked for $α"))
    indn = indmin(abs.(n .- WATSON_U2N_N))
    indα = indmin(abs.(α .- WATSON_U2N_ALPHAS))
    WATSON_U2N_TABLE[indn,indα]
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
