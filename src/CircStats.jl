"""
module CircStats contains routines for investigating the statistic of circular
data.

By default, all input and output from routines is in radians, but in general passing
`true` as the last argument to a routine will change this to degrees.
"""
module CircStats

using Printf

import Distributions

export
    # Summary stats
    cdist,
    cmean,
    cmedian,
    cresultant,
    cstd,
    cvariance,
    # Distributions
    von_mises_cdf,
    von_mises_pdf,
    # Fitting/estimation
    fit_vonMises,
    # Hypothesis testing
    V_test,
    watson_U2,
    watson_U2n

include("Datasets.jl")
using .Datasets

include("hypothesistests.jl")

"""
    cdist(a, b, degrees::Bool=false) -> angle

Return the angular distance from `a` to `b` (b - a) in the forward
direction; hence if `b` is 'behind' `a`, `distance` will be negative.

Angles are confined to be the smaller possible angle, so are in the range
[-π:π], or [-180°:180°].

Angles are in radians, unless `degrees` == true.
"""
cdist(a, b) = mod(b - a + pi, 2π) - pi
cdist(a, b, degrees::Bool) =
    degrees ? rad2deg(cdist(deg2rad(a), deg2rad(b))) : cdist(a, b)

"""
    cmean(a::Array, degrees::Bool=false) -> mean

Return the mean angle from the set of angles in `a`.

Angles are in radians, unless `degrees` == true.
"""
cmean(a) = atan(sum(sin.(a)), sum(cos.(a)))
cmean(a, degrees::Bool) = degrees ? rad2deg(cmean(deg2rad.(a))) : cmean(a)

"""
    cmedian(a::Array, degrees::Bool=false, axial=false) -> median

Return the median angle from the set of angles `a`.

Angles are in radians, unless `degrees` == true.
"""
function cmedian(a, degrees::Bool=false, axial::Bool=false)
    medians = Vector{Float64}()
    A = Float64.(sort(a))
    n = length(A)
    degrees && (A .= deg2rad.(A))
    axial && (A .= 2A)
    # Compute trial bisectors, which are either the points themselves or adjacent means
    p = Array{Float64}(undef, n)
    if iseven(n)
        for i = 1:n
            j = (i + 1 <= n) ? i + 1 : 1
            p[i] = cmean([A[i], A[j]])
        end
    else
        p[:] .= A[:]
    end
    # Try all possible diameters
    for i = 1:n
        # Count points on either side of diameter
        n_plus = sum(cdist.(p[i], A) .> 0)
        n_minus = sum(cdist.(p[i], A) .< 0)
        if n_plus == n_minus
            # Determine which side of circle is correct direction by counting the number
            # within pi/2 of each of the two opposing possible medians
            if sum(abs.(cdist.(p[i], A)) .<= pi/2) > sum(abs.(cdist.(p[i], A)) .> pi/2)
                push!(medians, A[i])
            else
                push!(medians, (A[i] + 2pi)%2pi - pi)
            end
        end
    end
    # If there is more than one median, take the mean thereof (Otieno & Anderson-Cook
    # (2003), Journal of Modern Applied Statistical Methods, 2(1), 168-176)
    median = if length(medians) > 1
        cmean(medians)
    elseif length(medians) == 1
        medians[1]
    else
        error("Zero medians found.  Are data axial but axial!=true?")
    end
    degrees && (median = rad2deg(median))
    axial && (median /= 2)
    median
end

"""
    cresultant(a, degrees=false) -> R
    cresultant(a, w, degrees=false) -> Rc

Return the resultant vector length, `R`, from a set of angles, `a`.

If data are binned by binwidth `w`, an unbiased estimate of `R`, `Rc` is returned
when `w` is supplied.

Angles are in radians, unless `degrees` == true.
"""
cresultant(a, degrees::Bool=false) = degrees ?
    sqrt(sum(sin.(deg2rad.(a)))^2 + sum(cos.(deg2rad.(a)))^2)/length(a) :
    sqrt(sum(sin.(a))^2 + sum(cos.(a))^2)/length(a)
cresultant(a, w::Real, degrees::Bool=false) = degrees ?
    cresultant(a, true)*deg2rad(w)/(2sin(deg2rad(w)/2)) : cresultant(a)*w/(2sin(w/2))

"""
    cstd(a, degrees=false) -> σ

Return the standard deviation, `σ`, for a set of angles `a`.

Angles are in radians, unless `degrees == true`.
"""
cstd(a, degrees::Bool=false) = sqrt(-2*log(cresultant(a, degrees)))

"""
    cvariance(a, degrees=false) -> σ²

Return the circular variance, `σ²`, of a set of angles `a`.

Angles are in radians, unless `degrees == true`.
"""
function cvariance(a, degrees::Bool=false)
    a_mean = cmean(a, degrees)
    if degrees
        1 - sum(cos.(deg2rad.(a .- a_mean)))/length(a)
    else
        1 - sum(cos.(a .- a_mean))/length(a)
    end
end

"""
    von_mises_pdf(a, µ, κ, degrees=false) -> p

Return the Von Mises probability density function at `a`, for a von Mises
distribution with circular mean `µ` and concentration `κ`.

Angles are in radians, unless `degrees` == true.
"""
von_mises_pdf(a, mu::Real, k::Real, degrees::Bool=false) = degrees ?
                                    exp(k*cos(deg2rad(a - mu)))/(2pi*besseli(0, k)) :
                                    exp(k*cos(a - mu))/(2pi*besseli(0, k))

"""
    von_mises_cdf(a, μ, κ, degrees=false) -> cdf

Return the von Mises cumulative distribution function at
`a` for a distribution with mean `μ` and concentration `κ`.
"""
function von_mises_cdf(a, μ, κ, degrees=false)
    degrees && (μ = deg2rad(μ))
    d = Distributions.VonMises(μ, κ)
    Distributions.cdf(d, degrees ? deg2rad(a) : a)
end

"""
    fit_vonMises(θ, degrees=false; axial=false) -> μ, κ

Return the maximum-likelihood values of the mean `μ` and concentration
`κ` for the von Mises distribution which fits the set of angles `θ`.
"""
function fit_vonMises(θ, degrees=false; axial=false)
    degrees && (θ = deg2rad.(θ))
    axial && (θ = 2θ)
    μ = cmean(θ)
    κ = estimate_kappa_vonMises(θ, μ)
    degrees ? (rad2deg(μ), κ) : (μ, κ)
end

"""
    estimate_kappa_vonMises(θ, μ) -> κ

Return the maximum likelihood estimation of the concentration
parameter `κ` for the best-fitting von Mises distribution, given
a set of angles `θ` in radians.

Uses the polynomial approximation given by Best and Fisher (1981).
**N.B.** This may not be reliable when R̄ is small (e.g., < 0.7).
"""
function estimate_kappa_vonMises(θ, μ)
    R̄ = mean(cos.(θ .- μ))
    if 0 <= R̄ < 0.53
        2R̄ + R̄^3 + 5R̄^5/6
    elseif 0.53 <= R̄ < 0.85
        -0.4 + 1.39R̄ + 0.43/(1 - R̄)
    elseif R̄ >= 0.85
        1/(R̄^3 - 4R̄^2 + 3R̄)
    end
end

end # module
