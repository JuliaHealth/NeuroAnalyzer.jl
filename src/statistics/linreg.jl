export linreg
export infcrit

"""
    linreg(x, y)

Linear regression between two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Notes

To predict, use: `new_x = DataFrame(x = [3.5, 7]); predict(lr, new_x)

# Returns

Named tuple containing:
- `lr::StatsModels.TableRegressionModel`: model
- `radj::Flpoat64`: R^2
- `c::Vector{Float64}`: coefficients
- `se::Vector{Float64}`: standard error for coefficients
- `aic::Float64`:: Akaike’s Information Criterion (AIC)
- `bic::Float64`:: Bayesian Information Criterion (BIC)
- `lf::Vector{Float64}`: linear fit (plot(x, lf))
"""
function linreg(x::AbstractVector, y::AbstractVector)

    df = DataFrame(:x=>x, :y=>y)
    lr = GLM.lm(@formula(y ~ x), df)
    radj = r2(lr)
    c = coef(lr)
    se = stderror(lr)
    aic, bic = infcrit(lr)
    lf = MultivariateStats.predict(lr)

    return (lr=lr, radj=radj, c=c, se=se, aic=aic, bic=bic, lf=lf)

end

"""
    infcrit(m)

Calculate Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression model `m`.

# Arguments

- `m::StatsModels.TableRegressionModel`: linear regression model

# Returns

Named tuple containing:
- `aic::Float64`
- `bic::Float64`
"""
function infcrit(m::T) where {T<:StatsModels.TableRegressionModel}

    k = length(coef(m)) - 1
    n = length(MultivariateStats.predict(m))
    aic = 2 * k - 2 * log(r2(m))
    bic = k * log(n) - 2 * log(r2(m))

    return (aic=aic, bic=bic)

end
