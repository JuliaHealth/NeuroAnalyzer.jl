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
- `c::Vector{Float64}`: coefficients
- `se::Vector{Float64}`: standard error for coefficients
- `R2::Float64`: R²
- `R2adj::Float64`: R² adjusted
- `aic::Float64`:: Akaike’s Information Criterion (AIC)
- `bic::Float64`:: Bayesian Information Criterion (BIC)
- `lf::Vector{Float64}`: linear fit (plot(x, lf))
"""
function linreg(x::AbstractVector, y::AbstractVector)::@NamedTuple{lr::StatsModels.TableRegressionModel, c::Vector{Float64}, se::Vector{Float64}, R2::Float64, R2adj::Float64, aic::Float64, bic::Float64, lf::Vector{Float64}}

    @assert length(x) == length(y) "Lengths of x and y must be equal."

    df = DataFrame(:x=>x, :y=>y)
    lr = GLM.lm(@formula(y ~ x), df)
    c = GLM.coef(lr)
    se = stderror(lr)
    R2, R2adj, aic, bic = infcrit(lr)
    lf = GLM.predict(lr)

    return (lr=lr, c=c, se=se, R2=R2, R2adj=R2adj, aic=aic, bic=bic, lf=lf)

end

"""
    infcrit(m)

Calculate R², R² adjusted, Akaike’s Information Criterion (AIC) and Bayesian Information Criterion (BIC) for a linear regression model `m`.

# Arguments

- `m::StatsModels.TableRegressionModel`: linear regression model

# Returns

Named tuple containing:
- `R2::Float64`
- `R2adj::Float64`
- `bic::Float64`
- `bic::Float64`
"""
function infcrit(m::T)::@NamedTuple{R2::Float64, R2adj::Float64, aic::Float64, bic::Float64} where {T<:StatsModels.TableRegressionModel}

    k = length(GLM.coef(m)) - 1
    n = length(GLM.predict(m))
    R2 = GLM.r2(m)
    R2adj = 1 - (1 - R2) * ((n - 1) / (n - k - 1))

    L = GLM.loglikelihood(m)

    aic = 2 * k - 2 * L
    n/k < 40 && (aic = aic + (2 * k * (k + 1)) / (n - k - 1))
    bic = k * log(n) - 2 * L

    return (R2=R2, R2adj=R2adj, aic=aic, bic=bic)

end
