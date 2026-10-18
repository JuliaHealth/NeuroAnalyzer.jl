export linreg
export infcrit

"""
    linreg(x, y)

Fit a simple linear regression model `y ~ x`.

# Arguments

- `x::AbstractVector`: predictor vector; must have the same length as `y` and contain at least 3 elements
- `y::AbstractVector`: response vector

# Returns

Named tuple:

- `lr::StatsModels.TableRegressionModel`: fitted model object
- `c::Vector{Float64}`: coefficients `[intercept, slope]`
- `se::Vector{Float64}`: standard errors of the coefficients
- `R2::Float64`: coefficient of determination R²
- `R2adj::Float64`: adjusted R²
- `aic::Float64`: Akaike Information Criterion (AICc-corrected when `n/k < 40`)
- `bic::Float64`: Bayesian Information Criterion
- `lf::Vector{Float64}`: fitted values (use `plot(x, lf)` to visualise)

# Notes

To predict at new x-values:

```julia
new_x = DataFrame(x = [3.5, 7])
predict(lr, new_x)
```
"""
function linreg(
    x::AbstractVector, y::AbstractVector,
)::@NamedTuple{
    lr::StatsModels.TableRegressionModel,
    c::Vector{Float64},
    se::Vector{Float64},
    R2::Float64,
    R2adj::Float64,
    aic::Float64,
    bic::Float64,
    lf::Vector{Float64},
}
    # validate
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length."))
    length(x) >= 3 || throw(ArgumentError("x and y must contain at least 3 elements."))

    df = DataFrame(:x => x, :y => y)
    lr = GLM.lm(@formula(y ~ x), df)
    c = GLM.coef(lr)
    se = stderror(lr)
    R2, R2adj, aic, bic = infcrit(lr)
    lf = GLM.predict(lr)

    return (; lr, c, se, R2, R2adj, aic, bic, lf)
end

"""
    infcrit(m)

Calculate R², adjusted R², AIC (with small-sample correction), and BIC for a fitted linear regression model.

The small-sample corrected AIC (AICc) is applied when `n / k < 40`: `AICc = AIC + 2k(k+1) / (n − k − 1)`, where `k` is the number of estimated parameters (coefficients, including the intercept, and the residual variance).

# Arguments

- `m::StatsModels.TableRegressionModel`: fitted linear regression model

# Returns

Named tuple:

- `R2::Float64`: coefficient of determination R²
- `R2adj::Float64`: adjusted R²; penalizes for the number of predictors
- `aic::Float64`: AIC (AICc-corrected when `n / k < 40`)
- `bic::Float64`: BIC

# Notes

- `k` = number of estimated parameters (`dof(m)`: coefficients + residual variance).
- AIC = `2k − 2L`; BIC = `k × ln(n) − 2L`, where L = log-likelihood (as in `GLM.aic()`, `GLM.bic()` and R).
- AICc correction is applied when the sample-to-parameter ratio `n/k < 40`.
"""
function infcrit(
    m::T,
)::@NamedTuple{
    R2::Float64,
    R2adj::Float64,
    aic::Float64,
    bic::Float64,
} where {T <: StatsModels.TableRegressionModel}
    # number of estimated parameters (coefficients incl. intercept + residual variance)
    k = GLM.dof(m)
    # number of observations
    n = GLM.nobs(m)
    n > k + 1 || throw(
        ArgumentError(
            "Model has too few observations relative to parameters (n must be > k + 1).",
        ),
    )

    R2 = GLM.r2(m)
    R2adj = GLM.adjr2(m)
    # small-sample AICc when n/k < 40
    aic = n / k < 40 ? GLM.aicc(m) : GLM.aic(m)
    bic = GLM.bic(m)

    return (; R2, R2adj, aic, bic)
end
