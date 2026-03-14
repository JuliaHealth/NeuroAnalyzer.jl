export mcc
export f1
export mscr

"""
    mcc(; <keyword arguments>)

Assess classification model performance using the Matthews Correlation Coefficient (MCC).

MCC is a balanced metric that accounts for all four confusion-matrix cells, making it reliable even for imbalanced classes.

Computed as: `MCC = (tp × tn − fp × fn) / √((tp+fp)(tp+fn)(tn+fp)(tn+fn))`

MCC’s value ranges from -1 to 1, depending on:

- A score of -1 denotes a complete discrepancy between expected and actual classes.
- 0 is equivalent to making an entirely arbitrary guess.
- Total agreement between expected and actual classes is indicated by a score of 1.

# Arguments

- `tp::Int64`: number of true positives; must be ≥ 0
- `tn::Int64`: number of true negatives; must be ≥ 0
- `fp::Int64`: number of false positives; must be ≥ 0
- `fn::Int64`: number of false negatives; must be ≥ 0

# Returns

- `Float64`: MCC ∈ [−1, 1]; returns `0.0` when the denominator is zero (degenerate classifier — all predictions in one class)

# Throws

- `ArgumentError`: if any argument is negative or the total count is zero

# References

https://finnstats.com/index.php/2022/09/06/assess-performance-of-the-classification-model/

# See also

[`f1`](@ref), [`mscr`](@ref)
"""
function mcc(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)::Float64

    !(tp >= 0) && throw(ArgumentError("tp must be ≥ 0."))
    !(tn >= 0) && throw(ArgumentError("tn must be ≥ 0."))
    !(fp >= 0) && throw(ArgumentError("fp must be ≥ 0."))
    !(fn >= 0) && throw(ArgumentError("fn must be ≥ 0."))
    !(tp + tn + fp + fn > 0) && throw(ArgumentError("Total count (tp + tn + fp + fn) must be > 0."))

    denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    denom == 0 && return 0.0
    return (tp * tn - fp * fn) / denom

end

export f1

"""
    f1(; <keyword arguments>)

Assess classification model performance using the F1-score.

The F1-score is the harmonic mean of precision and recall: `F1 = 2 × (precision × recall) / (precision + recall)`

# Arguments
- `tp::Int64`: number of true positives; must be ≥ 0
- `tn::Int64`: number of true negatives; must be ≥ 0
- `fp::Int64`: number of false positives; must be ≥ 0
- `fn::Int64`: number of false negatives; must be ≥ 0

# Returns

Named tuple:

- `f1::Float64`: F1-score ∈ [0, 1]
- `p::Float64`: precision `tp / (tp + fp)`
- `r::Float64`: recall `tp / (tp + fn)`

# Throws
- `ArgumentError`: if any argument is negative, `tp + fp == 0` (undefined precision), or `tp + fn == 0` (undefined recall)

# References

https://www.statology.org/what-is-a-good-f1-score/

# See also

[`mcc`](@ref), [`mscr`](@ref)
"""
function f1(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)::@NamedTuple{f1::Float64, p::Float64, r::Float64}

    !(tp >= 0) && throw(ArgumentError("tp must be ≥ 0."))
    !(tn >= 0) && throw(ArgumentError("tn must be ≥ 0."))
    !(fp >= 0) && throw(ArgumentError("fp must be ≥ 0."))
    !(fn >= 0) && throw(ArgumentError("fn must be ≥ 0."))
    !(tp + fp > 0) && throw(ArgumentError("tp + fp must be > 0 (precision is undefined when no positives are predicted)."))
    !(tp + fn > 0) && throw(ArgumentError("tp + fn must be > 0 (recall is undefined when there are no actual positives)."))

    prec = tp / (tp + fp)
    rec  = tp / (tp + fn)

    # when both precision and recall are 0 the harmonic mean is undefined;
    # return 0.0 by convention rather than NaN
    f1_score = (prec + rec) == 0 ? 0.0 : 2 * (prec * rec) / (prec + rec)

    # note: `tn` is accepted for API consistency but is not used in F1 calculation
    return (f1=f1_score, p=prec, r=rec)

end

"""
    mscr(; <keyword arguments>)

Assess classification model performance using the misclassification rate and accuracy.

- Misclassification rate: `(fp + fn) / (tp + tn + fp + fn)`
- Accuracy: `1 − misclassification rate`

# Arguments
- `tp::Int64`: number of true positives; must be ≥ 0
- `tn::Int64`: number of true negatives; must be ≥ 0
- `fp::Int64`: number of false positives; must be ≥ 0
- `fn::Int64`: number of false negatives; must be ≥ 0

# Returns

Named tuple:

- `mr::Float64`: Misclassification rate ∈ [0, 1].
- `acc::Float64`: Accuracy ∈ [0, 1].

# Throws
- `ArgumentError`: if any argument is negative or total count is zero

# References

https://www.statology.org/misclassification-rate/

# See also

[`mcc`](@ref), [`f1`](@ref)
"""
function mscr(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)::@NamedTuple{mr::Float64, acc::Float64}

    !(tp >= 0) && throw(ArgumentError("tp must be ≥ 0."))
    !(tn >= 0) && throw(ArgumentError("tn must be ≥ 0."))
    !(fp >= 0) && throw(ArgumentError("fp must be ≥ 0."))
    !(fn >= 0) && throw(ArgumentError("fn must be ≥ 0."))
    n = tp + tn + fp + fn
    !(n > 0) && throw(ArgumentError("tp + tn + fp + fn must be > 0."))

    mr  = (fp + fn) / n
    return (mr=mr, acc=1 - mr)

end
