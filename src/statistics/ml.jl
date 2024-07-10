export mcc
export f1
export mscr

"""
    mcc(; <keyword arguments>)

Assess performance of the classification model using Matthews correlation coefficient (MCC).

MCC’s value ranges from -1 to 1, depending on:
- a score of -1 denotes a complete discrepancy between expected and actual classes
- 0 is equivalent to making an entirely arbitrary guess
- total agreement between expected and actual classes is indicated by a score of 1

# Arguments

- `tp::Int64`: number of true positives
- `tn::Int64`: number of true negatives
- `fp::Int64`: number of false positives
- `fn::Int64`: number of false negatives

# Returns

- `mcc::Float64`

# Source

https://finnstats.com/index.php/2022/09/06/assess-performance-of-the-classification-model/
"""
function mcc(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)

    return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

end

export f1

"""
    f1(; <keyword arguments>)

Assess performance of the classification model using F1-score. F1-score value ranges from 0 to 1.

# Arguments

- `tp::Int64`: number of true positives
- `tn::Int64`: number of true negatives
- `fp::Int64`: number of false positives
- `fn::Int64`: number of false negatives

# Returns

Named tuple containing:
- `f1::Float64`: F1-score
- `p::Float64`: precision
- `r::Float64`: recall

# Source

https://www.statology.org/what-is-a-good-f1-score/
"""
function f1(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)

    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f1 = 2 * (p * r) / (p + r)

    return (f1=f1, p=p, r=r)

end

"""
    mscr(; <keyword arguments>)

Assess performance of the classification model using misclassification rate.

# Arguments

- `tp::Int64`: number of true positives
- `tn::Int64`: number of true negatives
- `fp::Int64`: number of false positives
- `fn::Int64`: number of false negatives

# Returns

Named tuple containing:
- `mr::Float64`: misclassification rate
- `acc::Float64`: accuracy

# Source

https://www.statology.org/misclassification-rate/
"""
function mscr(; tp::Int64, tn::Int64, fp::Int64, fn::Int64)

    mr = (fp + fn) / (tp + tn + fp + fn)
    acc = 1 - mr

    return (mr=mr, acc=acc)

end
