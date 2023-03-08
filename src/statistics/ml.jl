export mcc

"""
    mcc(tp, tn, fp, fn)

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
function mcc(tp::Int64, tn::Int64, fp::Int64, fn::Int64)
    return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
end
