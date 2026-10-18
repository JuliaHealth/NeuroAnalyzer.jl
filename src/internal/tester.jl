function _get_lag(n::Int64 = 10)::Float64
    n = LinRange(0.001, 2.0, n)
    lag = zeros(length(n))
    for idx in eachindex(n)
        t = n[idx]
        lag[idx] = t - @elapsed sleep(t)
    end
    lag = mean(lag)
    return abs(lag)
end

function _mtime(f::String; n::Int64 = 10)::Float64
    # parse once; time `n` evaluations
    ex = Meta.parse(f)
    lag = zeros(n)
    for idx = 1:n
        lag[idx] = @elapsed eval(ex)
    end
    lag = mean(lag)
    return abs(lag)
end

function _delay(n::Real, l::Real)::Nothing
    n < 0.003 && _warn("For n < 0.003 delay will not be accurate.")
    !(l >= 0) && throw(ArgumentError("l must be ≥ 0."))
    # compensate for the scheduler lag `l`; a lag longer than the delay means no sleep at all
    n = max(n - l, 0)
    n > 0 && sleep(n)
    return nothing
end

function _check_accuracy()::Nothing
    # measured once and used to compensate every lag-corrected delay
    l = _get_lag()
    ac = zeros(100)
    n = LinRange(0.01, 1.0, 100)
    for idx = 1:100
        t = n[idx]
        tt = @elapsed sleep(t)
        t1 = round(t - tt; digits = length(string(t)))
        tt = @elapsed _delay(t, l)
        t2 = round(t - tt; digits = length(string(t)))
        ac[idx] = t2 - t1
    end

    @show mean(ac)
    @show median(ac)

    return nothing
end
