function _get_lag(n::Int64=10)::Float64
    n = LinRange(0.001, 2.0, n)
    lag = zeros(length(n))
    for idx in eachindex(n)
        t = n[idx]
        lag[idx] = t - @elapsed sleep(t)
    end
    lag = mean(lag)
    return abs(lag)
end

function _mtime(f::String; n::Int64=10)::Float64
    lag = zeros(length(n))
    for idx in eachindex(n)
        t = n[idx]
        lag[idx] = @elapsed eval(Meta.parse(f))
    end
    lag = mean(lag)
    return abs(lag)
end

function _delay(n::Real, l::Real)::Nothing
    n < 0.003 && _warn("For n < 0.003 delay will not be accurate.")
    @assert l >= 0 "l must be ≥ 0."
    n = l < n ? n - l : l - n
    sleep(n)
    return nothing
end

function _check_accuracy()::Nothing
    ac = zeros(100)
    n = LinRange(0.01, 1.0, 100)
    for idx in 1:100
        t = n[idx]
        tt = @elapsed sleep(t)
        t1 = round(t - tt, digits=length(string(t)))
        tt = @elapsed asleep(t, l)
        t2 = round(t - tt, digits=length(string(t)))
        # t2 = t - @elapsed wait(t, l)
        ac[idx] = t2 - t1
    end

    @show mean(ac)
    @show median(ac)

    return nothing
end