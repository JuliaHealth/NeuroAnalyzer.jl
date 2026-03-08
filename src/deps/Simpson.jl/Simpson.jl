module Simpson

export simpson

function basic_simpson(y::AbstractVector, x::Union{AbstractVector, Nothing}=nothing, start::Real=1, stop::Real=length(y)-2, dx::Real=1.0)

    slice0 = start:2:stop
    slice1 = start+1:2:stop+1
    slice2 = start+2:2:stop+2

    if isnothing(x)
        # even-spaced Simpson's rule
        integral = sum(collect(y[slice0] .+ 4 .* y[slice1] .+ y[slice2]))
        integral *= dx / 3.0
    else
        # account for possibly different spacings
        # Simpson's rule changes a bit
        h = diff(x)
        sl0 = start:2:stop
        sl1 = start+1:2:stop+1
        h0 = h[sl0]
        h1 = h[sl1]
        hsum = h0 .+ h1
        hprod = h0 .* h1
        h0divh1 = h0 ./ h1
        tmp = @. hsum / 6.0 .* (y[slice0] * (2 - 1.0 / h0divh1) +
                                y[slice1] * (hsum * hsum / hprod) +
                                y[slice2] * (2 - h0divh1))
        integral = sum(tmp)
    end

    return convert(Float64, integral[1])
end

"""
    simpson(y, x, dx, even)

Integrate y(x) using samples and the composite Simpson's rule.

If x is nothing, spacing of dx is assumed.

If there are an even number of samples, N, then there are an odd number of intervals (N-1), but Simpson's rule requires an even number of intervals. The parameter 'even' controls how this is handled.

The code is based on SciPy v1.7.1: https://github.com/scipy/scipy/blob/v1.7.1/scipy/integrate/_quadrature.py

# Arguments

- `y::AbstractVector` : the vector to be integrated.
- `x::Union{AbstractVector, Nothing}=nothing`, optional : the vector at which `y` is sampled.
- `dx::Real`, optional : spacing of integration points along axis of `x`. Only used when `x` is nothing. Default is 1.0.
- `even::Symbol[:avg, :first, :last]`, optional
    :avg, default : Average two results:
    1) use the first N-2 intervals with a trapezoidal rule on the last interval and 
    2) use the last N-2 intervals with a trapezoidal rule on the first interval.
    :first : Use Simpson's rule for the first N-2 intervals with a trapezoidal rule on the last interval.
    :last : Use Simpson's rule for the last N-2 intervals with a trapezoidal rule on the first interval.

# Notes

For an odd number of samples that are equally spaced the result is exact if the function is a polynomial of order 3 or less. If the samples are not equally spaced, then the result is exact only if the function is a polynomial of order 2 or less.

# Examples
```jldoctest
julia> x = 0:9
julia> y = 0:9
julia> simpson(x, y)
40.5

julia> y = x .^ 3
julia> simpson(y, x)
1642.5

julia> simpson(y, x, even=:first)
1644.5

julia> simpson(y, x, even=:last)
1640.5
```
"""
function simpson(y::AbstractVector, x::Union{AbstractVector, Nothing}=nothing; dx::Real=1.0, even::Symbol=:avg)

    isnothing(x) || (length(x) != length(y) && throw(ArgumentError("If given, length of x must be the same as y.")))
    even in (:avg, :last, :first) || throw(ArgumentError("""Parameter "even" must be :avg, :last, or :first."""))

    N = length(y)
    last_dx = dx
    first_dx = dx

    if N % 2 == 0
        val = 0.0
        integral = 0.0

        # compute using Simpson's rule on first intervals
        if even in (:avg, :first)
            isnothing(x) || (last_dx = x[end] - x[end-1])
            val += 0.5 * last_dx * (y[end] + y[end-1])
            integral = basic_simpson(y, x, 1, N-3, dx)
        end

        # compute using Simpson's rule on last set of intervals
        if even in (:avg, :last)
            isnothing(x) || (first_dx = x[2] - x[1])
            val += 0.5 * first_dx * (y[2] + y[1])
            integral += basic_simpson(y, x, 2, N-2, dx)
        end

        if even === :avg
            val /= 2.0
            integral /= 2.0
        end

        integral = integral + val
    else
        integral = basic_simpson(y, x, 1, N-2, dx)
    end

    return integral
end

end # module