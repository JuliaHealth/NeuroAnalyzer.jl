"""
Gets the position value n in array x
"""
function vsearch(x, n)
    _, n_idx = findmin(abs.(x .- n))
    return n_idx
end

"""
Converts cartographic coordinates to polar
"""
function cart2pol(x, y)
    rho = hypot(x, y)
    theta = atan(y, x)
    return rho, theta
end

"""
Converts polar coordinates to cartographic
"""
function pol2cart(theta, rho)
    x = rho * cos(theta)
    y = rho * sin(theta)
    return x, y
end

function simpson(x, dx)
    return int
end

"""
Calculates absolute band power
"""
function band_power(psd, freqs, f1, f2)
    # find intersecting values in frequency vector
    idx_delta = np.logical_and(freqs >= f1, freqs <= f2)
    # dx: frequency resolution
    p = simpson(psd[idx_delta], dx=freqs[1]-freqs[0])
    return p

"""
Unity based (min..max) data scaling: 0..1
"""
function min_max(x)
    x_min = findmin(x)
    x_max = findmax(x)
    x_s = (x - x_min) / (x - x_max)
    return x_s
end

"""
Returns the derivative of the signal with length same as the signal
"""
function signal_derivative(x)
    x_deriv = np.diff(x)
    x_deriv = np.append(x_deriv, x_deriv[-1])
    return x_deriv
end