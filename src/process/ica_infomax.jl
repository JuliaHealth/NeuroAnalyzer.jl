# (extended) Infomax ICA (method = :infomax of ica_decompose()): a port of runica() of EEGLAB
#
# 1. Bell AJ, Sejnowski TJ. An information-maximization approach to blind separation and blind
#    deconvolution. Neural Comput. 1995;7(6):1129-1159.
# 2. Lee TW, Girolami M, Sejnowski TJ. Independent component analysis using an extended infomax
#    algorithm for mixed subgaussian and supergaussian sources. Neural Comput.
#    1999;11(2):417-441.
# 3. Delorme A, Makeig S. EEGLAB: an open source toolbox for analysis of single-trial EEG
#    dynamics including independent component analysis. J Neurosci Methods. 2004;134(1):9-21.

# rank of a signal (channels × samples): eigenvalues of the covariance > 1e-7 × the largest
function _data_rank(s::AbstractMatrix)::Int64
    x = s .- mean(s; dims = 2)
    ev = LinearAlgebra.eigvals(LinearAlgebra.Symmetric(x * x' ./ (size(x, 2) - 1)))
    return count(>(1e-7 * maximum(ev)), ev)
end

# runica() defaults
const _INFOMAX_ANNEALDEG = 60.0
const _INFOMAX_ANNEALSTEP = 0.90
const _INFOMAX_EXTANNEAL = 0.98
const _INFOMAX_MAX_WEIGHT = 1e8
const _INFOMAX_RESTART_FAC = 0.9
const _INFOMAX_MIN_LRATE = 1e-6
const _INFOMAX_MAX_KURTSIZE = 6000
const _INFOMAX_SIGNSBIAS = 0.02
const _INFOMAX_EXTMOMENTUM = 0.5
const _INFOMAX_SIGNCOUNT_THRESHOLD = 25
const _INFOMAX_SIGNCOUNT_STEP = 2

"""
    _ica_infomax(s; <keyword arguments>)

Decompose a signal into Independent Components (ICs) using the (extended) Infomax algorithm (as `runica()` of EEGLAB).

The signal is demeaned, reduced to `n` principal components (if `n` < number of channels) and sphered (`2 × cov^(-1/2)`); the unmixing matrix is trained by natural-gradient learning on random blocks of samples (logistic nonlinearity; with `extended = true`, `tanh` with the sign of the kurtosis of each component, so that sub-Gaussian sources, e.g. line noise, are also separated), with annealing of the learning rate and restarts with a lower learning rate if the weights blow up.

# Arguments

- `s::AbstractMatrix`: input signal, shape `(channels, samples)`
- `n::Int64`: number of independent components to extract
- `iter::Int64=512`: maximum number of training steps (passes through the data)
- `extended::Bool=true`: if `true`, use extended Infomax (sub- and super-Gaussian sources)
- `rng::Int64=42`: seed for the random order of the samples

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
"""
function _ica_infomax(
    s::AbstractMatrix;
    n::Int64,
    iter::Int64 = 512,
    extended::Bool = true,
    rng::Int64 = 42,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
}
    # validate
    nch, frames = size(s)
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    n <= nch || throw(ArgumentError("n must be ≤ number of channels."))
    iter >= 1 || throw(ArgumentError("iter must be ≥ 1."))
    frames > n || throw(ArgumentError("Number of samples must be > number of components."))

    rg = Random.MersenneTwister(rng)

    # demean
    x = Matrix{Float64}(s)
    x .-= mean(x; dims = 2)

    # PCA reduction to n components
    pca_proj = if n < nch
        e = LinearAlgebra.eigen(LinearAlgebra.Symmetric(x * x' ./ (frames - 1)))
        E = e.vectors[:, sortperm(e.values; rev = true)[1:n]]
        x = E' * x
        E
    else
        nothing
    end

    # sphering
    c = LinearAlgebra.Symmetric(x * x' ./ (frames - 1))
    ev = LinearAlgebra.eigvals(c)
    minimum(ev) > 1e-7 * maximum(ev) || throw(
        ArgumentError("Data rank ($(_data_rank(s))) < n ($n); use n ≤ data rank."),
    )
    sphere = 2.0 .* inv(sqrt(c))
    sphere = Matrix{Float64}((sphere .+ sphere') ./ 2)
    x = sphere * x

    # training parameters (runica defaults)
    lrate0 = 0.00065 / log(n)
    block = ceil(Int64, min(5 * log(frames), 0.3 * frames))
    annealstep = extended ? _INFOMAX_EXTANNEAL : _INFOMAX_ANNEALSTEP
    nochange = n > 32 ? 1e-7 : 1e-6
    degconst = 180 / π
    kurtsize = min(_INFOMAX_MAX_KURTSIZE, frames)

    _info("Infomax ICA$(extended ? " (extended)" : ""): $n components, block = $block samples, learning rate = $(round(lrate0; sigdigits = 3))")

    w = Matrix{Float64}(LinearAlgebra.I, n, n)
    lrate = lrate0
    converged = false
    step = 0

    progbar = Progress(iter; dt = 1, barlen = 20, color = :white, enabled = progress_bar)
    while true
        # (re)start
        w = Matrix{Float64}(LinearAlgebra.I, n, n)
        bias = zeros(n)
        signs = ones(n)
        extended && (signs[1] = -1.0)
        oldsigns = copy(signs)
        old_kk = zeros(n)
        signcount = 0
        extblocks = 1
        blockno = 0
        oldw = copy(w)
        olddelta = zeros(n * n)
        oldchange = 0.0
        step = 0
        blowup = false

        u = zeros(n, block)
        y = zeros(n, block)
        g = zeros(n, n)
        bi = block * Matrix{Float64}(LinearAlgebra.I, n, n)
        tmp = zeros(n, n)

        while step < iter
            perm = randperm(rg, frames)
            lastt = (frames ÷ block - 1) * block + 1
            for t in 1:block:lastt
                xb = @view x[:, perm[t:(t + block - 1)]]
                # u = w × xb + bias
                LinearAlgebra.mul!(u, w, xb)
                u .+= bias
                if extended
                    y .= tanh.(u)
                    # g = block × I - signs × y × u' - u × u'
                    LinearAlgebra.mul!(g, y, u')
                    g .*= signs
                    LinearAlgebra.mul!(g, u, u', -1.0, -1.0)
                    g .+= bi
                    bias .+= lrate .* vec(sum(-2.0 .* y; dims = 2))
                else
                    y .= 1.0 ./ (1.0 .+ exp.(-u))
                    # g = block × I + (1 - 2y) × u'
                    y .= 1.0 .- 2.0 .* y
                    LinearAlgebra.mul!(g, y, u')
                    g .+= bi
                    bias .+= lrate .* vec(sum(y; dims = 2))
                end
                LinearAlgebra.mul!(tmp, g, w)
                w .+= lrate .* tmp
                blockno += 1

                # kurtosis signs of the components (extended Infomax)
                if extended && blockno % extblocks == 0
                    idx = kurtsize < frames ? rand(rg, 1:frames, kurtsize) : (1:frames)
                    pa = w * @view(x[:, idx])
                    m2 = vec(mean(pa .^ 2; dims = 2)) .^ 2
                    m4 = vec(mean(pa .^ 4; dims = 2))
                    kk = m4 ./ m2 .- 3.0
                    kk = _INFOMAX_EXTMOMENTUM .* old_kk .+ (1.0 - _INFOMAX_EXTMOMENTUM) .* kk
                    old_kk .= kk
                    signs .= sign.(kk .+ _INFOMAX_SIGNSBIAS)
                    if signs == oldsigns
                        signcount += 1
                    else
                        signcount = 0
                    end
                    oldsigns .= signs
                    if signcount >= _INFOMAX_SIGNCOUNT_THRESHOLD
                        extblocks = trunc(Int64, extblocks * _INFOMAX_SIGNCOUNT_STEP)
                        signcount = 0
                    end
                end

                if maximum(abs, w) > _INFOMAX_MAX_WEIGHT || !all(isfinite, w)
                    blowup = true
                    break
                end
            end

            if blowup
                break
            end

            # weight change and annealing
            step += 1
            delta = vec(w .- oldw)
            change = sum(abs2, delta)
            angledelta = 0.0
            if step > 2
                angledelta = acos(clamp(LinearAlgebra.dot(delta, olddelta) / sqrt(change * oldchange), -1.0, 1.0))
            end
            oldw .= w
            if step > 2 && degconst * angledelta > _INFOMAX_ANNEALDEG
                lrate *= annealstep
                olddelta .= delta
                oldchange = change
            elseif step == 1
                olddelta .= delta
                oldchange = change
            end
            progress_bar && next!(progbar)

            if step > 2 && change < nochange
                converged = true
                break
            end
        end

        blowup || break
        # restart with a lower learning rate
        lrate0 *= _INFOMAX_RESTART_FAC
        lrate = lrate0
        lrate < _INFOMAX_MIN_LRATE &&
            throw(ErrorException("Infomax ICA failed: weights blow up even with a learning rate < $(_INFOMAX_MIN_LRATE)."))
        _info("Weights blow up: restarting with learning rate = $(round(lrate; sigdigits = 3))")
    end
    progress_bar && finish!(progbar)

    converged ? _info("Converged after $step steps") :
    _warn("Infomax ICA did not converge after $iter steps; consider a larger iter.")

    # unmixing of the (demeaned) channel signals and mixing matrix
    unmix = w * sphere
    isnothing(pca_proj) || (unmix = unmix * pca_proj')
    xs = Matrix{Float64}(s)
    xs .-= mean(xs; dims = 2)
    ic = unmix * xs
    ic_mw = n == nch ? inv(unmix) : LinearAlgebra.pinv(unmix)

    return (; ic, ic_mw = Matrix{Float64}(ic_mw))
end
