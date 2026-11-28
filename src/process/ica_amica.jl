# AMICA (method = :amica of ica_decompose()): the AMICA program installed by the user (e.g.
# amica15ub) is run on the signal, as runamica15() of EEGLAB, and its output is read
# (_amica_unmix(), see import_amica())
#
# Palmer JA, Kreutz-Delgado K, Makeig S. AMICA: An Adaptive Mixture of Independent Component
# Analyzers with Shared Components. Technical report, Swartz Center for Computational
# Neuroscience, UCSD; 2011.

"""
    _ica_amica(s; <keyword arguments>)

Decompose a signal into Independent Components (ICs) with the AMICA program.

The signal is written to a temporary binary file (Float32, frame by frame), the parameter file is written (the defaults of `runamica15()` of EEGLAB), the AMICA program is run and its output is read; the temporary files are deleted (the output folder is kept if `outdir` is given).

# Arguments

- `s::AbstractMatrix`: input signal, shape `(channels, samples)`
- `n::Int64`: number of independent components (AMICA `pcakeep`)
- `iter::Int64=2000`: maximum number of iterations (AMICA `max_iter`)
- `amica::String="amica15ub"`: AMICA program (name in the `PATH` or path)
- `models::Int64=1`: number of ICA models (AMICA `num_models`); the components of model 1 are returned
- `threads::Int64=Sys.CPU_THREADS`: number of threads (AMICA `max_threads`)
- `outdir::Union{Nothing, String}=nothing`: folder for the AMICA output (kept, e.g. for `import_amica()` of other models); default: a temporary folder (deleted)

# Returns

Named tuple:

- `ic::Matrix{Float64}`: independent components, shape (n, samples)
- `ic_mw::Matrix{Float64}`: mixing matrix (component weights at the channels), shape (channels, n)
"""
function _ica_amica(
    s::AbstractMatrix;
    n::Int64,
    iter::Int64 = 2000,
    amica::String = "amica15ub",
    models::Int64 = 1,
    threads::Int64 = Sys.CPU_THREADS,
    outdir::Union{Nothing, String} = nothing,
)::@NamedTuple{
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
}
    # validate
    nch, frames = size(s)
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    n <= nch || throw(ArgumentError("n must be ≤ number of channels."))
    iter >= 1 || throw(ArgumentError("iter must be ≥ 1."))
    models >= 1 || throw(ArgumentError("amica_models must be ≥ 1."))
    threads >= 1 || throw(ArgumentError("amica_threads must be ≥ 1."))
    bin = isfile(amica) ? abspath(amica) : Sys.which(amica)
    isnothing(bin) && throw(
        ArgumentError(
            "AMICA program $amica not found: install AMICA (e.g. amica15ub, from the AMICA plugin of EEGLAB or https://github.com/sccn/amica) and give its path with amica = \"…\".",
        ),
    )

    tmp = mktempdir()
    try
        out = isnothing(outdir) ? joinpath(tmp, "amicaout") : abspath(outdir)
        out = endswith(out, "/") ? out : out * "/"
        data_file = joinpath(tmp, "data.fdt")
        param_file = joinpath(tmp, "input.param")
        log_file = joinpath(tmp, "amica.log")

        # signal: Float32, channels × frames (frame by frame)
        write(data_file, htol.(Float32.(vec(Matrix(s)))))

        # parameters (order and defaults of runamica15(); files must precede field_dim)
        params = Pair{String, Any}[
            "files" => data_file,
            "outdir" => out,
            "num_models" => models,
            "num_mix_comps" => 3,
            "pdftype" => 0,
            "block_size" => 128,
            "max_iter" => iter,
            "num_samples" => 1,
            "data_dim" => nch,
            "field_dim" => frames,
            "field_blocksize" => 1,
            "share_comps" => 0,
            "share_start" => 100,
            "comp_thresh" => 0.99,
            "share_iter" => 100,
            "lrate" => 0.1,
            "minlrate" => 1e-8,
            "lratefact" => 0.5,
            "rholrate" => 0.05,
            "rho0" => 1.5,
            "minrho" => 1.0,
            "maxrho" => 2.0,
            "rholratefact" => 0.5,
            "kurt_start" => 3,
            "num_kurt" => 5,
            "kurt_int" => 1,
            "do_newton" => 1,
            "newt_start" => 50,
            "newt_ramp" => 10,
            "newtrate" => 1.0,
            "do_reject" => 0,
            "numrej" => 3,
            "rejsig" => 3.0,
            "rejstart" => 2,
            "rejint" => 3,
            "max_threads" => threads,
            "writestep" => 10,
            "write_nd" => 0,
            "write_LLt" => 1,
            "decwindow" => 1,
            "max_decs" => 3,
            "update_A" => 1,
            "update_c" => 1,
            "update_gm" => 1,
            "update_alpha" => 1,
            "update_mu" => 1,
            "update_beta" => 1,
            "invsigmax" => 100.0,
            "invsigmin" => 0.0,
            "do_rho" => 1,
            "load_rej" => 0,
            "load_W" => 0,
            "load_c" => 0,
            "load_gm" => 0,
            "load_alpha" => 0,
            "load_mu" => 0,
            "load_beta" => 0,
            "load_rho" => 0,
            "load_comp_list" => 0,
            "do_mean" => 1,
            "do_sphere" => 1,
            "doPCA" => 1,
            "pcakeep" => n,
            "pcadb" => 30.0,
            "byte_size" => 4,
            "doscaling" => 1,
            "scalestep" => 1,
            "do_opt_block" => 0,
        ]
        open(param_file, "w") do io
            for (k, v) in params
                println(io, k, " ", v)
            end
        end

        _info("AMICA: $n components of $nch channels, $models model(s), at most $iter iterations ($bin)")

        # the AMICA output (one line per iteration) is shown if verbose
        cmd = `$bin $param_file`
        ok = success(
            pipeline(cmd; stdout = NeuroAnalyzer.verbose ? stdout : log_file, stderr = log_file),
        )
        if !ok || !isfile(joinpath(out, "W"))
            msg = isfile(log_file) ? join(last(readlines(log_file), 10), "\n") : ""
            throw(ErrorException("AMICA failed (see the output above)$(isempty(msg) ? "" : ":\n" * msg)"))
        end

        um = _amica_unmix(out; model = 1)
        x = Matrix{Float64}(s) .- um.mean
        ic = um.unmix * x
        ic_mw = um.nw == nch ? inv(um.unmix) : LinearAlgebra.pinv(um.unmix)
        isnothing(outdir) || _info("AMICA output saved in $out")

        return (; ic, ic_mw = Matrix{Float64}(ic_mw))
    finally
        rm(tmp; recursive = true, force = true)
    end
end
