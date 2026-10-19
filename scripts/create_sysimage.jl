using PackageCompiler
using Pkg

@info "Activating packages.."
Pkg.activate(".")
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl#main")
Pkg.instantiate()

@info "Generating sysimage.."

# "generic" builds a sysimage that runs on any CPU (portable, e.g. for
# distribution/CI); switch to "native" if this is only ever built and
# run on the same machine, for the best performance
cpu_target = "generic"

# if a script exercising NeuroAnalyzer's real usage patterns is placed next
# to this file as na_precompile.jl, PackageCompiler will trace it to cache
# actual method compilations (not just package loading) into the sysimage
precompile_file = joinpath(@__DIR__, "na_precompile.jl")

if isfile(precompile_file)
    @info "Using precompile execution file: $precompile_file"
    create_sysimage(["NeuroAnalyzer"],
                     sysimage_path="NeuroAnalyzer.so",
                     cpu_target=cpu_target,
                     precompile_execution_file=precompile_file)
else
    @warn "No na_precompile.jl found next to this script — the sysimage " *
          "will only cache package-load precompilation, not method calls " *
          "from real usage. See the comment above precompile_file."
    create_sysimage(["NeuroAnalyzer"],
                     sysimage_path="NeuroAnalyzer.so",
                     cpu_target=cpu_target)
end

@info "Done."

exit()
