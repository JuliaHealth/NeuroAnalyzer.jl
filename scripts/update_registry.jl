using Pkg

# usage: julia update_registry.jl <registry folder>
if length(ARGS) != 1
    @error "Usage: julia update_registry.jl <registry folder>"
    exit(1)
end
reg_path = expanduser(ARGS[1])
if !isdir(reg_path)
    @error "The registry folder $reg_path does not exist."
    exit(1)
end

@info "Activating LocalRegistry package"
Pkg.add("LocalRegistry")
using LocalRegistry

@info "Updating registry"
# NeuroAnalyzer.jl is the parent folder of this script
Pkg.activate(dirname(@__DIR__))
register(registry = reg_path, push = false)
@info "Done"
