using Pkg

@info "Activating PackageCompatUI package"
Pkg.add("PackageCompatUI")
using PackageCompatUI

# NeuroAnalyzer.jl is the parent folder of this script
Pkg.activate(dirname(@__DIR__))
compat_ui()
