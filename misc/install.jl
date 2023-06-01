using Pkg
Pkg.add("Revise")
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl")
Pkg.instantiate()
Pkg.resolve()
Pkg.update()
