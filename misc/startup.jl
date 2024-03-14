# copy this file to ~/.julia/config

@info "Loading NeuroAnalyzer"
using Revise
using NeuroAnalyzer

ENV["JULIA_REVISE"] = "manual"
# set this according to your system
ENV["JULIA_NUM_THREADS"] = 24
# set this according to your system
ENV["JULIA_NUM_PRECOMPILE_TASKS"] = 24