module Experiment

using JuMP, COSMO
using Statistics, LinearAlgebra, Random
using Plots
using MATLAB
using FFTW
# Write your package code here.
include("data.jl")
include("position.jl")
include("vector_matrix.jl")
include("solver.jl")
include("basic_function.jl")
include("plot_map.jl")

end
