# Run once per machine, from the QGTestbed folder:
#
#   julia setup.jl            # model environment (what the runs need)
#   julia setup.jl plotting   # additionally: separate plotting environment in ./plotting
#
# The model environment contains no plotting packages, so compute nodes never
# compile them. Plot scripts run with `julia --project=plotting ...`.
using Pkg

Pkg.activate(@__DIR__)
Pkg.add(["GeophysicalFlows", "FourierFlows", "JLD2", "Statistics", "LinearAlgebra", "Printf", "Random"])
# For NVIDIA GPUs uncomment:
# Pkg.add("CUDA")
Pkg.instantiate()
Pkg.precompile()
println("QGTestbed model environment ready.")

if "plotting" in ARGS
    Pkg.activate(joinpath(@__DIR__, "plotting"))
    # Plots.jl (GR backend) rather than Makie: Makie cannot be precompiled with Julia 1.11.0
    Pkg.add(["Plots", "FFTW", "JLD2", "Statistics", "Printf"])
    Pkg.instantiate()
    Pkg.precompile()
    println("Plotting environment ready: julia --project=plotting scripts/plot_vorticity.jl ...")
end
