# Run once per machine (and again after updating the code), from the QGTestbed folder:
#
#   julia setup.jl            # model environment (what the runs need)
#   julia setup.jl plotting   # additionally: separate plotting environment in ./plotting
#
# The package ships without a Project.toml so that unpacking an update never
# overwrites your installed environment. This script creates Project.toml (if
# missing), adds the dependencies, and precompiles QGTestbed, so batch jobs
# don't have to compile it themselves.
#
# The model environment contains no plotting packages, so compute nodes never
# compile them. Plot scripts run with `julia --project=plotting ...`.
using Pkg

const HEADER = """
name = "QGTestbed"
uuid = "56dfe223-75fe-4787-9176-14535bf19f60"
authors = ["Kai Rasmus Keller"]
version = "0.1.0"
"""

proj = joinpath(@__DIR__, "Project.toml")
if !isfile(proj) || !occursin("name = \"QGTestbed\"", read(proj, String))
    # missing, or created by a bare Pkg.add without the package header
    old = isfile(proj) ? read(proj, String) : ""
    write(proj, HEADER * "\n" * old)
    println("Wrote package header to Project.toml")
end

Pkg.activate(@__DIR__)
Pkg.add(["GeophysicalFlows", "FourierFlows", "JLD2", "Statistics", "LinearAlgebra", "Printf", "Random"])
# Remove plotting packages that may have been added to the model environment earlier
for p in ("CairoMakie", "Makie", "Plots", "FFTW")
    haskey(Pkg.project().dependencies, p) && Pkg.rm(p)
end
# For NVIDIA GPUs uncomment:
# Pkg.add("CUDA")
Pkg.instantiate()
Pkg.precompile()
@eval using QGTestbed                       # compile QGTestbed itself now, not in the jobs
println("QGTestbed model environment ready.")

if "plotting" in ARGS
    Pkg.activate(joinpath(@__DIR__, "plotting"))
    # Plots.jl (GR backend) rather than Makie: Makie cannot be precompiled with Julia 1.11.0
    Pkg.add(["Plots", "FFTW", "JLD2", "Statistics", "Printf"])
    Pkg.instantiate()
    Pkg.precompile()
    println("Plotting environment ready: julia --project=plotting scripts/plot_vorticity.jl ...")
end
