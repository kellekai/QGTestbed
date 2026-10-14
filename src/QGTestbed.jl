"""
QGTestbed — a latitude-sweep testbed for mesoscale eddy closures in a
2-layer, doubly periodic, β-plane QG model (GeophysicalFlows.jl / MultiLayerQG).

Every closure is an `AbstractClosure` that returns a PV tendency (in Fourier
space, per layer) and optionally carries prognostic state (e.g. eddy energy)
that is stepped once per model time step.

Closures implemented
- `NoClosure`                 : reference, scale-selective spectral filter only
- `NEMODefault`               : NEMO-style biharmonic viscosity + constant GM
- `GEOMETRIC`                 : Mak et al. (2023), eddy-energy constrained GM, optional "splitting"
- `MEKEGMBackscatter`         : Jansen et al. (2019), MEKE with GM + negative-Laplacian backscatter
- `ReynoldsBackscatter`       : Perezhogin & Glazunov (2023), Reynolds-stress (KEB) model on PV
"""
module QGTestbed

using LinearAlgebra
using Printf
using Random
using Statistics
using JLD2
using FourierFlows
using FourierFlows: TwoDGrid
using GeophysicalFlows
using GeophysicalFlows: MultiLayerQG

include("cases.jl")
include("spectral.jl")
include("context.jl")
include("closures/common.jl")
include("closures/nemo_default.jl")
include("closures/geometric.jl")
include("closures/meke_jansen.jl")
include("closures/reynolds_backscatter.jl")
include("problem.jl")
include("diagnostics.jl")
include("coarsegrain.jl")
include("run.jl")

export CaseConfig, physical_parameters, latitude_table
export NoClosure, NEMODefault, GEOMETRIC, MEKEGMBackscatter, ReynoldsBackscatter, make_closure
export build_problem, run_simulation, truth_resolution
export coarsegrain_pv, subgrid_forcing, apriori_skill

end # module
