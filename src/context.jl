# ----------------------------------------------------------------------------
# ModelContext: everything a closure needs to know about the model, plus
# template arrays (so closures allocate on the same device as the model).
# ----------------------------------------------------------------------------

struct ModelContext{G, P, T, M}
    grid   :: G
    params :: P
    pp     :: NamedTuple       # physical parameters from physical_parameters()
    f₀     :: T
    β      :: T
    g′     :: T
    H      :: Vector{T}         # layer thicknesses [H1, H2]
    Htot   :: T
    F      :: Vector{T}         # stretching coefficients [F1, F2]
    ΔU     :: T                 # imposed U1 − U2
    U      :: Vector{T}         # imposed layer velocities
    Ld     :: T
    Δx     :: T
    mask   :: M                 # 2/3-rule dealiasing mask on (kr, l)
    real3  :: AbstractArray     # template (nx, ny, 2)
    cplx3  :: AbstractArray     # template (nkr, nl, 2)
    real2  :: AbstractArray     # template (nx, ny)
    cplx2  :: AbstractArray     # template (nkr, nl)
end

function ModelContext(prob, pp)
    grid, params = prob.grid, prob.params
    T = eltype(grid)
    nkr, nl = size(grid.Krsq)
    mask_cpu = ones(T, nkr, nl)
    # no aliased wavenumbers when dealiasing is off (aliased_fraction = 0)
    (grid.kralias === nothing || isempty(grid.kralias)) || (mask_cpu[grid.kralias, :] .= 0)
    (grid.lalias === nothing || isempty(grid.lalias)) || (mask_cpu[:, grid.lalias] .= 0)
    mask = device_array(grid.device)(mask_cpu)
    real3 = similar(prob.vars.q)
    cplx3 = similar(prob.sol)
    real2 = similar(prob.vars.q, grid.nx, grid.ny)
    cplx2 = similar(prob.sol, nkr, nl)
    return ModelContext(grid, params, pp, T(pp.f₀), T(pp.β), T(pp.g′),
                        T[pp.H1, pp.H2], T(pp.H), T[pp.F1, pp.F2], T(pp.ΔU), T.(pp.U),
                        T(pp.Ld), T(grid.dx), mask, real3, cplx3, real2, cplx2)
end

"""
    layerfields(sol, m)

Diagnose resolved fields (anomalies about the imposed zonal flow) from the
PV state `sol = q̂`. Returns spectral (`qh, ψh, ζh`) and physical
(`q, ψ, ζ, u, v`) fields of size (nx, ny, 2).
"""
function layerfields(sol, m::ModelContext)
    qh = copy(sol)
    ψh = similar(qh)
    MultiLayerQG.streamfunctionfrompv!(ψh, qh, m.params, m.grid)
    ζh = @. -m.grid.Krsq * ψh
    uh = @. -im * m.grid.l  * ψh
    vh = @.  im * m.grid.kr * ψh
    return (; qh, ψh, ζh, uh, vh,
              q = tophys(qh, m), ψ = tophys(ψh, m), ζ = tophys(ζh, m),
              u = tophys(uh, m), v = tophys(vh, m))
end
