# ----------------------------------------------------------------------------
# Latitude-dependent physical parameters (Option A: one doubly periodic
# β-plane run per latitude).
#
# Stratification (g′, H1, H2) is held fixed, so the deformation radius scales
# as 1/f₀. The reference latitude is chosen so that the 49°N member reproduces
# the pyqg "eddy" configuration used by Ross et al. (2023) and by
# Perezhogin & Glazunov (2023): β ≈ 1.5e-11 m⁻¹ s⁻¹, Ld = 15 km,
# δ = H1/H2 = 0.25, U1 − U2 = 0.025 m/s, bottom drag 5.787e-7 s⁻¹.
# That anchors their tuned C_R = 7 to one member of the sweep.
# ----------------------------------------------------------------------------

const Ω_EARTH = 7.2921e-5     # s⁻¹
const R_EARTH = 6.371e6       # m

coriolis(φ) = 2Ω_EARTH * sind(φ)
beta_param(φ) = 2Ω_EARTH * cosd(φ) / R_EARTH

Base.@kwdef struct CaseConfig
    lat::Float64     = 49.0      # latitude of the β-plane (°N)
    L::Float64       = 1.0e6     # domain size (m), square and doubly periodic
    H1::Float64      = 500.0     # upper layer thickness (m)
    H2::Float64      = 2000.0    # lower layer thickness (m)
    lat_ref::Float64 = 49.0      # latitude where Ld = Ld_ref
    Ld_ref::Float64  = 15.0e3    # deformation radius at lat_ref (m)
    ΔU_ref::Float64  = 0.025     # U1 − U2 at lat_ref (m/s)
    shear::Symbol    = :fixed_criticality   # or :fixed_shear
    μ::Float64       = 5.787e-7  # linear bottom drag (s⁻¹)
end

"""
    physical_parameters(c::CaseConfig)

Return a NamedTuple with f₀, β, g′, layer depths, deformation radius `Ld`,
stretching coefficients F1, F2, imposed shear ΔU, the layer velocities `U`,
and the supercriticality ξ = ΔU F2 / β (instability requires ξ > 1, the lower
layer PV gradient β − F2 ΔU must change sign).

`shear = :fixed_criticality` keeps ξ equal to its value at `lat_ref`, so every
latitude is equally unstable and differences come from the scales.
`shear = :fixed_shear` keeps ΔU fixed; low latitudes may then be stable.
"""
function physical_parameters(c::CaseConfig)
    H = c.H1 + c.H2
    f_ref, β_ref = coriolis(c.lat_ref), beta_param(c.lat_ref)
    g′ = c.Ld_ref^2 * f_ref^2 * H / (c.H1 * c.H2)
    F2_ref = f_ref^2 / (g′ * c.H2)
    ξ_ref = c.ΔU_ref * F2_ref / β_ref

    f₀, β = coriolis(c.lat), beta_param(c.lat)
    Ld = sqrt(g′ * c.H1 * c.H2 / H) / f₀
    F1 = f₀^2 / (g′ * c.H1)
    F2 = f₀^2 / (g′ * c.H2)
    ΔU = c.shear === :fixed_criticality ? ξ_ref * β / F2 :
         c.shear === :fixed_shear       ? c.ΔU_ref :
         error("shear must be :fixed_criticality or :fixed_shear")
    ξ = ΔU * F2 / β
    c_rossby = β * Ld^2                          # long Rossby wave speed (m/s)
    return (; lat = c.lat, f₀, β, g′, H1 = c.H1, H2 = c.H2, H, Ld, F1, F2,
              ΔU, U = [ΔU, 0.0], ξ, c_rossby, L = c.L, μ = c.μ)
end

"""
    truth_resolution(pp; points_per_Ld = 4, nmin = 256)

Grid points for the eddy-resolving "truth" run: at least `points_per_Ld`
grid points per deformation radius, rounded up to a power of two.
"""
function truth_resolution(pp; points_per_Ld = 4, nmin = 256)
    n = ceil(Int, points_per_Ld * pp.L / pp.Ld)
    return max(nmin, nextpow(2, n))
end

"""
    latitude_table(lats = 20:10:70; coarse_nx = (32, 64, 128), kwargs...)

Print the parameters of each sweep member and Δx/Ld for the coarse grids.
"""
function latitude_table(lats = 20:10:70; coarse_nx = (32, 64, 128), kwargs...)
    @printf("%5s %10s %10s %8s %8s %6s %8s", "lat", "f0", "beta", "Ld[km]", "dU[m/s]", "xi", "nx_true")
    for n in coarse_nx
        @printf(" %9s", "dx/Ld@$n")
    end
    println()
    for φ in lats
        pp = physical_parameters(CaseConfig(; lat = φ, kwargs...))
        @printf("%5.0f %10.3e %10.3e %8.1f %8.3f %6.2f %8d", φ, pp.f₀, pp.β, pp.Ld / 1e3, pp.ΔU, pp.ξ, truth_resolution(pp))
        for n in coarse_nx
            @printf(" %9.2f", (pp.L / n) / pp.Ld)
        end
        println()
    end
end
