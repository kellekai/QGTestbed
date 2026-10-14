# ----------------------------------------------------------------------------
# Diagnostics shared by truth and coarse runs
# ----------------------------------------------------------------------------

"Domain-mean KE per layer, PE of the interface, thickness-weighted total."
function energy_summary(prob, m)
    KE, PE = MultiLayerQG.energies(prob)
    return (; KE1 = KE[1], KE2 = KE[2], PE = PE[1], Etot = sum(KE) + sum(PE))
end

"""
    closure_energy_rate(Fqh, f, m)

Rate of change of resolved total energy (vertical mean, m² s⁻³) due to a
PV tendency: dE/dt = −Σ (H_j/H) ⟨ψ_j F_j⟩. Negative = dissipation,
positive = backscatter.
"""
function closure_energy_rate(Fqh, f, m)
    F = tophys(Fqh, m)
    r = 0.0
    for j in 1:2
        r -= m.H[j] / m.Htot * mean(view(f.ψ, :, :, j) .* view(F, :, :, j))
    end
    return r
end

"""
    ke_spectra(f, m)

Isotropic KE spectrum per layer, E(k) in shells of width dk = 2π/L, normalized
so that Σ E dk = ½⟨|∇ψ_j|²⟩ (domain mean). Independent of the grid size, so
truth and coarse runs are directly comparable. (Earlier versions used the
unnormalized FFT amplitudes, which scale with the number of grid points.)
Computed on the CPU.
"""
function ke_spectra(f, m)
    ψh = Array(f.ψh)
    K2 = Array(m.grid.Krsq)
    nx, ny = m.grid.nx, m.grid.ny
    w = fill(2.0, size(K2)); w[1, :] .= 1; iseven(nx) && (w[end, :] .= 1)   # rfft half-plane weights
    dk = 2π / m.grid.Lx
    nb = nx ÷ 2
    k = collect((0:nb) .* dk)
    out = Vector{Vector{Float64}}()
    for j in 1:2
        E = 0.5 .* K2 .* abs2.(view(ψh, :, :, j)) ./ (nx * ny)^2 .* w
        spec = zeros(nb + 1)
        for (idx, e) in pairs(E)
            b = round(Int, sqrt(K2[idx]) / dk)
            b <= nb && (spec[b + 1] += e)
        end
        push!(out, spec ./ dk)
    end
    return k, out
end

"Eddy meridional PV flux ⟨v_j q_j⟩ per layer (the quantity GM-type closures target)."
pv_flux(f) = [mean(view(f.v, :, :, j) .* view(f.q, :, :, j)) for j in 1:2]
