# Helpers shared by the plotting/analysis scripts (plotting environment only).
# Included with `include(joinpath(@__DIR__, "plotting_utils.jl"))`.
using FFTW, Statistics

const YEAR = 365 * 86400.0

"""
    ke_spectrum_from_vorticity(Z, L) -> (k, E)

Isotropic KE spectrum from a relative-vorticity field ζ (n × n, periodic,
side L): E(k) in shells of width dk = 2π/L with Σ E dk = ½⟨u² + v²⟩.
Grid independent: runs at different resolutions can be compared directly.
"""
function ke_spectrum_from_vorticity(Z, L)
    n = size(Z, 1)
    Zh = rfft(Float64.(Z))                       # (n÷2+1, n)
    kr = 2π .* (0:n÷2) ./ L
    l  = 2π .* fftfreq(n, n) ./ L
    K2 = kr .^ 2 .+ (l .^ 2)'
    w = fill(2.0, size(Zh)); w[1, :] .= 1; w[end, :] .= 1   # rfft half-plane weights
    E = ifelse.(K2 .> 0, 0.5 .* abs2.(Zh) ./ K2 ./ n^4 .* w, 0.0)
    dk = 2π / L
    nb = n ÷ 2
    spec = zeros(nb + 1)
    for (idx, e) in pairs(E)
        b = round(Int, sqrt(K2[idx]) / dk)
        b <= nb && (spec[b + 1] += e)
    end
    return collect((0:nb) .* dk), spec ./ dk
end

"""
    mean_spectrum(vort, times, L, spinup_years) -> (k, E) or nothing

Time-mean upper-layer KE spectrum over the vorticity snapshots after the spin-up.
"""
function mean_spectrum(vort, times, L, spinup_years)
    sel = findall(>(spinup_years * YEAR), times)
    isempty(sel) && return nothing
    k, S = nothing, nothing
    for j in sel
        k, s = ke_spectrum_from_vorticity(vort[j], L)
        S = S === nothing ? s : S .+ s
    end
    return k, S ./ length(sel)
end

"Largest wavenumber kept by the 2/3-rule dealiasing on an nx grid of side L."
dealias_cutoff(nx, L) = (2 / 3) * π * nx / L

# --- scoring against the truth (same definitions as scale_awareness.jl) -------------

"""
    load_summary(path; spinup) -> NamedTuple or nothing

Time-mean upper-layer KE spectrum and upper-layer PV flux (resolved +
parametrized) after `spinup` years, plus grid and metadata, from one run file.
"""
function load_summary(path; spinup)
    jldopen(path, "r") do f
        haskey(f, "vorticity_upper") || return nothing
        pp = f["params"]
        res = mean_spectrum(f["vorticity_upper"], f["vorticity_times"], pp.L, spinup)
        res === nothing && return nothing
        ts = f["timeseries"]
        haskey(ts, "vq_resolved_1") || return nothing
        late = ts["t"] .> spinup * YEAR
        vq = mean(ts["vq_resolved_1"][late]) + mean(ts["vq_param_1"][late])
        meta = haskey(f, "metadata") ? f["metadata"] : Dict{String, Any}()
        return (; k = res[1], S = res[2], vq1 = vq, L = pp.L, Ld = pp.Ld, lat = pp.lat,
                  nx = f["nx"], closure = f["closure"], metadata = meta)
    end
end

"""
    truth_reference(dir, lat; spinup = 5.0) -> NamedTuple or nothing

Truth spectrum and PV flux at a latitude, averaged over all truth seeds in `dir`.
"""
function truth_reference(dir, lat; spinup = 5.0)
    tag = string(round(Int, lat), "N_nx")
    files = filter(f -> startswith(f, tag) && occursin(r"_truth_seed\d+\.jld2$", f), readdir(dir))
    sums = filter(!isnothing, [load_summary(joinpath(dir, f); spinup) for f in files])
    isempty(sums) && return nothing
    n = minimum(length(s.S) for s in sums)
    return (k = sums[1].k[1:n], S = mean(s.S[1:n] for s in sums), vq1 = mean(s.vq1 for s in sums),
            nseeds = length(sums))
end

"""
    run_metrics(r, t; max_decades = 6.0)

spectral error : RMS log10(E_run/E_truth) for 0 < k ≤ k_c/2 (capped at `max_decades`)
ke_ratio       : resolved KE (k ≤ k_c) of the run / truth − 1
flux_ratio     : upper-layer PV flux (resolved + parametrized) / truth − 1
"""
function run_metrics(r, t; max_decades = 6.0)
    kc = dealias_cutoff(r.nx, r.L)
    n = min(length(r.S), length(t.S))
    half = [i for i in 2:n if r.k[i] <= kc / 2 && t.S[i] > 0]
    full = [i for i in 2:n if r.k[i] <= kc]
    floor_ = t.S[half] .* 10.0^(-max_decades)
    spec_err = sqrt(mean(log10.(max.(r.S[half], floor_) ./ t.S[half]) .^ 2))
    ke_ratio = sum(r.S[full]) / sum(t.S[full]) - 1
    flux_ratio = r.vq1 / t.vq1 - 1
    return (; spec_err, ke_ratio, flux_ratio)
end

"Weighted cost: w_s·spectral error + w_k·|KE error| + w_f·|flux error|."
cost(m, w) = w[1] * m.spec_err + w[2] * abs(m.ke_ratio) + w[3] * abs(m.flux_ratio)
