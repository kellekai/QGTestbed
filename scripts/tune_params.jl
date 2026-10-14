# Tunable parameters of each closure, by the names used on the command line of
# tune_plan.jl. Each name maps to the closure field it sets and to a conversion
# from the convenient unit typed on the command line to the field's SI value.
# No packages needed; included by tune_plan.jl and tune_member.jl.

const TUNABLE = Dict(
    "nemo" => Dict(
        "Uv"           => (:Uv, identity),               # m/s, bilaplacian ahm = Uv Δx³/12
        "Ue"           => (:Ue, identity),               # m/s, GM aei = Ue Le/2
        "Le_km"        => (:Le, x -> 1e3x),              # km
    ),
    "geometric" => Dict(
        "alpha"        => (:α, identity),                # GEOMETRIC efficiency α
        "lambda_days"  => (:λ, d -> 1 / (86400d)),       # eddy-energy damping time (days)
        "etaE"         => (:ηE, identity),               # m²/s, eddy-energy diffusivity
        "E0"           => (:E0, identity),               # m³/s², minimum eddy energy
        "L_split_km"   => (:L_split, x -> 1e3x),         # km, splitting filter scale
        "Uv"           => (:Uv, identity),               # m/s, bilaplacian viscosity
    ),
    "meke" => Dict(
        "c"            => (:c, identity),                # GM/backscatter/MEKE-diffusion coefficient
        "c_smag"       => (:c_smag, identity),           # biharmonic Smagorinsky coefficient
        "tau_vis_days" => (:τ_vis, d -> 86400d),         # background biharmonic time scale (days)
        "c_d"          => (:c_d, identity),              # bottom-drag coefficient of MEKE
        "U_bg"         => (:U_bg, identity),             # m/s
    ),
    "reynolds" => Dict(
        "C_R"          => (:C_R, identity),              # Reynolds-stress backscatter coefficient
        "FGR"          => (:FGR, identity),              # filter-to-grid ratio
    ),
)
TUNABLE["geometric_nosplit"] = TUNABLE["geometric"]

"Check names and convert `name => value` pairs into closure keyword arguments."
function closure_kwargs_from(closure::AbstractString, pairs)
    haskey(TUNABLE, closure) || error("closure '$closure' has no tunable parameters; choose one of $(sort(collect(keys(TUNABLE))))")
    table = TUNABLE[closure]
    kw = Pair{Symbol, Float64}[]
    for (name, value) in pairs
        haskey(table, name) || error("unknown parameter '$name' for $closure; available: $(join(sort(collect(keys(table))), ", "))")
        field, conv = table[name]
        push!(kw, field => conv(value))
    end
    return (; kw...)
end

"File-name tag for a parameter set, e.g. \"C_R=7\" or \"alpha=0.04_lambda_days=80\"."
tune_tag(pairs) = join(["$(n)=$(v)" for (n, v) in sort(collect(pairs); by = first)], "_")
