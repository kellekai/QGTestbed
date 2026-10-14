# ----------------------------------------------------------------------------
# Spectral helpers. Arrays with 3 dims are (x, y, layer); 2 dims are (x, y).
# Transforms always copy their input because c2r FFTs may overwrite it.
# These allocate; fine for a first version, optimise with preallocated
# buffers once the closures are validated.
# ----------------------------------------------------------------------------

function toh(A::AbstractArray{<:Real, 3}, m)
    Ah = similar(m.cplx3)
    MultiLayerQG.fwdtransform!(Ah, copy(A), m.params)
    return Ah
end

function toh(A::AbstractArray{<:Real, 2}, m)
    Ah = similar(m.cplx2)
    mul!(Ah, m.grid.rfftplan, copy(A))
    return Ah
end

function tophys(Ah::AbstractArray{<:Complex, 3}, m)
    A = similar(m.real3)
    MultiLayerQG.invtransform!(A, copy(Ah), m.params)
    return A
end

function tophys(Ah::AbstractArray{<:Complex, 2}, m)
    A = similar(m.real2)
    ldiv!(A, m.grid.rfftplan, copy(Ah))
    return A
end

# derivatives of a physical field, returned in physical space
∂x(A, m) = tophys(@.(im * m.grid.kr * $toh(A, m)), m)
∂y(A, m) = tophys(@.(im * m.grid.l  * $toh(A, m)), m)
∇²(A, m) = tophys(@.(-m.grid.Krsq * $toh(A, m)), m)

"Spectral divergence of a physical flux (Fx, Fy)."
divh(Fx, Fy, m) = @. im * m.grid.kr * $toh(Fx, m) + im * m.grid.l * $toh(Fy, m)

"Gaussian filter of width Δf applied in spectral space: exp(-k²Δf²/24)."
gaussian_transfer(Δf, m) = @. exp(-m.grid.Krsq * Δf^2 / 24)
gfilter(A, Δf, m) = tophys(gaussian_transfer(Δf, m) .* toh(A, m), m)

"""
    curl_of_momentum_laplacian(ν, u, v, m)

Vorticity tendency (spectral) of the momentum tendency ∇·(ν∇u), ∇·(ν∇v)
with spatially varying ν (2D, or broadcastable to the layers), and the local
kinetic-energy tendency density ν(|∇u|² + |∇v|²) with sign convention
"positive = removed from resolved KE" (so negative ν injects energy).
"""
function curl_of_momentum_laplacian(ν, u, v, m)
    ux, uy, vx, vy = ∂x(u, m), ∂y(u, m), ∂x(v, m), ∂y(v, m)
    tu = divh(ν .* ux, ν .* uy, m)                  # spectral ∇·(ν∇u)
    tv = divh(ν .* vx, ν .* vy, m)
    Fζh = @. im * m.grid.kr * tv - im * m.grid.l * tu
    dissipation = @. ν * (ux^2 + uy^2 + vx^2 + vy^2)
    return Fζh, dissipation
end

"""
    curl_of_momentum_biharmonic(ν4, u, v, m)

Vorticity tendency (spectral) of the momentum tendency −∇²(ν4 ∇²u), and the
local dissipation density ν4((∇²u)² + (∇²v)²) ≥ 0. The global energy loss is
exactly the integral of this density, also for spatially varying ν4.
"""
function curl_of_momentum_biharmonic(ν4, u, v, m)
    lu, lv = ∇²(u, m), ∇²(v, m)
    tu = @. m.grid.Krsq * $toh(ν4 .* lu, m)          # spectral of −∇²(ν4∇²u)
    tv = @. m.grid.Krsq * $toh(ν4 .* lv, m)
    Fζh = @. im * m.grid.kr * tv - im * m.grid.l * tu
    dissipation = @. ν4 * (lu^2 + lv^2)
    return Fζh, dissipation
end

"Thickness-weighted vertical average of a (nx, ny, 2) field."
layer_average(A, m) = (m.H[1] .* view(A, :, :, 1) .+ m.H[2] .* view(A, :, :, 2)) ./ m.Htot
