using LinearAlgebra, Einsum, Rho2sdf
ρₜ = 0.5

Xₑ = [15.03616 14.74912 11.28356 10.27432 15.03616 14.74912 11.28356 10.27432;
      45.04339 48.95985 48.78704 44.89852 45.04339 48.95985 48.78704 44.89852;
      2.5 2.5 2.5 2.5 5.0 5.0 5.0 5.0]

ρₑ = [0.375, 0.75, 0.5, 0.25, 0.375, 0.5, 0.0, 0.0]

using WGLMakie, ColorSchemes
WGLMakie.activate!()

x = y = z = -1.0:0.1:1.0


function ρ(Ξ::Vector{Float64})
    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = Rho2sdf.sfce(Ξ)
    return H⋅ρₑ
end

vol = [ρ([ix, iy, iz]) for ix in x, iy in y, iz in z]
fig, ax, _ = contour(x, y, z, vol;
    colorrange = (0.0, 1.0),
    levels = 10,
    colormap = :Egypt, transparency = true,
    figure = (; resolution = (1200, 800)),
    axis = (;
        type = Axis3,
        perspectiveness = 0.5,
        azimuth = 2.19,
        elevation = 0.57,
        aspect = (1, 1, 1)
        )
    )
Colorbar(fig)
fig



# x = points[:, v]
x = [7.5, 45.0, 4.0]

xₚ = [0.0, 0.0, 0.0]
n = [0.0, 0.0, 0.0]
Ξ = [0.0, 0.0, 0.0]
λ = 1.0
Ξ_tol = 1e-2
Ξ_norm = 2*Ξ_tol
r_tol = 1e-3
r_norm = 2*r_tol
niter = 100
iter = 1


while (Ξ_norm ≥ Ξ_tol && iter ≤ niter)# || r_norm ≥ r_tol)
    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = Rho2sdf.sfce(Ξ)
    ρ = H⋅ρₑ
    dρ_dΞ = d¹N_dξ¹' * ρₑ
    ΔΞ₁ = dρ_dΞ[1] \ (ρₜ - ρ)
    Ξ[1] = Ξ[1] + ΔΞ₁
    println("Ξ = ", Ξ)
    Ξ_norm = norm(ΔΞ₁)
    println("Ξ_norm = ", Ξ_norm)
    iter += 1
end


Ξ_norm = 2*Ξ_tol
iter = 1
while (Ξ_norm ≥ Ξ_tol && iter ≤ niter)# || r_norm ≥ r_tol)

    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = Rho2sdf.sfce(Ξ)

    xₚ = Xₑ * H
    dx_dΞ = Xₑ * d¹N_dξ¹
    dρ_dΞ = d¹N_dξ¹' * ρₑ

    d²ρ_dΞ² = zeros(Float64,3,3)
    for k in 1:length(H)
        d²ρ_dΞ² +=  ρₑ[k] * d²N_dξ²[k,:,:]
    end

    norm_dρ_dΞ = norm(dρ_dΞ)
    n = dρ_dΞ / norm_dρ_dΞ

    dn_dΞ = zeros(Float64,3,3)
    @einsum dn_dΞ[i,j] := d²ρ_dΞ²[i,j] / norm_dρ_dΞ - (dρ_dΞ[i] * d²ρ_dΞ²[j,k] * dρ_dΞ[k] ) / norm_dρ_dΞ^(3/2)

    dd_dΞ = zeros(Float64,3)
    @einsum dd_dΞ[i] := -dx_dΞ[i,k] * n[k] + (x[k] - xₚ[k]) * dn_dΞ[k,i]


    dL_dΞ = zeros(Float64,3)
    @einsum dL_dΞ[i] := dd_dΞ[i] + λ*dρ_dΞ[i]

    ρ = H⋅ρₑ
    dL_dλ = ρ - ρₜ

    d²x_dΞ² = zeros(Float64,3,3,3)
    @einsum d²x_dΞ²[i,j,k] := Xₑ[i,m] * d²N_dξ²[m,j,k]

    d³ρ_dΞ³ = zeros(Float64,3,3,3)
    for k in 1:length(H)
        d³ρ_dΞ³ +=  ρₑ[k] * d³N_dξ³[k,:,:,:]
    end

    d²n_dΞ² = zeros(Float64,3,3,3)
    @einsum d²n_dΞ²[i,j,k] := d³ρ_dΞ³[i,j,k] / norm_dρ_dΞ - (d²ρ_dΞ²[i,j] * d²ρ_dΞ²[k,m] * dρ_dΞ[m]) / norm_dρ_dΞ^(3/2) - d²ρ_dΞ²[i,j] * d²ρ_dΞ²[k,m] * dρ_dΞ[m] / norm_dρ_dΞ^(3/2) + dρ_dΞ[i] * (d²ρ_dΞ²[j,m] * d²ρ_dΞ²[m,k] + d³ρ_dΞ³[j,k,m] * dρ_dΞ[m]) / norm_dρ_dΞ^(3/2) + 3 * (dρ_dΞ[i] * dρ_dΞ[m] * d²ρ_dΞ²[m,j] * dρ_dΞ[l] * d²ρ_dΞ²[l,k]) / norm_dρ_dΞ^(5/2)

    d²d_dΞ² = zeros(Float64,3,3)
    @einsum d²d_dΞ²[i, j] := -d²x_dΞ²[i,j,k] * n[k] - 2 * dx_dΞ[i,k] * dn_dΞ[k,j] + (x[k] - xₚ[k]) * d²n_dΞ²[k,i,j]
    d²L_dΞ² = d²d_dΞ² + d²ρ_dΞ²*λ

    d²L_dΞdλ = dρ_dΞ
    d²L_dλ² = 0.0

    K = [d²L_dΞ²   d²L_dΞdλ;
         d²L_dΞdλ' d²L_dλ²]
    r = [dL_dΞ; dL_dλ]

    r_norm = norm(r)

    println("K = ",K)
    println("r = ",r)
    ΔΞ_and_Δλ = K \ -r
    println("ΔΞ_and_Δλ = ",ΔΞ_and_Δλ)

    max_abs_Ξ = maximum(abs.(ΔΞ_and_Δλ[1:end-1]))
    if (max_abs_Ξ > 0.5)
        ΔΞ_and_Δλ[1:end-1] = 0.01 * ΔΞ_and_Δλ[1:end-1] / max_abs_Ξ
    end


    Ξ = Ξ + ΔΞ_and_Δλ[1:end-1]
    λ = λ + ΔΞ_and_Δλ[end]

    Ξ_norm = norm(ΔΞ_and_Δλ)

    iter = iter + 1
end
