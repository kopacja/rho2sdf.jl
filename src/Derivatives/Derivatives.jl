module Derivatives

export AnalyticalDerivations, NumericalDerivations

using Einsum
using Statistics
using LinearAlgebra
using Rho2sdf.ShapeFunctions


function AnalyticalDerivations(
    Ξ::Vector{Float64},
    Xₑ::Matrix{Float64},# Vector{Float64} ??
    ρₑ::Vector{Float64}, # Float64 ??
    λ::Float64,
    ρₜ::Float64,
    x::Vector{Float64}, # ?!? co to je?
)

    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ) # tvarové funkce a jejich derivace

    xₚ = Xₑ * H # Xₑ je souřadnice uzlů, H tvarové funkce -> souřadnice bodu KONTROLA!!
    dx_dΞ = Xₑ * d¹N_dξ¹
    dρ_dΞ = d¹N_dξ¹' * ρₑ

    d²ρ_dΞ² = zeros(Float64, 3, 3)
    for k = 1:length(H)
        d²ρ_dΞ² += ρₑ[k] * d²N_dξ²[k, :, :]
    end

    d³ρ_dΞ³ = zeros(Float64, 3, 3, 3)
    for k = 1:length(H)
        d³ρ_dΞ³ += ρₑ[k] * d³N_dξ³[k, :, :, :]
    end

    norm_dρ_dΞ = norm(dρ_dΞ) # ok
    n = dρ_dΞ / norm_dρ_dΞ # ok

    dn_dΞ = zeros(Float64, 3, 3)
    @einsum dn_dΞ[i, j] := # dve tečky když matice neni alokovaná
        d²ρ_dΞ²[i, j] / norm_dρ_dΞ -
        (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) /
        norm_dρ_dΞ^3  # ok

    dd_dΞ = zeros(Float64, 3)
    @einsum dd_dΞ[i] :=
        -dx_dΞ[i, k] * n[k] +
        (x[k] - xₚ[k]) * dn_dΞ[k, i] # ok


    dL_dΞ = zeros(Float64, 3)
    @einsum dL_dΞ[i] := dd_dΞ[i] + λ * dρ_dΞ[i] # ok

    ρ = H ⋅ ρₑ # hustota v bodě
    dL_dλ = ρ - ρₜ # ok

    d²x_dΞ² = zeros(Float64, 3, 3, 3)
    @einsum d²x_dΞ²[i, j, k] :=
        Xₑ[i, m] * d²N_dξ²[m, j, k] # ok


    d²n_dΞ² = zeros(Float64, 3, 3, 3)
    @einsum d²n_dΞ²[i, j, k] :=
        d³ρ_dΞ³[i, j, k] / norm_dρ_dΞ -
        (d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m]) /
        norm_dρ_dΞ^3 -
        d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m] /
        norm_dρ_dΞ^3 +
        dρ_dΞ[i] * (
            d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
            d³ρ_dΞ³[j, k, m] * dρ_dΞ[m]
        ) / norm_dρ_dΞ^3 +
        3 * (
            dρ_dΞ[i] *
            dρ_dΞ[m] *
            d²ρ_dΞ²[m, j] *
            dρ_dΞ[l] *
            d²ρ_dΞ²[l, k]
        ) / norm_dρ_dΞ^5 # ok

    d²d_dΞ² = zeros(Float64, 3, 3)
    @einsum d²d_dΞ²[i, j] :=
        -d²x_dΞ²[i, j, k] * n[k] -
        2 * dx_dΞ[i, k] * dn_dΞ[k, j] +
        (x[k] - xₚ[k]) * d²n_dΞ²[k, i, j]

    # d²L_dΞ² = d²d_dΞ² + d²ρ_dΞ² * λ
    d²L_dΞ² = d²ρ_dΞ² * λ
    d²L_dΞdλ = dρ_dΞ
    d²L_dλ² = 0.0

    K = [
        d²L_dΞ² d²L_dΞdλ
        d²L_dΞdλ' d²L_dλ²
    ]
    r = [dL_dΞ; dL_dλ] # r1 ... r4 vyčíslím v bode xi_temp

    return K, r
end


function NumericalDerivations(
    Ξ::Vector{Float64},
    Xₑ::Matrix{Float64},# Vector{Float64} ??
    ρₑ::Vector{Float64}, # Float64 ??
    λ::Float64,
    ρₜ::Float64,
    x::Vector{Float64}, # ?!? co to je?
)

    Ξ_tmp = zeros(Float64, 4) # bod pro linearizaci ξ₁,ξ₂,ξ₃,λ
    Ξ_tmp = [Ξ; λ]
    K_diff = zeros(Float64, 4, 4)
    h = 1e-6
    sign = [1 -1]
    ϵ = [h -h]

    for kk in 1:4
        K_diff_col = zeros(Float64, 4)
        for ll in 1:2
            ΔΞ_tmp = zeros(Float64, 4)
            ΔΞ_tmp[kk] = ϵ[ll]

            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ_tmp + ΔΞ_tmp)

            xₚ = Xₑ * H
            dx_dΞ = Xₑ * d¹N_dξ¹

            dρ_dΞ = d¹N_dξ¹' * ρₑ

            d²ρ_dΞ² = zeros(Float64, 3, 3)
            for m in 1:length(H)
                d²ρ_dΞ² += ρₑ[m] * d²N_dξ²[m, :, :]
            end


            norm_dρ_dΞ = norm(dρ_dΞ)
            n = dρ_dΞ / norm_dρ_dΞ

            dn_dΞ = zeros(Float64, 3, 3)
            @einsum dn_dΞ[i, j] := d²ρ_dΞ²[i, j] / norm_dρ_dΞ - (dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k]) / norm_dρ_dΞ^3

            dd_dΞ = zeros(Float64, 3)
            @einsum dd_dΞ[i] := -dx_dΞ[i, k] * n[k] + (x[k] - xₚ[k]) * dn_dΞ[k, i]

            dL_dΞ = zeros(Float64, 3)
            # @einsum dL_dΞ[i] := dd_dΞ[i] + (Ξ_tmp[end]+ΔΞ_tmp[end])*dρ_dΞ[i]
            @einsum dL_dΞ[i] := (Ξ_tmp[end] + ΔΞ_tmp[end]) * dρ_dΞ[i]

            ρ = H ⋅ ρₑ
            dL_dλ = ρ - ρₜ

            r_tmp = [dL_dΞ; dL_dλ]

            K_diff_col = K_diff_col + sign[ll] .* r_tmp
        end
        K_diff[kk, :] = K_diff_col ./ (2 * h)
    end


    return K_diff
end



end
