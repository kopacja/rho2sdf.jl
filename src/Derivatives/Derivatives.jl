module Derivatives

export AnalyticalDerivations, NumericalDerivations

using Einsum
using Statistics
using LinearAlgebra
using Rho2sdf.ShapeFunctions

include("Single_derivatives.jl")


function AnalyticalDerivations(
    Ξ::Vector{Float64},
    Xₑ::Matrix{Float64}, # Souřadnice uzlů jednoho elementu
    ρₑ::Vector{Float64}, # Uzlové hustoty
    λ::Float64,
    ρₜ::Float64, # hraniční hustota
    x::Vector{Float64}, # uzly pravidelné sítě
)

    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ) # tvarové funkce a jejich derivace
    (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, Ξ)

    norm_dρ_dΞ = norm(dρ_dΞ) # ok
    n = dρ_dΞ / norm_dρ_dΞ # ok

    xₚ = Xₑ * H # Xₑ je souřadnice uzlů, H tvarové funkce -> souřadnice bodu KONTROLA!!

    ### 1. Derivace ###
    dx_dΞ = Xₑ * d¹N_dξ¹ # ok

    dn_dΞ = dn_dΞ_compute(d²ρ_dΞ², norm_dρ_dΞ, dρ_dΞ) # ok

    dd_dΞ = dd_dΞ_compute(dx_dΞ, n, x, xₚ, dn_dΞ) # ok

    dL_dΞ = zeros(Float64, 3)
    @einsum dL_dΞ[i] = dd_dΞ[i] + λ * dρ_dΞ[i] # ok = dd_dΞ + λ * dρ_dΞ (alternativní zápis)
    # dL_dΞ = dd_dΞ + λ .* dρ_dΞ 

    ρ = H ⋅ ρₑ # hustota v bodě
    dL_dλ = ρ - ρₜ # ok

    ### 2. Derivace ###
    d²x_dΞ² = d²x_dΞ²_compute(Xₑ, d²N_dξ²) # ok

    d²n_dΞ² = d²n_dΞ²_compute(d³ρ_dΞ³, norm_dρ_dΞ, d²ρ_dΞ², dρ_dΞ) # ok

    d²d_dΞ² = d²d_dΞ²_compute(x, xₚ, d²x_dΞ², n, dx_dΞ, dn_dΞ, d²n_dΞ²) # ok 

    ### Lagrange ###
    d²L_dΞ² = d²d_dΞ² + λ .* d²ρ_dΞ²
    # d²L_dΞ² = λ .* d²ρ_dΞ²
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

    for m in 1:4
        K_diff_col = zeros(Float64, 4)
        for p in 1:2 # = length(ϵ)
            ΔΞ_tmp = zeros(Float64, 4)
            ΔΞ_tmp[m] = ϵ[p]

            H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ_tmp + ΔΞ_tmp)
            (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³) = ρ_derivatives(ρₑ, Ξ_tmp + ΔΞ_tmp)

            norm_dρ_dΞ = norm(dρ_dΞ)
            n = dρ_dΞ / norm_dρ_dΞ

            xₚ = Xₑ * H

            ### 1. Derivace ###
            dx_dΞ = Xₑ * d¹N_dξ¹ # ok

            dn_dΞ = dn_dΞ_compute(d²ρ_dΞ², norm_dρ_dΞ, dρ_dΞ)

            dd_dΞ = dd_dΞ_compute(dx_dΞ, n, x, xₚ, dn_dΞ)

            ### Lagrange ###
            # dL_dΞ = zeros(Float64, 3)
            ΔΞ_Ξ_tmp = Ξ_tmp[end] + ΔΞ_tmp[end]
            # @einsum dL_dΞ[i] := dd_dΞ[i] + ΔΞ_Ξ_tmp *dρ_dΞ[i]
            ##### @einsum dL_dΞ[i] := (Ξ_tmp[end] + ΔΞ_tmp[end]) * dρ_dΞ[i] # ZDE JE CHYBA!!
            dL_dΞ = dd_dΞ + ΔΞ_Ξ_tmp .* dρ_dΞ
            # dL_dΞ = dρ_dΞ .* ΔΞ_Ξ_tmp

            ρ = H ⋅ ρₑ
            dL_dλ = ρ - ρₜ

            r_tmp = [dL_dΞ; dL_dλ]

            K_diff_col = K_diff_col + sign[p] .* r_tmp
            println("K_diff_col:", sign[p] .* r_tmp)
        end
        K_diff[m, :] = K_diff_col ./ (2 * h)
        println("K_diff:", K_diff)
    end

    return K_diff
end

end
