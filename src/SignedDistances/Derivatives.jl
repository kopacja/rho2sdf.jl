"""
Derivatives for SDF
"""


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
    dx_dΞ = d¹N_dξ¹' * Xₑ' # ok

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

    d²L_dΞdλ = dρ_dΞ
    d²L_dλ² = 0.0

    K = [
        d²L_dΞ² d²L_dΞdλ
        d²L_dΞdλ' d²L_dλ²
    ]
    r = [dL_dΞ; dL_dλ] # r1 ... r4 vyčíslím v bode xi_temp

    return K, r, n
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
            dx_dΞ = d¹N_dξ¹' * Xₑ'  # ok

            dn_dΞ = dn_dΞ_compute(d²ρ_dΞ², norm_dρ_dΞ, dρ_dΞ)

            dd_dΞ = dd_dΞ_compute(dx_dΞ, n, x, xₚ, dn_dΞ)

            ### Lagrange ###
            # dL_dΞ = zeros(Float64, 3)
            λ_tmp = Ξ_tmp[end] + ΔΞ_tmp[end]
            # @einsum dL_dΞ[i] := dd_dΞ[i] + ΔΞ_Ξ_tmp *dρ_dΞ[i]
            ##### @einsum dL_dΞ[i] := (Ξ_tmp[end] + ΔΞ_tmp[end]) * dρ_dΞ[i] # ZDE JE CHYBA!!
            dL_dΞ = dd_dΞ + λ_tmp .* dρ_dΞ

            ρ = H ⋅ ρₑ
            dL_dλ = ρ - ρₜ

            r_tmp = [dL_dΞ; dL_dλ]

            K_diff_col = K_diff_col + sign[p] .* r_tmp
        end
        K_diff[m, :] = K_diff_col ./ (2 * h)
    end

    return K_diff
end


#### Single derivatives ####

function ρ_derivatives(ρₑ::Vector, Ξ::Vector)
    
    H, d¹N_dξ¹, d²N_dξ², d³N_dξ³ = sfce(Ξ) # tvarové funkce a jejich derivace

    # dρ_dΞ = zeros(Float64, 3)
    # @einsum dρ_dΞ[i] = ρₑ[k] * d¹N_dξ¹[k, i]
    dρ_dΞ = d¹N_dξ¹' * ρₑ

    # d²ρ_dΞ² = zeros(Float64, 3, 3)
    # @einsum d²ρ_dΞ²[i, j] = ρₑ[k] * d²N_dξ²[k, i, j]

    d²ρ_dΞ² = zeros(Float64, 3, 3)
    for k = 1:length(H)
        d²ρ_dΞ² += ρₑ[k] * d²N_dξ²[k, :, :]
    end
    
    """
    Zvláštní je že oba druhy výpočtu generují
    rozdílné výsledky ale matici K to nezmění.
    """

    # d³ρ_dΞ³ = zeros(Float64, 3, 3, 3)
    # @einsum d³ρ_dΞ³[i, j, k] = ρₑ[m] * d³N_dξ³[m, i, j, k]
    
    d³ρ_dΞ³ = zeros(Float64, 3, 3, 3)
    for k = 1:length(H)
        d³ρ_dΞ³ += ρₑ[k] * d³N_dξ³[k, :, :, :]
    end

    return (dρ_dΞ, d²ρ_dΞ², d³ρ_dΞ³)
end

function dd_dΞ_compute(
    dx_dΞ::Matrix,
    n::Vector,
    x::Vector,
    xₚ::Vector,
    dn_dΞ::Matrix
)

    X = x - xₚ
    dd_dΞ = zeros(Float64, 3)
    @einsum dd_dΞ[i] = -dx_dΞ[i, k] * n[k] + X[k] * dn_dΞ[k, i]

    return dd_dΞ
end


function dn_dΞ_compute(
    d²ρ_dΞ²::Matrix,
    norm_dρ_dΞ::Float64,
    dρ_dΞ::Vector)

    Z = zeros(Float64, 3, 3)
    @einsum Z[i, j] = dρ_dΞ[i] * (dρ_dΞ[k] * d²ρ_dΞ²[k, j])
    
    dn_dΞ = d²ρ_dΞ² ./ norm_dρ_dΞ - Z ./ norm_dρ_dΞ^3 
    
    return dn_dΞ
end


function d²x_dΞ²_compute(
    Xₑ::Matrix,
    d²N_dξ²::Array
)

    d²x_dΞ² = zeros(Float64, 3, 3, 3)
    @einsum d²x_dΞ²[i, j, k] = Xₑ[k, m] * d²N_dξ²[m, j, i] # ok testováno

    return d²x_dΞ²
end


function d²n_dΞ²_compute(
    d³ρ_dΞ³::Array,
    norm_dρ_dΞ::Float64,
    d²ρ_dΞ²::Matrix,
    dρ_dΞ::Vector,
)

    Z1 = zeros(Float64, 3, 3, 3)
    @einsum Z1[i, j, k] = d²ρ_dΞ²[i, k] * (d²ρ_dΞ²[j, m] * dρ_dΞ[m])

    Z2 = zeros(Float64, 3, 3, 3)
    @einsum Z2[i, j, k] = d²ρ_dΞ²[i, j] * (dρ_dΞ[m] * d²ρ_dΞ²[m, k]) +
                           dρ_dΞ[i] * (d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k]) +
                           dρ_dΞ[i] * (dρ_dΞ[m] * d³ρ_dΞ³[m, j, k])

    Z3 = zeros(Float64, 3, 3, 3)
    @einsum Z3[i, j, k] = dρ_dΞ[i] * (dρ_dΞ[m] * d²ρ_dΞ²[m, j]) * (dρ_dΞ[l] * d²ρ_dΞ²[l, k])

    d²n_dΞ² = d³ρ_dΞ³ ./ norm_dρ_dΞ -
              Z1 ./ norm_dρ_dΞ^3 -
              Z2 ./ norm_dρ_dΞ^3 +
              3 .* Z3 ./ norm_dρ_dΞ^5

    return d²n_dΞ²
end

function d²d_dΞ²_compute(
    x::Vector,
    xₚ::Vector,
    d²x_dΞ²::Array,
    n::Vector,
    dx_dΞ::Matrix,
    dn_dΞ::Matrix,
    d²n_dΞ²::Array,
)

    X = x - xₚ
    
    d²d_dΞ² = zeros(Float64, 3, 3)
#    @einsum d²d_dΞ²[i, j] = -d²x_dΞ²[i, j, k] * n[k] -
#        2 * dx_dΞ[j, k] * dn_dΞ[k, i] +
#        X[k] * d²n_dΞ²[k, i, j]

@einsum d²d_dΞ²[i, j] = -d²x_dΞ²[i, j, k] * n[k] - 
         dx_dΞ[j, k] * dn_dΞ[k, i] -
         dx_dΞ[i, k] * dn_dΞ[k, j] +
         X[k] * d²n_dΞ²[k, i, j]

    return d²d_dΞ²
end

