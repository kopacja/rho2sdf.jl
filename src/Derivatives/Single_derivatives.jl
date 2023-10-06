"""
Derivatives for SDF
"""

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
    @einsum dd_dΞ[i] = -dx_dΞ[i, k] * n[k] + X[k] * dn_dΞ[k, i] # ok

    return dd_dΞ
end


function dn_dΞ_compute(
    d²ρ_dΞ²::Matrix,
    norm_dρ_dΞ::Float64,
    dρ_dΞ::Vector)

    Z = zeros(Float64, 3, 3)
    # @einsum Z[i, j] = dρ_dΞ[i] * d²ρ_dΞ²[j, k] * dρ_dΞ[k] # tyto dva zápisy jsou identické
    @einsum Z[i, j] = dρ_dΞ[i] * dρ_dΞ[k] * d²ρ_dΞ²[k, j]   # tyto dva zápisy jsou identické

    dn_dΞ = d²ρ_dΞ² ./ norm_dρ_dΞ - Z ./ norm_dρ_dΞ^3  # ok testovano

    return dn_dΞ
end


function d²x_dΞ²_compute(
    Xₑ::Matrix,
    d²N_dξ²::Array
)

    d²x_dΞ² = zeros(Float64, 3, 3, 3)
    @einsum d²x_dΞ²[i, j, k] = Xₑ[i, m] * d²N_dξ²[m, j, k] # ok testováno

    return d²x_dΞ²
end


function d²n_dΞ²_compute(
    d³ρ_dΞ³::Array,
    norm_dρ_dΞ::Float64,
    d²ρ_dΞ²::Matrix,
    dρ_dΞ::Vector,
)

    Z1 = zeros(Float64, 3, 3, 3)
    @einsum Z1[i, j, k] = d²ρ_dΞ²[i, j] * d²ρ_dΞ²[k, m] * dρ_dΞ[m] # ok 

    # Z2 = zeros(Float64, 3, 3, 3)
    # @einsum Z2[i, j, k] := dρ_dΞ[i] * (d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
    # d³ρ_dΞ³[j, k, m] * dρ_dΞ[m]) # alternativní zápis (nechápu) dle Honzy

    Z2 = zeros(Float64, 3, 3, 3)
    @einsum Z2[i, j, k] = d²ρ_dΞ²[i, j] * dρ_dΞ[l] * d²ρ_dΞ²[l, k] +
                           dρ_dΞ[i] * d²ρ_dΞ²[j, m] * d²ρ_dΞ²[m, k] +
                           dρ_dΞ[i] * dρ_dΞ[m] * d³ρ_dΞ³[j, k, m] # zápis přesně dle studie

    Z3 = zeros(Float64, 3, 3, 3)
    @einsum Z3[i, j, k] = dρ_dΞ[i] * dρ_dΞ[m] * d²ρ_dΞ²[m, j] * dρ_dΞ[l] * d²ρ_dΞ²[l, k] # ok

    d²n_dΞ² = d³ρ_dΞ³ ./ norm_dρ_dΞ -
              Z1 ./ norm_dρ_dΞ^3 -
              Z2 ./ norm_dρ_dΞ^3 +
              3 .* Z3 ./ norm_dρ_dΞ^5 # ok testovano

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
    @einsum d²d_dΞ²[i, j] = -d²x_dΞ²[i, j, k] * n[k] -
        2 * dx_dΞ[i, k] * dn_dΞ[k, j] +
        X[k] * d²n_dΞ²[k, i, j] # ok

    return d²d_dΞ²
end

