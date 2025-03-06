const tile_ref = [
    SVector(1.0, 0.0, 0.0),   # 0
    SVector(2.0, 2.0, 0.0),   # 1
    SVector(1.0, 4.0, 0.0),   # 2
    SVector(3.0, 4.0, 0.0),   # 3
    SVector(1.0, 0.0, 4.0),   # 4
    SVector(3.0, 0.0, 4.0),   # 5
    SVector(2.0, 1.0, 2.0),   # 6
    SVector(0.0, 2.0, 1.0),   # 7
    SVector(0.0, 2.0, 3.0),   # 8
    SVector(2.0, 2.0, 4.0),   # 9
    SVector(1.0, 4.0, 4.0),   # 10
    SVector(3.0, 4.0, 4.0),   # 11
    SVector(0.0, 4.0, 2.0),   # 12
    SVector(2.0, 3.0, 2.0),   # 13
    SVector(2.0, 5.0, 2.0),   # 14
    SVector(4.0, 2.0, 1.0),   # 15
    SVector(4.0, 2.0, 3.0),   # 16
    SVector(5.0, 4.0, 4.0),   # 17
    SVector(4.0, 4.0, 2.0),   # 18
    SVector(0.0, 2.0, 5.0),   # 19
    SVector(4.0, 2.0, 5.0),   # 20
    SVector(0.0, 0.0, 2.0),   # 21
    SVector(5.0, 0.0, 4.0),   # 22
    SVector(4.0, 0.0, 2.0),   # 23
    SVector(3.0, 0.0, 0.0),   # 24
    SVector(5.0, 0.0, 0.0),   # 25
    SVector(5.0, 4.0, 0.0)    # 26
]

# Definice konektivity tetraedrů (indexy 0-based, převedeme na 1-based při zpracování)
const tetra_connectivity = [
    [3,4,15,14], # 0
    [3,15,13,14],
    [6,21,17,10],
    [6,17,23,24],
    [12,14,17,10],
    [1,25,2,7],
    [21,12,18,17],
    [4,14,16,19],
    [4,15,14,19],
    [14,16,2,4],
    [1,7,8,22],
    [7,16,25,2],
    [9,7,5,22],
    [8,7,9,22],
    [14,12,17,19],
    [12,21,10,17],
    [14,9,13,8],
    [8,3,14,2],
    [14,17,16,19],
    [17,14,7,10],
    [16,7,25,24],
    [17,6,7,24],
    [11,15,14,13],
    [4,3,2,14],
    [9,14,7,8],
    [9,14,11,10],
    [2,8,1,7],
    [14,8,2,7],
    [12,15,19,14],
    [19,27,16,4], # zde se vyskytuje index 27 (0-based), který po posunu dává 28 – předpokládáme, že jde o drobnou nepřesnost v literatuře
    [17,12,18,19],
    [14,9,11,13],
    [14,12,11,10],
    [3,8,14,13],
    [6,17,7,10],
    [5,9,20,10],
    [14,9,7,10],
    [11,9,10,20],
    [7,9,5,10],
    [17,6,23,21],
    [6,7,5,10],
    [15,12,11,14],
    [16,14,2,7],
    [7,17,24,16],
    [26,24,16,25],
    [14,16,17,7]
]

@inline function compute_hex8_shape!(
    # @inline function compute_shape_functions!(
    N::MVector{8,Float64},
    d¹N_dξ¹::MMatrix{8,3,Float64},
    ξ₁::Float64,
    ξ₂::Float64,
    ξ₃::Float64
  )
    # Předpočítání základních výrazů - compiler může optimalizovat do registrů
    ξ₁m1 = ξ₁ - 1
    ξ₁p1 = ξ₁ + 1
    ξ₂m1 = ξ₂ - 1
    ξ₂p1 = ξ₂ + 1
    ξ₃m1 = ξ₃ - 1
    ξ₃p1 = ξ₃ + 1
  
    # Předpočítání společných součinů
    t1 = ξ₁m1 * ξ₂m1  # Pro uzly 1 a 5
    t2 = ξ₁p1 * ξ₂m1  # Pro uzly 2 a 6
    t3 = ξ₁p1 * ξ₂p1  # Pro uzly 3 a 7
    t4 = ξ₁m1 * ξ₂p1  # Pro uzly 4 a 8
  
    # Konstantní koeficient
    coef = 0.125
  
    # Výpočet hodnot shape funkcí - přímo do výstupního vektoru
    N[1] = -coef * t1 * ξ₃m1
    N[2] = coef * t2 * ξ₃m1
    N[3] = -coef * t3 * ξ₃m1
    N[4] = coef * t4 * ξ₃m1
    N[5] = coef * t1 * ξ₃p1
    N[6] = -coef * t2 * ξ₃p1
    N[7] = coef * t3 * ξ₃p1
    N[8] = -coef * t4 * ξ₃p1
  
    # Předpočítání společných výrazů pro derivace
    d1_coef = coef * ξ₃m1
    d1_coefp = coef * ξ₃p1
  
    # Derivace podle ξ₁
    d¹N_dξ¹[1, 1] = -d1_coef * ξ₂m1
    d¹N_dξ¹[2, 1] = d1_coef * ξ₂m1
    d¹N_dξ¹[3, 1] = -d1_coef * ξ₂p1
    d¹N_dξ¹[4, 1] = d1_coef * ξ₂p1
    d¹N_dξ¹[5, 1] = d1_coefp * ξ₂m1
    d¹N_dξ¹[6, 1] = -d1_coefp * ξ₂m1
    d¹N_dξ¹[7, 1] = d1_coefp * ξ₂p1
    d¹N_dξ¹[8, 1] = -d1_coefp * ξ₂p1
  
    # Derivace podle ξ₂
    d¹N_dξ¹[1, 2] = -d1_coef * ξ₁m1
    d¹N_dξ¹[2, 2] = d1_coef * ξ₁p1
    d¹N_dξ¹[3, 2] = -d1_coef * ξ₁p1
    d¹N_dξ¹[4, 2] = d1_coef * ξ₁m1
    d¹N_dξ¹[5, 2] = d1_coefp * ξ₁m1
    d¹N_dξ¹[6, 2] = -d1_coefp * ξ₁p1
    d¹N_dξ¹[7, 2] = d1_coefp * ξ₁p1
    d¹N_dξ¹[8, 2] = -d1_coefp * ξ₁m1
  
    # Derivace podle ξ₃
    d¹N_dξ¹[1, 3] = -coef * t1
    d¹N_dξ¹[2, 3] = coef * t2
    d¹N_dξ¹[3, 3] = -coef * t3
    d¹N_dξ¹[4, 3] = coef * t4
    d¹N_dξ¹[5, 3] = coef * t1
    d¹N_dξ¹[6, 3] = -coef * t2
    d¹N_dξ¹[7, 3] = coef * t3
    d¹N_dξ¹[8, 3] = -coef * t4
  end
  
  

function hex8_shape(Ξ::Vector{Float64})
    ξ₁ = Ξ[1]
    ξ₂ = Ξ[2]
    ξ₃ = Ξ[3]
    # Předpočítání základních výrazů - compiler může optimalizovat do registrů
    ξ₁m1 = ξ₁ - 1
    ξ₁p1 = ξ₁ + 1
    ξ₂m1 = ξ₂ - 1
    ξ₂p1 = ξ₂ + 1
    ξ₃m1 = ξ₃ - 1
    ξ₃p1 = ξ₃ + 1
    # Předpočítání společných součinů
    t1 = ξ₁m1 * ξ₂m1  # Pro uzly 1 a 5
    t2 = ξ₁p1 * ξ₂m1  # Pro uzly 2 a 6
    t3 = ξ₁p1 * ξ₂p1  # Pro uzly 3 a 7
    t4 = ξ₁m1 * ξ₂p1  # Pro uzly 4 a 8
    # Konstantní koeficient
    coef = 0.125
    N = zeros(Float64, 8)
    # Výpočet hodnot shape funkcí - přímo do výstupního vektoru
    N[1] = -coef * t1 * ξ₃m1
    N[2] = coef * t2 * ξ₃m1
    N[3] = -coef * t3 * ξ₃m1
    N[4] = coef * t4 * ξ₃m1
    N[5] = coef * t1 * ξ₃p1
    N[6] = -coef * t2 * ξ₃p1
    N[7] = coef * t3 * ξ₃p1
    N[8] = -coef * t4 * ξ₃p1
    return N
  end
  
