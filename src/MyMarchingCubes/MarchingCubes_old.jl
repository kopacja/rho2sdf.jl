module MarchingCubes

export marchingCubes3D
   
include("Tables.jl")

function marchingCubes3D(
    dists::Vector{Float64},
    Nv::Vector{Int64},
    big::Float64,
)::Vector{Float64}

    idx = findall(x -> x < 0.0, dists)
    while true
        if isempty(idx)
            break
        end
        idx_new = Vector{Int64}()
        sameSign = false
        for I in idx
            I₁ = mod(mod(I - 1, Nv[1] * Nv[2]) + 1 - 1, Nv[1]) + 1
            I₂ = Int(ceil((mod(I - 1, Nv[1] * Nv[2]) + 1) / Nv[1]))
            I₃ = Int(ceil(I / (Nv[1] * Nv[2])))
            for i = I₁-1:I₁+1
                for j = I₂-1:I₂+1
                    for k = I₃-1:I₃+1
                        if minimum([i, j, k]) <= 0 ||
                           minimum(Nv - [i, j, k]) < 0
                            continue
                        end
                        I_adj = Int(
                            (k - 1) .* Nv[1] .* Nv[2] .+ (j - 1) .* Nv[1] .+
                            (i - 1) .+ 1,
                        )

                        if (sign(dists[I_adj]) == sign(dists[I]))
                            sameSign = true
                        end

                        if (abs(dists[I_adj] - big) < 1.0e-5)
                            dists[I_adj] = -big
                            push!(idx_new, I_adj)
                        end
                    end
                end
            end
        end
        if sameSign == false
            #dists[I] = -dists[I]
        end
        idx = idx_new
    end

    idx = findall(x -> x > 0.0, dists)
    while true
        if isempty(idx)
            break
        end
        idx_new = Vector{Int64}()

        sameSign = false
        for I in idx
            I₁ = mod(mod(I - 1, Nv[1] * Nv[2]) + 1 - 1, Nv[1]) + 1
            I₂ = Int(ceil((mod(I - 1, Nv[1] * Nv[2]) + 1) / Nv[1]))
            I₃ = Int(ceil(I / (Nv[1] * Nv[2])))
            for i = I₁-1:I₁+1
                for j = I₂-1:I₂+1
                    for k = I₃-1:I₃+1
                        if minimum([i, j, k]) <= 0 ||
                           minimum(Nv - [i, j, k]) < 0
                            continue
                        end
                        I_adj = Int(
                            (k - 1) .* Nv[1] .* Nv[2] .+ (j - 1) .* Nv[1] .+
                            (i - 1) .+ 1,
                        )

                        if (sign(dists[I_adj]) == sign(dists[I]))
                            sameSign = true
                        end
                        if (abs(dists[I_adj] + big) < 1.0e-5)
                            dists[I_adj] = +big
                            push!(idx_new, I_adj)
                        end
                    end
                end
            end
        end
        if sameSign == false
            #dists[I] = -dists[I]
        end
        idx = idx_new
    end
    return dists
end



end
