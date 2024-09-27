function compute_sign_claude(points, X, INE, IEN, ρₙ, ρₜ, ngp, sfce, nel)
  println("sign assignment by claude")
  signs = -1 * ones(ngp)

  @threads for i in 1:ngp
    x = @view points[:, i]

    # Nalezení nejbližšího uzlu
    IDmin = argmin(mapslices(norm, (X .- x), dims=1))[2]
    none = length(INE[IDmin])
    ρₙₑ = ρₙ[IEN[:, INE[IDmin]]]

    #BUG: wtf is this:
    if maximum(ρₙₑ) < ρₜ
      continue
    end
    #BUG: wtf is this:
    if minimum(ρₙₑ) >= ρₜ
      signs[i] = 1.0
      continue
    end

    for j in 1:none
      el = INE[IDmin][j]
      Xₑ = @view X[:, IEN[:, el]]

      min_bounds, max_bounds = compute_aabb(Xₑ)

      #TODO: smazat break a jet cyklus přes všechny elementy a pak vybrat řešení
      # s nejmenší hodnotou lok souřadnic, to však musí být menší jak 1.1
      # možná chat by mi to mohl vysvětlit když mě na to upozorňoval
      #
      #TODO: pokud je bod na hranici dvou elmentů, vypočítat hustotu pro dvoje lok
      # souřadnice -> jsou tam velké skoky??? pokud ano tak vážený průměr asi

      if is_point_inside_aabb(x, min_bounds, max_bounds)
        NotConv, local_coords = find_local_coordinates(sfce, Xₑ, x)

        # Změněno: kontrola lokálních souřadnic
        if NotConv && all(-1.0001 .<= local_coords .<= 1.0001)
          # if NotConv && all(x -> isapprox(x, 1.0, atol=1e-6, rtol=1e-6), abs.(local_coords))
          H, _, _, _ = sfce(local_coords)
          ρₑ = @view ρₙ[IEN[:, el]]
          ρ = dot(H, ρₑ)

          if ρ >= ρₜ
            signs[i] = 1.0
            break  # Přidáno: ukončení cyklu po nalezení vnitřního bodu
          end
        end
      end
    end

  end

  return signs
end
