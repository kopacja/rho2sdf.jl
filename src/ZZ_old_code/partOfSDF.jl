
                        # As = IEN[ISN[sg][a], el] # globální číslo uzlu
                        # adj_els = INE[As] # všechny elementy které jsou součástí tohoto uzlu

                        # #################################x
                        # # hledání vnější stěny?
                        # for adj_el in adj_els # cyklus přes sousední sněny
                        #     common_adj_els = []
                        #     for adj_sg = 1:nes # cyklus přes sousední stěny elementu
                        #         # všechny uzly mají jen jeden společený elementy                               
                        #         common_adj_els = INE[IEN[mesh.ISN[adj_sg][1], adj_el]]
                        #         for b = 2:nsn
                        #             idx = findall(in(INE[IEN[ISN[adj_sg][b], adj_el]]), common_adj_els)
                        #             common_adj_els = common_adj_els[idx]
                        #         end

                        #         if (
                        #             length(common_adj_els) == 1 && # je to jen jeden element?
                        #             in(As, IEN[mesh.ISN[adj_sg], adj_el]) # je As na face ajd ele?
                        #         )
                        #             # println("Adjacent element")

                        #             adj_Xs = X[:, IEN[ISN[adj_sg], adj_el]]
                        #             adj_Xc = mean(adj_Xs, dims = 2) ## ???

                        #             as = indexin(As, IEN[mesh.ISN[adj_sg], adj_el])

                        #             a_prev = ((as[1] + nsn - 1 - 1) % nsn) + 1 # (% zbytek po dělení)
                        #             a_next = ((as[1] + nsn + 1 - 1) % nsn) + 1

                        #             x_prev = X[:, IEN[mesh.ISN[adj_sg][a_prev], adj_el]]
                        #             xs     = X[:, IEN[mesh.ISN[adj_sg][as], adj_el]]
                        #             x_next = X[:, IEN[mesh.ISN[adj_sg][a_next], adj_el]]

                        #             Xt = [x_prev, xs, x_next]
                        #             Xt = reduce(hcat, Xt)

                        #             Et = calculate_triangle_edges(Xt)

                        #             n = cross(Et[1], Et[2])
                        #             n = n / norm(n)

                        #         end
                        #     end
                        # end
                        # 
                        # #################################x
