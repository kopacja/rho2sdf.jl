using Rho2sdf

ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] # 1, 7
Ξ = [0., 1., 1.]
n = Rho2sdf.SignedDistances.RhoNorm(ρₑ, Ξ)


using Plots
# Zkontrolujte, jestli máte nainstalovaný backend GR pro Plots
gr()

function vykresliNormályRhoNorm()
    ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    
    # Definujeme rozsah pro souřadnice Xi
    x_range = -2.0:0.5:2.0
    y_range = -2.0:0.5:2.0
    z_range = -2.0:0.5:2.0
    
    xs, ys, zs, us, vs, ws = [], [], [], [], [], []
    
    for x in x_range
        for y in y_range
            for z in z_range
                Ξ = [x, y, z]
                normal = Rho2sdf.SignedDistances.RhoNorm(ρₑ, Ξ)  # Zde předpokládáme, že vrací [nx, ny, nz]
                
                # Uložení počátečních bodů a komponent normálových vektorů
                push!(xs, x)
                push!(ys, y)
                push!(zs, z)
                push!(us, normal[1])  # Komponenta nx
                push!(vs, normal[2])  # Komponenta ny
                push!(ws, normal[3])  # Komponenta nz
            end
        end
    end
    
    # Vykreslení normálových vektorů jako šipek
    quiver(xs, ys, zs, quiver=(us, vs, ws), legend=false, xlabel="x", ylabel="y", zlabel="z")
end

vykresliNormályRhoNorm()


using GLMakie
GLMakie.activate!()  # Pro GLFW, nebo jiný příkaz podle toho, který backend chcete použít

using GLMakie

function vykresliNormályRhoNormMakie()
    ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

    x_range = -1.0:0.2:1.0
    y_range = -1.0:0.2:1.0
    z_range = -1.0:0.2:1.0

    positions = Point3f[]
    directions = Vec3f[]

    for x in x_range, y in y_range, z in z_range
        Ξ = [x, y, z]
        normal = Rho2sdf.SignedDistances.RhoNorm(ρₑ, Ξ)
        push!(positions, Point3f(x, y, z))
        push!(directions, Vec3f(normal...))
    end

    fig = Figure()
    ax = Axis3(fig[1, 1])
    arrows!(ax, positions, directions, color=:blue, linewidth=2)
    fig
end

vykresliNormályRhoNormMakie()


using GLMakie

function vykresliNormályRhoNormGLMakie()
    ρₑ = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    
    # Definujeme rozsah pro souřadnice Xi
    x_range = -1.0:0.2:1.0
    y_range = -1.0:0.2:1.0
    z_range = -1.0:0.2:1.0
    
    xs, ys, zs, us, vs, ws = [], [], [], [], [], []
    
    for x in x_range
        for y in y_range
            for z in z_range
                Ξ = [x, y, z]
                normal = Rho2sdf.SignedDistances.RhoNorm(ρₑ, Ξ)  # Zde předpokládáme, že vrací [nx, ny, nz]
                
                # Uložení počátečních bodů a komponent normálových vektorů
                push!(xs, x)
                push!(ys, y)
                push!(zs, z)
                push!(us, normal[1])  # Komponenta nx
                push!(vs, normal[2])  # Komponenta ny
                push!(ws, normal[3])  # Komponenta nz
            end
        end
    end
    
    # Příprava pro vykreslení normálových vektorů jako šipek
    fig = Figure()
    ax = Axis3(fig[1, 1], xlabel="x", ylabel="y", zlabel="z")
    quiver!(ax, xs, ys, zs, us, vs, ws, arrow_size=0.00005, linewidth=0.002, color=:blue)
    fig
end
vykresliNormályRhoNormGLMakie()