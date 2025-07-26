using DynamicalSystems, Statistics, SharedArrays, OrdinaryDiffEq, DelimitedFiles, DataStructures
@inline @inbounds function ideia(u, par, t)
    eps1 = 1.0
    eps2 = 1.0
    alpha = 2.0
    beta = -1.5
    delta = 1.0
    w = par[1]
    a0 = 0.0
    b0 = 0.0
    
    k = 0.5*cos(w*t+a0*cos(b0*w*t))

    du1 = u[2]
    du2 = (eps1 - (u[1] + beta*u[3])^2) * u[2] - (u[1] + beta*u[3])
    du3 = u[4]
    du4 = (eps2 - (u[3] + alpha*u[1])^2) * u[4] - (1 + delta)*(u[3] + alpha*u[1]) + k

    return SVector{4}(du1, du2, du3, du4)
end


function scatter_attractors(attractors, basins, index)
   # for k in keys(attractors)
   #     x, y = columns(attractors[k])
   #     writedlm("attractor_$(k)_01.txt", [x y])
   # end
    writedlm("$(index)_basins_.txt", basins)
end
#alpha = SharedArray{Float64}(length(theta))
#N_attractor = SharedArray{Float64}(length(theta))
#Sb = SharedArray{Float64}(length(theta))
#Sbb = SharedArray{Float64}(length(theta))

omega = range(0.0,1.0,length=100)

for i =97:1:length(omega)
    x0     = range(-4.5, 4.5, length = 800)
    y0     = range(-4.5, 4.5, length = 800)

    u0 = [0.0,0.0,0.0,0.0]
    diffeq    = (alg = Vern9(), abstol = 1e-3, reltol = 1e-3)
    ds        = ContinuousDynamicalSystem(ideia,u0,[omega[i]];diffeq)
    psys = ProjectedDynamicalSystem(ds, [1, 3], [0.0, 0.0])    
    #psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0])    
    mapper = AttractorsViaRecurrences(psys, (x0, y0); sparse = false, Δt=0.01, consecutive_recurrences = 1000, attractor_locate_steps = 1000, consecutive_lost_steps = 100)
    basins, attractors = basins_of_attraction(mapper; show_progress = true)

    scatter_attractors(attractors, basins, i)
end
#eps,f_eps, alpha[i] = uncertainty_exponent(basins)
#N_attractor[i] = maximum(basins)
#Sb[i],Sbb[i] = basin_entropy(basins)
#@show i, alpha[i], N_attractor[i], Sb[i],Sbb[i]

#writedlm("ALPHA.dat",[alpha Sb Sbb])
