include("voronoi_kkt.jl")
K = 2
n = 1
m = 2
seeed = 1
p = 100 #penalty for unsatisfied demand
instance = generate_instance(n, m, seeed, demand_bound=1, agg_supply_bound=1)

instance.loc_I = [5 5]
instance.loc_J = [0 5; 10 5]
instance.W = 1
instance.D = 1
instance.pc = 1

println("Instance = $instance")

v, x, y, s, ξ, μ_2, μ_3 = solve_voronoi_kkt(K, instance)
println("v = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ\nμ_2 = $μ_2\nμ_3 = $μ_3")

#ξ[2,:,1]