include("voronoi_kkt.jl")
K= 2
n = 1
m = 2
seeed = 1
p = 10 #penalty for unsatisfied demand
instance = generate_instance(n, m, seeed)

println("Instance = $instance")

v, x, y, s, ξ = solve_voronoi_kkt(K, instance)
println("v = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ")

ξ[1,:,1]