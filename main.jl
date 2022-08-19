using Dates
include("voronoi_kkt.jl")
K = 3
n = 1
m = 3
seeed = 3
p = 100 #penalty for unsatisfied demand
instance = generate_instance(n, m, seeed)

# instance.loc_I = [5 5]
# instance.loc_J = [0 5; 10 5]
# instance.W = 1
# instance.D = 1
# instance.pc = 1

println("Instance = $instance")

z, v, x, y, s, ξ = solve_voronoi_kkt(K, instance)

open("results.txt","a") do io
    println(io,"\n$(now())\nInstance = $instance\nz= $z\nv = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ")
 end

println("z= $z\nv = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ")

#ξ[2,:,1]