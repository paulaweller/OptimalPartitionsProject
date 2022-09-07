using Dates
include("voronoi_kkt.jl")
include("solve_iter.jl")
K = 2
n = 1
m = 2
seeed = 2
p = 100 #penalty for unsatisfied demand
instance = generate_instance(n, m, seeed)

# instance.loc_I = [1 1]
# instance.loc_J = [5 1; 1 5]
# instance.W = 1
# instance.D = 1
# instance.pc = 1

println("Instance = $instance")

z, v, x, y, s, ξ = solve_voronoi_opt(K, instance)

# z_it, x_it, y_it, p_it, p_true = k_adapt_solution(2, instance)

# open("results.txt","a") do io
#     println(io,"\n$(now())\nInstance = $instance\nz= $z\nv = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ")
#  end

println("z= $z\nv = $v\n x= $x\n y= $y\n s = $s\n ξ = $ξ")

# int_scenarios = enum_uncset(5, instance.loc_J, 0.5)

# I = size(instance.loc_I, 1)
# J = size(instance.loc_J, 1)
# c = reshape([norm(instance.loc_I[i,:]-instance.loc_J[j,:]) for j in 1:J for i in 1:I],I,J)

# objectives_optpart = obs_obj(int_scenarios, find_plan, y, c)
# objectives_it = obs_obj(int_scenarios, find_plan, last(y_it)..., c)

# objectives_for_plot = hcat(objectives_optpart, objectives_it)

# box_plot(objectives_for_plot, ["v variable", "v fixed"])