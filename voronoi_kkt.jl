using JuMP, Gurobi, LinearAlgebra
include("helpers.jl")


"""
    solve_voronoi_kkt(K, inst)

    Solve the K-adaptable pre-allocation problem for instance "inst" with Voronoi-style partitioning.
"""

function solve_voronoi_kkt(K, inst)

    I_size = size(inst.loc_I, 1)
    J_size = size(inst.loc_J, 1)
    c = reshape([norm(inst.loc_I[i,:]-inst.loc_J[j,:]) for j in 1:J_size for i in 1:I_size],I_size,J_size)
    

    part_model = Model(Gurobi.Optimizer)

    # variables ---------------------------------

    # supply
    @variable(part_model, 0<=x[1:I_size])
    # allocation of supplies
    @variable(part_model, 0<=y[1:I_size,1:J_size,1:K])
    # unsatisfied demand
    @variable(part_model, 0<=s[1:J_size,1:K])
    # objectives of cells
    @variable(part_model, z[1:K])
    # overall objective
    @variable(part_model, z_obj)
    # voronoi points
    @variable(part_model, v[1:K])

    # lower-level variables
    # worst-case demand scenario, one per partition k and row j
    @variable(part_model, 0<= ξ[1:J_size,1:K,1:J_size]<=inst.D)
    # duals
    @variable(part_model, μ_1[1:J_size, 1:K, 1:J_size, 1:J_size]>= 0)
    @variable(part_model, μ_2[1:J_size, 1:K]>=0)
    @variable(part_model, μ_3[1:J_size,1:K, 1:K]>= 0)


    # upper-level constraints ----------------------

    # z_obj is worst of the K cell objectives
    @constraint(part_model, [k=1:K], z[k] <= z_obj)
    # define cell objectives
    @constraint(part_model, [k=1:K], z[k] >= p*sum(s[j,k] for j=1:m)+sum(c[i,j]*y[i,j,k] for i=1:I_size for j=1:J_size))
    # allocation must comply with supply storage
    @constraint(part_model, [i=1:I_size,k=1:K], sum(y[i,j,k] for j=1:J_size) <=x[i])
    # bounded supply
    @constraint(part_model, sum(x[i] for i=1:I_size)<=inst.agg_supply_bound)
    # demand must be satisfied
    @constraint(part_model, [j=1:J_size, k=1:K], sum(y[i,j,k] for i=1:I_size)+s[j,k] >= ξ[j,k,j])

    # lower-level constraints as KKT--------------------------
    f_diff = LinearAlgebra.UniformScaling(-1)
    
    g_1(j1, j2, ξ_1) = ξ_1[j1] - ξ_1[j2] - norm(loc_J[j1,:]-loc_J[j2,:],Inf)

    g_2(ξ_2) = sum(ξ_2[j] for j in 1:J_size) - inst.cont_perc*inst.D*J_size
    g_2_diff = ones(m)

    g_3(l,k,ξ_3) = (v[l].-v[k])*ξ_3.-0.5(v[l].-v[k])*(v[l]+v[k])
    g_3_diff(l,k) = v[l]-v[k]

    # stationary
    @constraint(part_model, [j=1:J_size, k=1:K], f_diff[j,:].+μ_2[j,k]g_2_diff.+sum(μ_3[j,k,l]g_3_diff(l,k) for l = 1:K) .== 0)

    # primal feasibility
    @constraint(part_model, [j=1:J_size, k=1:K, j1=1:J_size, j2=1:J_size], g_1(j1, j2, ξ[j,k]) <= 0)
    @constraint(part_model, [j=1:J_size, k=1:K], g_2(ξ[j,k]) <= 0)
    @constraint(part_model, [j=1:J_size, k=1:K, l=1:K], g_3(l,k,ξ[j,k]) <= 0)

    # complementary slack
    F(j,k) = [g_1.(1:J_size,1:J_size,ξ[j,k]),g_2(ξ[j,k]), g_3.(1:K, k, ξ[j,k])]
    M(j,k) = [μ_1[j,k, 1:J_size,1:J_size], μ_2[j,k], μ_3[j,k,1:K]]
    @constraint(part_model, [j=1:J_size,k=1:K], F(j,k) ⟂ M(j,k))

    # @constraint(part_model, [j=1:J_size, k=1:K], sum(μ_1[j,k,j1,j2]g_1(j1,j2,ξ[j,k]) for j1 in 1:J_size for j2 in 1:J_size) +
    #                                             μ_2[j,k]g_2(ξ[j,k]) +
    #                                             sum(μ_3[j,k,l]g_3(l,k,ξ[j,k]) for l in 1:J_size) == 0 )

    optimize!(part_model)

    return v,x,y
end

