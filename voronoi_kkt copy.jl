using JuMP, GAMS, LinearAlgebra
include("helpers.jl")


"""
    solve_kkt(K, inst)

    Solve the KKT conditions for the K-adaptable pre-allocation problem for instance "inst".
"""

function solve_voronoi_kkt(K, inst)

    I_size = size(inst.loc_I, 1)
    J_size = size(inst.loc_J, 1)
    @show c = reshape([norm(inst.loc_I[i,:]-inst.loc_J[j,:]) for j in 1:J_size for i in 1:I_size],I_size,J_size)
    

    part_model = Model(GAMS.Optimizer)
    set_optimizer_attribute(part_model, "SysDir", "C:\\GAMS\\37")
    # variables ---------------------------------

    # supply
    @variable(part_model, 0<=x[1:I_size])
    # allocation of supplies
    @variable(part_model, 0<=y[1:I_size,1:J_size,1:K])
    # unsatisfied demand
    @variable(part_model, 0<=s[1:J_size,1:K])
    # objectives of cells
    @variable(part_model, z[1:K]>= 0)
    # overall objective
    @variable(part_model, z_obj>= 0)
    # voronoi points
    @variable(part_model, 0<=v[1:J_size,1:K]<=inst.D)

    # lower-level variables
    # worst-case demand scenario, one per partition k and row j
    @variable(part_model, 0<= ξ[1:J_size,1:J_size,1:K]<=inst.D)
    # duals
    # @variable(part_model, μ_1[1:J_size, 1:K, 1:J_size, 1:J_size]>= 0)
    @variable(part_model, μ_2[1:J_size, 1:K]>=0)
    @variable(part_model, μ_3[1:J_size,1:K, 1:K]>= 0)


    # upper-level constraints ----------------------

    # z_obj is worst of the K cell objectives
    @constraint(part_model, [k=1:K], z[k] <= z_obj)
    # define cell objectives
    @constraint(part_model, [k=1:K], z[k] >= 100*sum(s[j,k] for j=1:m)+sum(c[i,j]*y[i,j,k] for i=1:I_size for j=1:J_size))
    # allocation must comply with supply storage
    @constraint(part_model, [i=1:I_size,k=1:K], sum(y[i,j,k] for j=1:J_size) <=x[i])
    # bounded supply
    @constraint(part_model, sum(x[i] for i=1:I_size)<= inst.W)
    # demand must be satisfied
    @constraint(part_model, [j=1:J_size, k=1:K], sum(y[i,j,k] for i=1:I_size)+s[j,k] >= ξ[j,j,k])

    # lower-level constraints as KKT--------------------------

    # gradient of objective j is -1*j'th unit vector
    #f_diff = -1(1:J_size .== j)
    
    # g_1(j1, j2, ξ_1) = ξ_1[j1] - ξ_1[j2] - norm(inst.loc_J[j1,:]-inst.loc_J[j2,:],Inf)
    # g_1_diff(j1,j2) = (1:J_size.== j1) .+ -1*(1:J_size.== j2)

    g_2(ξ_2) = sum(ξ_2[j] for j in 1:J_size) - inst.pc*inst.D*J_size
    g_2_diff = ones(m)

    g_3(l,k,ξ_3) = dot(v[:,l].-v[:,k],ξ_3) -0.5*dot(v[:,l].-v[:,k],v[:,l].+v[:,k])
    g_3_diff(l,k) = v[:,l].-v[:,k]

    # stationary
    @constraint(part_model, [j=1:J_size, k=1:K], -1*(1:J_size.== j)#.+sum(μ_1[j,k,j1,j2]*g_1_diff(j1,j2) for j1=1:J_size for j2=vcat([1:j1-1]...,[j1+1:J_size]...)) 
                                                .+ μ_2[j,k]*g_2_diff.+sum(μ_3[j,k,l]*g_3_diff(l,k) for l=vcat([1:k-1]...,[k+1:K]...)) .== 0)

    # primal feasibility
    # @constraint(part_model, [j=1:J_size, k=1:K, j1=1:J_size, j2=vcat([1:j1-1]...,[j1+1:J_size]...)], g_1(j1, j2, ξ[j,:,k]) <= 0)
    @constraint(part_model, [j=1:J_size, k=1:K], g_2(ξ[j,:,k]) <= 0)
    @constraint(part_model, [j=1:J_size, k=1:K, l=vcat([1:k-1]...,[k+1:K]...)], g_3(l,k,ξ[j,:,k]) <= 0)

    # complementary slack

    # @constraint(part_model, [j=1:J_size,k=1:K,j1=1:J_size,j2=vcat([1:j1-1]...,[j1+1:J_size]...)], g_1(j1,j2,ξ[j,:,k])⟂ μ_1[j,k,j1,j2])
    @constraint(part_model, [j=1:J_size,k=1:K], g_2(ξ[j,:,k])⟂ μ_2[j,k])
    @constraint(part_model, [j=1:J_size,k=1:K,l=vcat([1:k-1]...,[k+1:K]...)], g_3(l, k, ξ[j,:,k])⟂ μ_3[j,k,l])

    # @constraint(part_model, [j=1:J_size, k=1:K], sum(μ_1[j,k,j1,j2]g_1(j1,j2,ξ[j,k,:]) for j1 in 1:J_size for j2 in 1:J_size) +
    #                                             μ_2[j,k]g_2(ξ[j,k,:]) +
    #                                             sum(μ_3[j,k,l]g_3(l,k,ξ[j,k,:]) for l in 1:J_size) == 0 )

    @objective(part_model, Min, z_obj)

    open("model.txt","a") do io
        println(io,part_model)
     end
    optimize!(part_model)

    return value.(v), value.(x), value.(y), value.(s), value.(ξ)
end

