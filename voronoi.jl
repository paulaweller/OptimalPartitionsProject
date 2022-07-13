using JuMP, Gurobi, BilevelJuMP
include("helpers.jl")


"""
    solve_voronoi(K, inst)

    Solve the K-adaptable pre-allocation problem for instance "inst" with Voronoi-style partitioning.
"""

function solve_voronoi(K, inst)

    I = size(inst.loc_I, 1)
    J = size(inst.loc_J, 1)
    c = reshape([norm(inst.loc_I[i,:]-inst.loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    

    part_model = BilevelModel(Gurobi.Optimizer, mode=BilevelJuMP.SOS1Mode())

    # variables ---------------------------------

    # supply
    @variable(Upper(part_model), 0<=x[1:I])
    # allocation of supplies
    @variable(Upper(part_model), 0<=y[1:I,1:J,1:K])
    # unsatisfied demand
    @variable(Upper(part_model), 0<=s[1:J,1:K])
    # objectives of cells
    @variable(Upper(part_model), z[1:K])
    # overall objective
    @variable(Upper(part_model), z_obj)
    # uncertain demand
    @variable(Lower(part_model), ξ[1:J,1:K]<=inst.D)
    # voronoi points
    @variable(Upper(part_model), v[1:K])

    # upper-level constraints ----------------------

    # z_obj is worst of the K cell objectives
    @constraint(Upper(part_model), [k=1:K], z[k] <= z_obj)
    # define cell objectives
    @constraint(Upper(part_model), [k=1:K], z[k] >= p*sum(s[j,k] for j=1:m)+sum(c[i,j]*y[i,j,k] for i=1:I for j=1:J))

    # allocation must comply with supply storage
    @constraint(Upper(part_model), [i=1:I,k=1:K], sum(y[i,j,k] for j=1:J) <=x[i])

    # bounded supply
    @constraint(Upper(part_model), sum(x[i] for i=1:I)<=inst.agg_supply_bound)

    # demand must be satisfied
    @constraint(Upper(part_model), [j=1:J, k=1:K], sum(y[i,j,k] for i=1:I)+s[j,k] >= ξ[j,k])

    # lower-level constraints--------------------------

    # aggregated demand is bounded
    @constraint(Lower(part_model), sum(ξ[j,k] for j in 1:J) <= round(Int, inst.cont_perc*inst.D*J))
    # for each pair of demand points, add constraint that if locations are close, demand values must be close, too
    for j1 in 1:J
        for j2 in j1+1:J
            @constraint(Lower(part_model), ξ[j1]-ξ[j2] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))
            @constraint(Lower(part_model), ξ[j2]-ξ[j1] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))
        end
    end

    @constraint(Lower(part_model), [k=1:K,l=1:K], (v[l].-v[k]).*ξ[j,k] <= (v[l].-v[k]).*0.5.*(v[l].+v[k]))
end

