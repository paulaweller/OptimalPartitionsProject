using JuMP, Gurobi, LinearAlgebra, Dualization

    J_size = 2
    K = 2
    v = Matrix{Int64}(undef, J_size, K)

    sub_model = Model(Gurobi.Optimizer)

    @variable(sub_model, ξ_sub[1:J_size])

    # gradient of objective j is -1*j'th unit vector
    #f_diff = -1(1:J_size .== j)
    
    g_1(j1, j2, ξ_1) = ξ_1[j1] - ξ_1[j2] - 10

    g_2(ξ_2) = sum(ξ_2[j] for j in 1:J_size) - 1

    g_3(l,k,ξ_3) = dot(v[:,l].-v[:,k],ξ_3) -0.5*dot(v[:,l].-v[:,k],v[:,l].+v[:,k])

    # primal feasibility
    @constraint(sub_model, g1_con[j=1:J_size, k=1:K, j1=1:J_size, j2=vcat([1:j1-1]...,[j1+1:J_size]...)], g_1(j1, j2, ξ_sub) <= 0)
    @constraint(sub_model, g2_con[j=1:J_size, k=1:K], g_2(ξ_sub) <= 0)
    @constraint(sub_model, g3_con[j=1:J_size, k=1:K, l=vcat([1:k-1]...,[k+1:K]...)], g_3(l,k,ξ_sub) <= 0)

    dual_model = dualize(sub_model; dual_names = DualNames("dual", "dual"))

    open("dualmodel.txt","a") do io
        println(io,dual_model)
     end