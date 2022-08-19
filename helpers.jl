using LinearAlgebra, StatsPlots, Random


mutable struct AllocationInstance
    loc_I::Matrix{Int64}
    loc_J::Matrix{Int64}
    W::Int64
    D::Int64
    pc::Float64
end
#const GRB_ENV = Gurobi.Env()
"""
    generate_instance(I_inst, J_inst, seed, demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J), plot_loc=false)
Generate a problem instance with the given parameters. Returns: loc_I_inst, loc_J_inst, demand_bound, cont_perc, agg_supply_bound
Fields:
    I_inst              Number of service points
    J_inst              Number of demand points
    seed                For reproducibility
    demand_bound        Upper bound of demand at one demand point
    cont_perc           Percentage of damage caused by the contingency
    agg_supply_bound    Aggregated supply bound, default: maximal aggregated demand 
    plot_loc            Should the locations of service and demand points be plotted?
    loc_max             How large is the grid
"""
function generate_instance(I_inst, J_inst, seed; demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J_inst), plot_loc=false, loc_max=5)

    # seed for reproducibility
    Random.seed!(seed)
    
    # set of locations for service and demand points
    loc_set = []
    # locations of service points
    loc_I_inst = []
    # locations of demand points
    loc_J_inst = []

    # no two points can have the same location
    while length(loc_set) !== I_inst+J_inst  
        
        # randomly generate locations of service points as x,y-coordinates
        loc_I_inst = rand(1:loc_max,I_inst,2)
        # randomly generate locations of demand points as coordinates
        loc_J_inst = rand(1:loc_max,J_inst,2)
        # add to set (not list) of locations, i.e. a location won't be listed twice
        loc_set = union(union(loc_I_inst[i,:] for i in 1:I_inst),union(loc_J_inst[j,:] for j in 1:J_inst))
    end

    # sort locations by x-coordinate for simplicity
    loc_I_inst = sort(loc_I_inst, dims=1)
    loc_J_inst = sort(loc_J_inst, dims=1)

    instance = AllocationInstance(loc_I_inst,loc_J_inst,agg_supply_bound,demand_bound,cont_perc)
    # plot service and demand points
    if plot_loc == true

        scatter(loc_I_inst[:,1],loc_I_inst[:,2],
                    lims=[0,loc_max+1.2],
                    series_annotations = text.(1:J_inst, :inside),
                    markersize = 20,
                    lab="service point")
        display(scatter!(loc_J_inst[:,1],loc_J_inst[:,2], 
                            series_annotations = text.(1:J_inst, :inside),
                            shape = :square,
                            markersize = 15,
                            lab="demand point",
                            aspect_ratio=:equal
                            ))
    end

    return instance
end