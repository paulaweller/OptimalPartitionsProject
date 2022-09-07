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

"""
    enum_uncset(D, loc_J, pc)

Construct the integerscenarios in the uncertainty set for maximal demand D and affected percentage pc.
"""
function enum_uncset(D, loc_J, pc)

    # number of demand points
    J = size(loc_J, 1)

    # demand scenarios that satisfy the constraints
    uncset_enum = []

    # all possible demand configurations, i.e. {0, ..., D}^J
    it = Iterators.repeated(0:D,J)

    for d in Iterators.product(it...)

        # d must not affect more than the given percentage
        if sum(d) <= pc*D*J

            # d is a candidate
            add = true
            # every pair of demand points must be considered
            for (j1,j2) in Iterators.product(1:J,1:J)

                # if constraint is violated, don't add this demand scenario
                if (d[j1]-d[j2] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))  == false
                add = false
                end
            end
        else 
            add = false
        end
        if add == true
            push!(uncset_enum, d)
        end
    end
    return uncset_enum
end

"""
    find_plan(d_real, q, c)

Find the best of the plans in q for this uncertainty realization d_real.
"""
function find_plan(d_real, q, c)

    obj_val = 10^10
    best_q = []
    for p in 1:size(q, 3)
        obj_val_it = 100*sum(max(0, d_real[j]-sum(q[i,j,p] for i in 1:size(q,1))) for j in 1:size(q,2)) + sum(c[i,j]*q[i,j,p] for i in 1:size(q,1), j in 1:size(q,2))
        if obj_val_it < obj_val
            obj_val = obj_val_it
            best_q = q[:,:,p]
        end
    end
    return obj_val, best_q
end

"""
    obs_obj(uncertainty_set, find_plan, q_all, c)

Calculate actually observable objectives.
"""
function obs_obj(uncertainty_set, find_plan, q_all, c)

    obs_objectives = []

    # for every demand scenario
    for u in uncertainty_set
        # find the best plan
        obs_obj, obs_q = find_plan(u, q_all, c)
        # add it to the list
        push!(obs_objectives, obs_obj)
    end
    return obs_objectives
end

function box_plot_from_files(filenames, labelnames, n_datapoints)

    plot_data_obj = zeros(n_datapoints,length(filenames))
    count = 1
    for o_file in filenames
        # extract objective values from file
        o_values = open(o_file) do file
            obj = []
            for ln in eachline(file)
                val = parse(Float64, ln)
                push!(obj, val)
            end
        end
        plot_data_obj[:, count] =  o_values
        count = count+1
    end
    obj_plot = plot(xlabel="", ylabel="objective")
    boxplot!(obj_plot, labelnames, plot_data_obj, leg=false, linewidth=2,colour = [RGB(239/255, 53/255, 59/255) RGB(0/255, 94/255, 184/255) RGB(255/255, 205/255, 0/255) RGB(210/255, 22/255, 53/255)],linecolour= :match,fillalpha = 0.4)
    savefig(obj_plot, "results/boxplot_obj.pdf")
end

function box_plot(valuevector, labelnames)
    plot_data_obj = valuevector'
    obj_plot = plot(xlabel="", ylabel="objective")
    boxplot!(obj_plot, labelnames, plot_data_obj, leg=false, linewidth=2,colour = [RGB(239/255, 53/255, 59/255) RGB(0/255, 94/255, 184/255)],linecolour= :match,fillalpha = 0.4)
    savefig(obj_plot, "results/boxplot_obj.pdf")
end