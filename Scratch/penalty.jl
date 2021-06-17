#function my_incidence_matrix(G)
#    I = vcat([src(e) for e in edges(G)], [dst(e) for e in edges(G)])
#    J = vcat(collect(1:ne(G)), collect(1:ne(G)))
#    V = vcat(fill(-1, ne(G)), fill(1, ne(G)))
#    return sparse(I, J, V)
#end


"""
Solves mixed assignment where a set of links (given by bounded indices) are
prohibited for the UE vehicles

    `me_excluded_assignment_penalty(rn,
                                    ods,
                                    demands,
                                    γ,
                                    bounded_indices,
                                    r_array)` 
"""
function me_excluded_assignment_penalty(rn, ods, demands, γ, bounded_indices, r_array) 
    
    # Cost structure (standard)
    #A = my_incidence_matrix(rn.g)
    A = incidence_matrix(rn.g)
    a = rn.edge_params[:a]
    B = Diagonal(rn.edge_params[:b])
    
    # Vectors for conservation constraints (including OD terms)
    n_ods = length(ods)
    
    d_vects_hv = SparseVector{Float64,Int64}[]
    d_vects_av = SparseVector{Float64,Int64}[]
    for i in 1:n_ods
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]
        
        push!(d_vects_hv, (1-γ).*d_vec)
        push!(d_vects_av, γ.*d_vec)
    end
    
    # ME objective function (unpenalised)       
    me_obj_func(x,y) = (0.5*(x+y)'*B*(x+y) + a'*x + 0.5*a'*y)[1]
                             
    # Penalty function per link
    penalty_function(x, r, m) = r*(x^m)
    
    # penalty term for objective
    function pen_term(x, param, inds; power=2)
       ps = penalty_function.(x[inds], 0.5*a[inds], power)
       param*sum(ps)
    end
    
    # Penalised objective
    p_obj(x, y, p_param, b_indices) = me_obj_func(x,y) + 
                                      pen_term(x, p_param, b_indices)


    ### JuMP model construction
    
    #p_stap = Model(Ipopt.Optimizer)
    #p_stap = Model(Gurobi.Optimizer)
    p_stap = Model(() -> Gurobi.Optimizer(env))

    # OD specific link flows HV         
    @variable(p_stap, x[1:m,1:n_ods])
    @constraint(p_stap, non_neg_x, x .>= 0)

    # OD specific link flows AV
    @variable(p_stap, y[1:m,1:n_ods])
    @constraint(p_stap, non_neg_y, y .>= 0)

    #hv_flow = sum(x, dims=2)
    #av_flow = sum(y, dims=2)

    # Conservation constraints
    @constraint(p_stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(p_stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])


    # Unpenalised objective
    @objective(p_stap, Min, me_obj_func(x, y))
    
    ### 
    ### Solve unpenalised
    ###
    
    optimize!(p_stap)

    # Update
    sols_x = value.(x)
    sols_y = value.(y)

    iters_x = [sols_x[:]]
    iters_y = [sols_y[:]]
    
    ###
    ### Sequence of penalised sub-problems
    ###
    
    #r_array = [10, 100, 1000, 10000, 100000, 1000000]
    
    for r in r_array
        @objective(p_stap,
                   Min,
                   p_obj(x, y, r, bounded_indices))
        
        #for multiple ods another loop has to be added
        for k in 1:m
            set_start_value(x[k,1], sols_x[k])
            set_start_value(y[k,1], sols_y[k])
        end           
        
        optimize!(p_stap)
        
        sols_x = value.(x)
        sols_y = value.(y)

        push!(iters_x, sols_x[:])
        push!(iters_y, sols_y[:])
    end
    
    return iters_x, iters_y
    # Think about returning dual variables
end

"""
Run penalty method for network with restricted flow on some
links for a range of penetration rates
"""
function me_excluded_assignment_penalty_pr(rn, ods, demands, γ_array, bounded_indices, r_array)
    
    xx = []
    yy = []
    agg_a = []
    for γ in γ_array
        x, y = me_excluded_assignment_penalty(rn,
                                              ods,
                                              demands,
                                              γ,
                                              bounded_indices,
                                              r_array)
        push!(xx, x)
        push!(yy, y)
        push!(agg_a, x+y)
    end
    
    return xx, yy, agg_a
end

