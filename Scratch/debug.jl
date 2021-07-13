# Optimiser choice
optimiser = :gurobi

# Penalty parameters
penalty_steps = 8
r_array = [10^i for i in 1:penalty_steps]

# Penetration rate range
γ_init = 0.1
γ_step = 0.4
γ_stop = 0.9
# Penetration rate
γ_array = γ_init:γ_step:γ_stop

# Demand range
d_start = 0.005
d_step =  0.005
d_stop =  0.02
# d range
demand_range = d_start:d_step:d_stop

# MH parameters
μ = 4.0
p_splice = 0.75

# perhaps should be a function of the network size
# Maybe the same number of paths should be added?
num_metropolis_steps = 10000
stagger = 10

# Number of path samples
samples = 5

# Dirs and net files
ens_dir = "./TestData"
nn = 5

######################################################
######################################################

files_in_dir = readdir(ens_dir)
network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
net_files = sort(network_files, lt=natural)

# Use number of nets given, or all of them
if nn > 0 && nn < length(net_files)
    net_files = net_files[1:nn]
else
    parsed_args["num_nets"] = length(net_files)
end
N = length(net_files)


########################################################
########################################################

if optimiser == :gurobi
    using Gurobi
    # Use single licence token
    if !(@isdefined env)
        const env = Gurobi.Env()
    end
elseif optimiser == :ipopt
    using Ipopt
end

###########################################################
###########################################################

file = joinpath(ens_dir, net_files[1])

net = loadgraph(file, MGFormat()) # Undirected
rn = skel2rn(net) # Directed
a = rn.edge_params[:a]
b = rn.edge_params[:b]
sorted_edges = collect(edges(rn.g))
O = node_closest_to(net, [0,0]) # lower left-most node
D = node_closest_to(net, [1,1]) # upper right-most node

set_samples = excl_edge_set_sample(net,
                                   sorted_edges,
                                   O,
                                   D,
                                   μ,
                                   p_splice,
                                   weight_func,
                                   num_metropolis_steps,
                                   stagger,
                                   samples)
