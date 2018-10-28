#=====
BAF complex structure inference using a genetic algorithm
=====#

doc = """Run genetic algorithm to infer BAF complex structure.

Usage:
  pbaf_runSim.jl [--threshold=<alpha>] [--population=<N>] [--generations=<L>] [--rate=<P>] [--monitor=<step>] [--output=<file>]
  pbaf_runSim.jl -h | --help
  pbaf_runSim.jl --version

Options:
  -h --help            Show this screen.
  --version            Show version.
  --threshold=<alpha>  Cutoff for significance of enrichment [default: 0.05].
  --population=<N>     Number of graphs in population [default: 25].
  --generations=<L>    Number of iterations of the algorithm [default: 200].
  --rate=<P>           Probability of having a mutation for each graph, for each mutation type, at every generation [default: 0.026].
  --monitor=<step>     Time step at which the performance of the population should be monitored [default: 50].
  --output=<file>.     Name of the output file [default: outputRun].
"""

using DocOpt, CSV, DataFrames, StatsBase, LightGraphs, Distributions, JLD, HDF5

args = docopt(doc, version=v"1.0.0")

#=== Parameters ===#
global const ALPHA = parse(Float32, args["--threshold"])

# Run parameters
const N = parse(Int32, args["--population"]) # Number of graphs [500, 1000]
const L = parse(Int32, args["--generations"]) # Number 0f iterations [minimum 2000/1000 needed, 5000, 10000,25000]
const P = parse(Float32, args["--rate"]) # Probability of mutation [0.01275, 0.026]
# Expect 10% of graphs mutated per generation

# How often should we keep track of the system's state?
monitorStep = parse(Int32, args["--monitor"])

# Where to output the results of this run
outputFileName = String(args["--output"])

# Prepare data
colnames = ["Units"    
 "ACTB"     
 "ARID1A.10"
 "ARID1A.3" 
 "ARID1B"   
 "ARID2"    
 "BCL11A"   
 "BCL11B"   
 "BCL7A"    
 "BCL7B"    
 "BRD7"     
 "BRD9"     
 "DPF1"     
 "DPF2"     
 "DPF3"     
 "PBRM1"    
 "PHF10"    
 "SMARCA2"  
 "SMARCC1"  
 "SMARCC2"  
 "SMARCD1"  
 "SMARCD2"  
 "SMARCD3"]
brgData = CSV.read("BRG1-data.csv"; delim='\t', header=colnames, datarow=2)
foreach(x -> brgData[x] = log2.(brgData[x]), names(brgData[:,2:end]))

brgPval = CSV.read("BRG1-pval.csv"; delim='\t', header=colnames, datarow=2)
brgPval[1] = brgData[1]

# We now remove variations where the fold change is not significantly greater than zero.
for i in 2:length(brgData)
    for j in 1:length(brgData[i])
        if brgPval[j,i] > ALPHA
            brgData[j,i] = 0
        end
    end
end

# BAF and PBAF complex structures

# Join ARID1A.10 and ARID1A.3
delete!(brgData, Symbol("ARID1A.3"))
rename!(brgData, Symbol("ARID1A.10") => :ARID1A)

# Create interaction graph
studyBAFko = convert(Array{String,1}, names(brgData[2:end]))
studyBAFpd = convert(Array{String,1}, brgData[1])
unitDict = Dict(s => i for (i,s) in enumerate(sort(union(studyBAFko, studyBAFpd))))

effectGraph = SimpleDiGraph()
add_vertices!(effectGraph, length(unitDict))

# Store the sign of the log2-fold-change associated with each link
edgeTypes = Dict{Tuple, String}() 

# Parse each column
for x = names(brgData[:,2:end])
    for y = 1:length(brgData[x])
        if brgData[y,x] < 0
            add_edge!(effectGraph, unitDict[String(x)], unitDict[String(brgData[y,:Units])])
            edgeTypes[(unitDict[String(x)], unitDict[String(brgData[y,:Units])])] = "inhibits"
        elseif brgData[y,x] > 0
            add_edge!(effectGraph, unitDict[String(x)], unitDict[String(brgData[y,:Units])])
            edgeTypes[(unitDict[String(x)], unitDict[String(brgData[y,:Units])])] = "enhances"
        end
    end
    println()
end


# Which elements should we include in our structural model?
studyBAFunits = [k for k in keys(unitDict)]

# How many subunits are we considering?
const M = length(studyBAFunits)

# Define graph Julia struct
# Pulldown graphs contain the directed graph of activation/inhibition, the node and edges annotations.  
# Structure graphs contain the structural graph, the node annotations and the competition classes of each node.

mutable struct pulldownGraph
    graph::SimpleDiGraph
    nodes::Dict{Int64,String}
    edges::Dict{Tuple, String}
end

mutable struct structureGraph
    graph::SimpleGraph
    nodes::Dict{Int64,String}
    competition::Dict{Int64,Int64}
end

# Define constants used by the algorithm

const inhibitEdge = "inhibits"
const enhanceEdge = "enhances"

# Remember SMARCA4 index
const brgIndex = [i for i in 1:length(studyBAFunits) if studyBAFunits[i] == "SMARCA4"][1]
# SMARCA4 should not be the last subunit in the

# Link indices to unsorted list of BAF units
unitDictStudy = Dict(enumerate(studyBAFunits))
# Convert node indices from experimental graph to simulated graphs
convertUnitIndex = Dict(unitDict[v] => u for (u,v) in unitDictStudy)
observedEdges = Dict((convertUnitIndex[u[1]], convertUnitIndex[u[2]]) => v for 
        (u,v) in edgeTypes if u[1] in keys(convertUnitIndex) && u[2] in keys(convertUnitIndex))

studyBAFpdIndices = [ipd for (ipd, pd) in enumerate(studyBAFunits) if pd in studyBAFpd]
studyBAFkoIndices = [iko for (iko, ko) in enumerate(studyBAFunits) if ko in studyBAFko]

# Define graph functions

"""
Compute pulldown graph corresponding to a
structure graph given as argument
"""
function structureToPulldown(sGraph::structureGraph)
    # The structure graph must include all BAF subunits
    # @assert nv(sGraph.graph) == length(studyBAFunits)
    
    # Initialise a pulldownGraph
    # with the studied nodes and no edges
    pGraph = pulldownGraph(
        SimpleDiGraph(M),
        sGraph.nodes,
        Dict{Tuple, String}()
    )
    
    # Create dict from competitions between units                
    competitionDict = getCompetitionDict(sGraph.competition)
                    
    # For each unit knocked-out
    for iko = studyBAFkoIndices
        # Compute what units are still connected to SMARCA4
        pulledComponent = getPulledComponent(sGraph.graph, iko)
        
        # Check what would be observed for each pulled down subunit
        for ipd = studyBAFpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraph, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledComponent)
                        add_pulldown_edge!(inhibitEdge, pGraph, iko, ipd)
                        continue # Look at next pulldowned subunit
                    end
                    # The PD subunit is connected
                    if enhanceIfDisconnectedCompetition!(pGraph, pulledComponent,
                        competitionDict, ipd, iko)
                        # The subunit is enriched
                        continue # Look at next pulldowned subunit
                    end
                elseif !(ipd in pulledComponent)
                    # If a subunit is not in the component connected
                    # to SMARCA4, the KO will decrease the quantity of
                    # this subunit that will be pulled-down
                    add_pulldown_edge!(inhibitEdge, pGraph, iko, ipd)
                    continue # Look at next pulldowned subunit
                else
                    # The PD subunit is connected
                    enhanceIfDisconnectedCompetition!(pGraph, pulledComponent,
                        competitionDict, ipd, iko)
                    continue # Look at next pulldowned subunit
                end
            end
        end        
    end
    
    return(pGraph)
end

"""
Add a link to a pulldownGraph
"""
function add_pulldown_edge!(edgeType::String, pGraph::pulldownGraph, from::Int64, to = from)
    add_edge!(pGraph.graph, from, to)
    pGraph.edges[(from, to)] = edgeType
end
                        
"""
Create a dictionary associating a subunit with its competitors
"""
function getCompetitionDict(competition::Dict{Int64,Int64})
    competitionDF = DataFrame(Int64, M, 2)
    for i in 1:M
        competitionDF[i,1] = i
        competitionDF[i,2] = competition[i]
    end
    names!(competitionDF, [:Key, :Value])
    
    competitionDict = Dict{Int64, Array}()
    for df in groupby(competitionDF, :Value)
        for value in df[:Key]
            competitionDict[value] = [i for i in df[:Key] if i != value]
        end
    end
    
    return(competitionDict)
end

"""
Predict enrichment if a KO disconnect a competitor
of a subunit
"""
function enhanceIfDisconnectedCompetition!(pGraph::pulldownGraph, 
        pulledComponent::Array{Int64,1}, competitionDict::Dict{Int64, Array},
        ipd::Int64, iko::Int64)
    # For the KOed subunit
    if ipd in competitionDict[iko]
        add_pulldown_edge!(enhanceEdge, pGraph, iko, ipd)
        return(true) # An edge has been added
    end    
    # For all non-KOed subunit
    for inc = (j for j in 1:(M-1) if !(j in pulledComponent))
        if inc == iko
            # If the subunit has the index 'iko' it is
            # actually the last subunit, that has been
            # swapped with the KOed subunit
            inc = M
        end
        if ipd in competitionDict[inc]
            add_pulldown_edge!(enhanceEdge, pGraph, iko, ipd)
            return(true) # An edge has been added
        end
    end
    return(false) # No edge has been added
end

"""
Return a list of all subunits still connected
to SMARCA4 after a given KO is performed
"""        
function getPulledComponent(graph::LightGraphs.SimpleGraphs.SimpleGraph{Int64}, iko::Int64)
    perturbGraph = copy(graph)
    rem_vertex!(perturbGraph, iko)
    pulledComponent = Array{Int64,1}
    for component in connected_components(perturbGraph) if brgIndex in component
        return(component)
    end end
end
                        
"""
Enforce the connectivity of a structureGraph
"""
function connectGraph!(sGraph::structureGraph)
    while !is_connected(sGraph.graph)
        mutateAddEdge!(sGraph)
end end
                        
"""
Attribute random competition classes for subunits not
yet present in competition dictionary of a structureGraph
"""
function randomCompetitionGraph!(sGraph::structureGraph)
    graph = sGraph.graph
    competition = sGraph.competition
    
    for i = 1:M
        # Get the index of all subunit not in the competition dict
        u = map(reverse, unitDictStudy)[studyBAFunits[i]]
        if !(u in keys(competition))
            # Assign random competition class
            competition[u] = rand(1:M)
        end
    end
end

# Define mutation functions

"""
Mutate a single structure graph
The keywords contain the mutation parameters:
    p_add: add edge probability
    p_del: del edge probability
    p_swp: swap edge probability
    p_cmp: competition class probability
"""
function mutateStructureGraph!(sGraph::structureGraph; 
        p_add = 0.1, p_del = p_add, p_swp = p_add, p_cmp = p_add)
    # Store exit codes of individual mutation functions
    status = 0
    
    # Determine which mutations to perform
    doMutate = rand(4) .< [p_add, p_del, p_swp, p_cmp]
    
    if doMutate[1]
        status += mutateAddEdge!(sGraph)
    end

    if doMutate[2]
        status += mutateDelEdge!(sGraph.graph)
    end

    if doMutate[3]
        status += mutateSwapEdges!(sGraph)
    end

    if doMutate[4]
        status += mutateCompetitors!(sGraph)
    end

    return(status)
end
  
"""
Add an edge to a structure graph
"""
function mutateAddEdge!(sGraph::structureGraph)
    graph = sGraph.graph
    competition = sGraph.competition
    N = nv(graph)
    
    if ne(graph) >= N*(N-1)/2
        # The graph is already complete
        return(1)
    else
        while true
            (a,b) = ceil.(N*rand(2))
            if (a != b) && (add_edge!(graph, a, b))
                # Do not allow self loop
                # Do not allow links between competitors
                if (competition[a] == competition[b])
                    rem_edge!(graph, Int64(a), Int64(b))
                    return(1)
                end
                # Exit if edge sucessfully added
                return(0)
            end
        end
    end
end

"""
Remove an edge to a structure graph
"""
function mutateDelEdge!(graph::LightGraphs.SimpleGraphs.SimpleGraph)
    edgesList = [e for e in edges(graph)]
    edgesIndicesOrder = randperm(length(edgesList))
    for edgeIndex in edgesIndicesOrder
        edgeToRemove = edgesList[edgeIndex]
        rem_edge!(graph, edgeToRemove)
        if is_connected(graph)
            return(0)
        else
            # So structure graph should be kept connected
            # Therefore we put back in the removed edge
            add_edge!(graph, edgeToRemove)
        end
    end
    
    # No edge can be removed without diconnecting the graph
    return(1)
end

"""
Swap edges in a structure graph
"""
function mutateSwapEdges!(sGraph::structureGraph)
    graph = sGraph.graph
    competition = sGraph.competition
    
    edgesList = [e for e in edges(graph)]
    edgesIndicesOrder = randperm(length(edgesList))
    
    for (indexIndex, edgeIndex) = enumerate(edgesIndicesOrder)
        edge1 = edgesList[edgeIndex]
        edge2 = edgesList[edgesIndicesOrder[1+(indexIndex % length(edgesList))]]
        # Ensure that no self link will be created
        if Tuple(edge1)[1] != Tuple(edge2)[2] && Tuple(edge2)[1] != Tuple(edge1)[2]
            # Start by deleting the old edges
            rem_edge!(graph, edge1)
            rem_edge!(graph, edge2)
            # Then add the new ones if not linking competitors
            if competition[Tuple(edge1)[1]] != competition[Tuple(edge2)[2]]
                add_edge!(graph, Tuple(edge1)[1], Tuple(edge2)[2])
            end
            if competition[Tuple(edge2)[2]] != competition[Tuple(edge1)[2]]
                add_edge!(graph, Tuple(edge2)[1], Tuple(edge1)[2])
            end
            if is_connected(graph)
                return(0)
            else
                # So structure graph should be kept connected
                # Therefore we put back in the removed edges
                add_edge!(graph, edge1)
                add_edge!(graph, edge2)
                # NB: extra edges will stay if any
            end
        end
    end
    
    # No edges can be swapped without diconnecting the graph
    return(1)
end

"""
Mutate competing nodes
"""
function mutateCompetitors!(sGraph::structureGraph)
    graph = sGraph.graph
    competition = sGraph.competition
    
    # Select node to change competition class
    nodeComp = rand(1:nv(graph))
    # Select new competition class
    newComp = rand(1:nv(graph))
    for n = neighbors(graph, nodeComp)
        if competition[n] == newComp
            # Changing the competition class would lead to linked competitors
            return(1)
        end
    end
    competition[nodeComp] = newComp
    
    return(0)
end

# Genetic algorithm module

"""
Compute loss for a given structure
compared to observation
"""
function observedLoss(sGraph::structureGraph,
    details::Bool = false)
    pGraph = structureToPulldown(sGraph)
    
    intersectEdges = intersect(pGraph.edges, observedEdges)
    unionEdges = union(pGraph.edges, observedEdges)
    
    if details
        # Return array with Jaccard index
        # length of union and length of  
        return([length(intersectEdges) / length(unionEdges), length(intersectEdges), length(pGraph.edges)])
    else
        # Return Jaccard index
        return([length(intersectEdges) / length(unionEdges)])
    end
end

"""
Generate in place the new generation of 
structure graphs based on their fitness.
Return the fitness array.
"""
function reproduceGeneration!(pop::Array{structureGraph,1},
    details::Bool = false)
    jaccard = map(x -> observedLoss(x,details), pop)
    fitness = map(x -> x[1], jaccard)
    fitness ./= sum(fitness)
    
    sumFitness = sum(fitness) 
    if sumFitness != 1
        fitness[end] += 1 - sumFitness
    end
    # Ensure the cumulative fitnesses is a probability distribution
    
    offspringPerGraph = rand(Multinomial(length(pop), fitness), 1)
    offspring = Array{structureGraph,1}(length(pop))
    
    offspringToFill = 1 # Which is the next index to be filled?
    for (ipop, noff) = enumerate(offspringPerGraph)
        for ioff = 1:noff
            offspring[offspringToFill] = deepcopy(pop[ipop])
            offspringToFill += 1
        end
    end
    
    # Ensure the best structure graph is kept
    bestGraphIndex = findmax(fitness)[2]
    if offspringPerGraph[bestGraphIndex] == 0
        # No offspring for the best graph
        # So we force one
        offspring[1] = deepcopy(pop[bestGraphIndex])
    end
    
    pop .= offspring
        
    return(jaccard)
end

"""
Generate the new generation of structure networks
"""
function newGeneration!(pop::Array{structureGraph,1},
        details::Bool = false;
        p_add = 0.1, p_del = p_add, p_swp = p_add, p_cmp = p_add, p_crs = p_add/10)
    # Fitness-based reproduction
    fitness = reproduceGeneration!(pop, details)
    
    # Mutate potentially each structure network
    map(x -> mutateStructureGraph!(x;
            p_add = p_add, p_del = p_del, p_swp = p_swp, p_cmp = p_cmp), pop)
    
    return(fitness)
end


# Run algorithm

# Array containing the graphs
pop = Array{structureGraph,1}(N)

# Max number of edges in a graph
maxEdges = Int64(M*(M-1)/2)

# Unit names dictionary
unitDict = Dict(e => studyBAFunits[e] for e in 1:M)

# Initialize population
map!(x -> structureGraph(
        Graph(M, rand(1:maxEdges)),
        unitDict,
        Dict(e => e for e in 1:M)),
    pop)

# Ensure connectivity
map(connectGraph!, pop)

# Iterate and store performance
quantileFitness = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileIntersect = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileSimulatedEdges = Array{Float16}(Int(ceil(L/monitorStep)), 5)
for i in 1:L
    if i % monitorStep == 1
        f = newGeneration!(pop, true, p_add = P)
        quantileFitness[Int(ceil(i/monitorStep)),:] = quantile(map(x -> x[1], f))
        quantileIntersect[Int(ceil(i/monitorStep)),:] = quantile(map(x -> x[2], f))
        quantileSimulatedEdges[Int(ceil(i/monitorStep)),:] = quantile(map(x -> x[3], f))
    else
        f = newGeneration!(pop, false, p_add = P)
    end
end

save(outputFileName*".jld","pop", pop,
    "fitness", quantileFitness, "intersect", quantileIntersect, "quantileSimulatedEdges", quantileSimulatedEdges)