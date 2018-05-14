#=====
BAF and PBAF complexes structures inference using a genetic algorithm
=====#

doc = """Run genetic algorithm to infer BAF complex structure.

Usage:
  bafpbaf_runSim.jl [--threshold=<alpha>] [--population=<N>] [--generations=<L>] [--rate=<P>] [--monitor=<step>] [--output=<file>]
  bafpbaf_runSim.jl -h | --help
  bafpbaf_runSim.jl --version

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
global const ALPHA = parse(Float64, args["--threshold"])

# Run parameters
const N = parse(Int32, args["--population"]) # Number of graphs [500, 1000]
const L = parse(Int32, args["--generations"]) # Number 0f iterations [minimum 2000/1000 needed, 5000, 10000,25000]
const P = parse(Float32, args["--rate"]) # Probability of mutation [0.01275, 0.026]
# Expect 10% of graphs mutated per generation

# How often should we keep track of the system's state?
monitorStep = parse(Int32, args["--monitor"])

# Where to output the results of this run
outputFileName = String(args["--output"])

colnames = ["Units"    
 "ACTB"     
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
 "SMARCA4.4"
 "SMARCA4.6"
 "SMARCC1"  
 "SMARCC2"  
 "SMARCD1"  
 "SMARCD2"  
 "SMARCD3"]
aridData = CSV.read("ARID1A-data.csv"; delim='\t', header=colnames, datarow=2)
foreach(x -> aridData[x] = log2.(aridData[x]), names(aridData[:,2:end]))

aridPval = CSV.read("ARID1A-pval.csv"; delim='\t', header=colnames, datarow=2)
aridPval[1] = aridData[1]

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

for i in 2:length(brgData)
    for j in 1:length(brgData[i])
        if brgPval[j,i] > ALPHA
            brgData[j,i] = 0
        end
    end
end

for i in 2:length(aridData)
    for j in 1:length(aridData[i])
        # Some values were stored as factors instead of floats, and could not be compared to ALPHA
        try
            if aridPval[j,i] > ALPHA
                aridData[j,i] = 0
            end
        catch e
            if isa(e, MethodError) # In case of type error when comparing the variable to ALPHA 
                if float(string(aridPval[j,i])) > ALPHA # Try converting the faulty variable
                    aridData[j,i] = 0
                end
            end
        end
    end
end

# Join SMARCA4.4 and SMARCA4.6
delete!(aridData, Symbol("SMARCA4.6"))
rename!(aridData, Symbol("SMARCA4.4") => :SMARCA4)

# Join ARID1A.10 and ARID1A.3
delete!(brgData, Symbol("ARID1A.3"))
rename!(brgData, Symbol("ARID1A.10") => :ARID1A)

# Define constants used by the algorithm
studyARIDko = convert(Array{String,1}, names(aridData[2:end]))
studyARIDpd = convert(Array{String,1}, aridData[1])
studyBRGko = convert(Array{String,1}, names(brgData[2:end]))
studyBRGpd = convert(Array{String,1}, brgData[1])
unitDictString = Dict(enumerate(sort(union(studyARIDko, studyARIDpd, studyBRGko, studyBRGpd))))
unitDict = map(reverse, unitDictString)
unitList = [string(u) for u in keys(unitDict)]
# How many subunits are we considering?
const M = length(unitDict)

edgeTypesARID = Dict{Float64, Dict{Tuple, String}}()

colnames = ["Units"    
 "ACTB"     
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
 "SMARCA4.4"
 "SMARCA4.6"
 "SMARCC1"  
 "SMARCC2"  
 "SMARCD1"  
 "SMARCD2"  
 "SMARCD3"]

for alpha = [0.1, 0.05, 0.01, 0.005]    
    aridData = CSV.read("ARID1A-data.csv"; delim='\t', header=colnames, datarow=2)
    
    foreach(x -> aridData[x] = log2.(aridData[x]), names(aridData[:,2:end]))

    for i in 2:length(aridData)
        for j in 1:length(aridData[i])
            # Some values were stored as factors instead of floats, and could not be compared to ALPHA
            try
                if aridPval[j,i] > alpha
                    aridData[j,i] = 0
                end
            catch e
                if isa(e, MethodError) # In case of type error when comparing the variable to ALPHA 
                    if float(string(aridPval[j,i])) > alpha # Try converting the faulty variable
                        aridData[j,i] = 0
                    end
                end
            end
        end
    end
    
    println(countnz(convert(Array,aridData)))

    # Join SMARCA4.4 and SMARCA4.6
    delete!(aridData, Symbol("SMARCA4.6"))
    rename!(aridData, Symbol("SMARCA4.4") => :SMARCA4)
    
    # Store the sign of the log2-fold-change associated with each link
    edgeTypesARID[alpha] = Dict{Tuple, String}()

    # Parse each column
    for x = names(aridData[:,2:end])
        for y = 1:length(aridData[x])
            if aridData[y,x] < 0
                edgeTypesARID[alpha][(unitDict[String(x)], unitDict[String(aridData[y,:Units])])] = "inhibits"
            elseif aridData[y,x] > 0
                edgeTypesARID[alpha][(unitDict[String(x)], unitDict[String(aridData[y,:Units])])] = "enhances"
            end
        end
    end
end

edgeTypesBRG = Dict{Float64, Dict{Tuple, String}}()

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

for alpha = [0.1, 0.05, 0.01, 0.005]    
    brgData = CSV.read("BRG1-data.csv"; delim='\t', header=colnames, datarow=2)
    
    foreach(x -> brgData[x] = log2.(brgData[x]), names(brgData[:,2:end]))

    for i in 2:length(brgData)
        for j in 1:length(brgData[i])
            # Some values were stored as factors instead of floats, and could not be compared to ALPHA
            try
                if brgPval[j,i] > alpha
                    brgData[j,i] = 0
                end
            catch e
                if isa(e, MethodError) # In case of type error when comparing the variable to ALPHA 
                    if float(string(brgPval[j,i])) > alpha # Try converting the faulty variable
                        brgData[j,i] = 0
                    end
                end
            end
        end
    end
    
    println(countnz(convert(Array,brgData)))
    
    # Join ARID1A.10 and ARID1A.3
    delete!(brgData, Symbol("ARID1A.3"))
    rename!(brgData, Symbol("ARID1A.10") => :ARID1A)
    
    # Store the sign of the log2-fold-change associated with each link
    edgeTypesBRG[alpha] = Dict{Tuple, String}()

    # Parse each column
    for x = names(brgData[:,2:end])
        for y = 1:length(brgData[x])
            if brgData[y,x] < 0
                edgeTypesBRG[alpha][(unitDict[String(x)], unitDict[String(brgData[y,:Units])])] = "inhibits"
            elseif brgData[y,x] > 0
                edgeTypesBRG[alpha][(unitDict[String(x)], unitDict[String(brgData[y,:Units])])] = "enhances"
            end
        end
    end
end

const inhibitEdge = "inhibits"
const enhanceEdge = "enhances"
# ARID1A and SMARCA4 should not be deleted in their respective immunoprecipitations
@assert !("ARID1A" in studyARIDko)
@assert !("SMARCA4" in studyBRGko)
    #=====
    When we delete a node from a lightgraph, the node to
    remove is swapped with the last node in the node list.
    To ensure that the index of ARID1A is stable, we make
    sur that it is never knocked out nor the last node.
    =====#
# Remember ARID1A and SMARCA4 indices
const brgIndex = [i for (i,u) in enumerate(unitList) if u == "SMARCA4"][1]
const aridIndex = [i for (i,u) in enumerate(unitList) if u == "ARID1A"][1]
# ARID1A and SMARCA4 should not be the last subunits in the lists
@assert aridIndex != M
@assert brgIndex != M

# Avoid to compute these sets in each iteration of the structureToPulldowns function
studyARIDpdIndices = [ipd for (ipd, pd) in unitDictString if pd in studyARIDpd]
studyARIDkoIndices = [iko for (iko, ko) in unitDictString if ko in studyARIDko]
studyBRGpdIndices = [ipd for (ipd, pd) in unitDictString if pd in studyBRGpd]
studyBRGkoIndices = [iko for (iko, ko) in unitDictString if ko in studyBRGko]
studyALLkoIndices = intersect(studyARIDkoIndices, studyBRGkoIndices)
studyALLpdIndices = intersect(studyARIDpdIndices, studyBRGpdIndices)
studyOnlyBRGpdIndices = setdiff(studyBRGpdIndices, studyARIDpdIndices)
studyOnlyARIDpdIndices = setdiff(studyARIDpdIndices, studyBRGpdIndices)
studyOnlyBRGkoIndices = setdiff(studyBRGkoIndices, studyARIDkoIndices)
studyOnlyARIDkoIndices = setdiff(studyARIDkoIndices, studyBRGkoIndices)

# Define graph Julia struct

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

# Define graph functions

"""
Compute pulldown graph corresponding to a
structure graph given as argument
"""
function structureToPulldowns(sGraph::structureGraph)
    # Initialise two pulldownGraphs
    # with the studied nodes and no edges
    pGraphARID = pulldownGraph(
        SimpleDiGraph(M),
        sGraph.nodes,
        Dict{Tuple, String}()
    )
    
    pGraphBRG = pulldownGraph(
        SimpleDiGraph(M),
        sGraph.nodes,
        Dict{Tuple, String}()
    )
    
    # Create dict from competitions between units                
    competitionDict = getCompetitionDict(sGraph.competition)
    
    # For units knocked-out in both precipitations
    for iko = studyALLkoIndices
        # Compute what units are still connected to ARID1A and SMARCA4
        pulledARIDComponent = getARIDPulledComponent(sGraph.graph, iko) 
        pulledBRGComponent = getBRGPulledComponent(sGraph.graph, iko) 
        
        # Check what would be observed for each pulled down subunit
        # Starting with the subunits pulled down in both precipitations
        for ipd = studyALLpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraphARID, ipd)
                add_pulldown_edge!(inhibitEdge, pGraphBRG, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledARIDComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                        competitionDict, ipd, iko)
                    end
                    if !(iko in pulledBRGComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                        competitionDict, ipd, iko)
                    end                    
                    continue # Look at next pulldowned subunit
                else
                    if !(ipd in pulledARIDComponent)
                        # If a subunit is not in the component connected
                        # to ARID1A, the KO will decrease the quantity of
                        # this subunit that will be pulled-down
                        add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                            competitionDict, ipd, iko)
                    end
                    if !(ipd in pulledBRGComponent)
                        # If a subunit is not in the component connected
                        # to ARID1A, the KO will decrease the quantity of
                        # this subunit that will be pulled-down
                        add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                            competitionDict, ipd, iko)
                    end
                end
            end
        end 
        
        # Then for ARID1A precipitation
        for ipd = studyOnlyARIDpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraphARID, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledARIDComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                        competitionDict, ipd, iko)
                    end               
                    continue # Look at next pulldowned subunit
                elseif !(ipd in pulledARIDComponent)
                    # If a subunit is not in the component connected
                    # to ARID1A, the KO will decrease the quantity of
                    # this subunit that will be pulled-down
                    add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                else
                    # The PD subunit is connected
                    enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                        competitionDict, ipd, iko)
                end
            end
        end
        
        # Finally for SMARCA4 precipitation
        for ipd = studyOnlyBRGpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraphBRG, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledBRGComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                        competitionDict, ipd, iko)
                    end               
                    continue # Look at next pulldowned subunit
                elseif !(ipd in pulledBRGComponent)
                    # If a subunit is not in the component connected
                    # to ARID1A, the KO will decrease the quantity of
                    # this subunit that will be pulled-down
                    add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                else
                    # The PD subunit is connected
                    enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                        competitionDict, ipd, iko)
                end
            end
        end
    end
    
    # Same for units knocked-out only in ARID1A precipitation
    for iko = studyOnlyARIDkoIndices
        # Compute what units are still connected to SMARCA4
        pulledARIDComponent = getARIDPulledComponent(sGraph.graph, iko) 
        
        # Check what would be observed for each pulled down subunit
        for ipd = studyARIDpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraphARID, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledARIDComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                        competitionDict, ipd, iko)
                    end                  
                    continue # Look at next pulldowned subunit
                else
                    if !(ipd in pulledARIDComponent)
                        # If a subunit is not in the component connected
                        # to ARID1A, the KO will decrease the quantity of
                        # this subunit that will be pulled-down
                        add_pulldown_edge!(inhibitEdge, pGraphARID, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphARID, pulledARIDComponent,
                            competitionDict, ipd, iko)
                    end
                end
            end
        end 
    end
    
    # Same for units knocked-out only in ARID1A precipitation
    for iko = studyOnlyBRGkoIndices
        # Compute what units are still connected to SMARCA4
        pulledBRGComponent = getBRGPulledComponent(sGraph.graph, iko) 
        
        # Check what would be observed for each pulled down subunit
        for ipd = studyBRGpdIndices
            if ipd == iko
                # The KOed subunit is inhibited
                add_pulldown_edge!(inhibitEdge, pGraphBRG, ipd)
            else
                # If the subunit is the last in the node list,
                # its index has been swapped with the deleted node
                if ipd == M
                    if !(iko in pulledBRGComponent)
                        add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                        competitionDict, ipd, iko)
                    end                  
                    continue # Look at next pulldowned subunit
                else
                    if !(ipd in pulledBRGComponent)
                        # If a subunit is not in the component connected
                        # to ARID1A, the KO will decrease the quantity of
                        # this subunit that will be pulled-down
                        add_pulldown_edge!(inhibitEdge, pGraphBRG, iko, ipd)
                    else
                        # The PD subunit is connected
                        enhanceIfDisconnectedCompetition!(pGraphBRG, pulledBRGComponent,
                            competitionDict, ipd, iko)
                    end
                end
            end
        end 
    end
    
    return(pGraphARID, pGraphBRG)
end

"""
Return a list of all subunits still connected
to ARID1A after a given KO is performed
"""        
function getARIDPulledComponent(graph::LightGraphs.SimpleGraphs.SimpleGraph{Int64}, iko::Int64)
    perturbGraph = copy(graph)
    rem_vertex!(perturbGraph, iko)
    pulledComponent = Array{Int64,1}
    for component in connected_components(perturbGraph) if aridIndex in component
        return(component)
    end end
end

"""
Return a list of all subunits still connected
to SMARCA4 after a given KO is performed
"""        
function getBRGPulledComponent(graph::LightGraphs.SimpleGraphs.SimpleGraph{Int64}, iko::Int64)
    perturbGraph = copy(graph)
    rem_vertex!(perturbGraph, iko)
    pulledComponent = Array{Int64,1}
    for component in connected_components(perturbGraph) if brgIndex in component
        return(component)
    end end
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
        # For all subunits not in the competition dict
        if !(i in keys(competition))
            # Continue until a competition class has been attributed
            while true
                # Assign random competition class
                newComp = rand(1:M)
                if all([competition[n] != newComp for n in intersect(neighbors(graph, i), keys(competition))])
                    competition[i] = rand(1:M)
                    break
                end
                # This competition would link interactors, try again
            end
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

"""
Cross-over between two structure graphs
"""
function crossOverGraphs!(sGraph1::structureGraph, sGraph2::structureGraph)
    return(1)
end

# Genetic algorithm module

"""
Compute loss for a given structure
compared to observation
"""
function observedLoss(sGraph::structureGraph,
    details::Bool = false, alpha::Float64 = ALPHA)
    pGraphARID, pGraphBRG = structureToPulldowns(sGraph)
    
    observedEdgesARID = edgeTypesARID[alpha]
    observedEdgesBRG = edgeTypesBRG[alpha]
    
    intersectEdgesARID = intersect(pGraphARID.edges, observedEdgesARID)
    intersectEdgesBRG = intersect(pGraphBRG.edges, observedEdgesBRG)
    unionEdgesARID = union(pGraphARID.edges, observedEdgesARID)
    unionEdgesBRG = union(pGraphBRG.edges, observedEdgesBRG)
    
    if details
        # Return array with Jaccard index
        # length of union and length of  
        return([(length(intersectEdgesARID)+length(intersectEdgesBRG)) / (length(unionEdgesARID)+length(unionEdgesBRG)),
                length(intersectEdgesARID), length(pGraphARID.edges), length(intersectEdgesBRG), length(pGraphBRG.edges)])
    else
        # Return Jaccard index
        return([(length(intersectEdgesARID)+length(intersectEdgesBRG)) / (length(unionEdgesARID)+length(unionEdgesBRG))])
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
    
    # Cross-over
#     if rand() < p_crs
#         sGraph1 = rand(pop)
#         sGraph2 = rand(pop)
#         if sGraph1 != sGraph2
#             crossOverGraphs!(sGraph1, sGraph2)
#         end
#     end
    
    return(fitness)
end

# Run genetic algorithm

# Max number of edges in a graph
maxEdges = Int64(M*(M-1)/2)

# Initialize population
pop = map(x -> structureGraph(
        Graph(M, rand(1:maxEdges)),
        copy(unitDictString),
        Dict(e => e for e in 1:M)),
    1:N)

# Ensure connectivity
map(connectGraph!, pop)

# Iterate and store performance
quantileFitness = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileIntersectARID = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileSimulatedEdgesARID = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileIntersectBRG = Array{Float16}(Int(ceil(L/monitorStep)), 5)
quantileSimulatedEdgesBRG = Array{Float16}(Int(ceil(L/monitorStep)), 5)
for i in 1:L
    if i % monitorStep == 1
        f = newGeneration!(pop, true, p_add = P)
        currentStep = Int(ceil(i/monitorStep))
        quantileFitness[currentStep,:] = quantile(map(x -> x[1], f))
        quantileIntersectARID[currentStep,:] = quantile(map(x -> x[2], f))
        quantileSimulatedEdgesARID[currentStep,:] = quantile(map(x -> x[3], f))
        quantileIntersectBRG[currentStep,:] = quantile(map(x -> x[4], f))
        quantileSimulatedEdgesBRG[currentStep,:] = quantile(map(x -> x[5], f))
    else
        f = newGeneration!(pop, false, p_add = P)
    end
end

save(outputFileName*".jld","pop", pop,
    "fitness", quantileFitness, "intersectARID", quantileIntersectARID, "quantileSimulatedEdgesARID", quantileSimulatedEdgesARID,
    "intersectBRG", quantileIntersectBRG, "quantileSimulatedEdgesBRG", quantileSimulatedEdgesBRG)
