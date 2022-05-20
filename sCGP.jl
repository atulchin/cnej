
## 0.6
## -- node network operates on vectors; functions act element-wise on input vectors; each vector is considered a single input
##  	-- vectorize all functions, get fitness from dot(output-target , output-target)
##		-- use dot syntax in execnode; change Node.out type; preallocate Node.out


module sCGP

struct CNF
	f::Function
	arity::Int
end


### nodes ###
abstract type AbstractNode end


### unweighted node ###
mutable struct Node <: AbstractNode
	func::CNF
	out::Vector{Float64}
	active::Bool
	inCodes::Vector{Int} ##positions of connecting nodes
	inNodes::Vector{AbstractNode}
end

Node(maxSize::Int) = Node( CNF((a...)->nothing,0), fill(NaN,maxSize), false, Vector{Int}(), Vector{AbstractNode}()) #out vector is preallocated

Node(f::CNF, o::Vector{Float64}, a::Bool, ic::Vector{Int}) = Node(f, o, a, ic, Vector{AbstractNode}())

function Node(maxSize::Int, f::CNF)
	n = Node(maxSize)
	n.func = f
	return n
end

#copies node data, but connection pointers have to be established in new context
copyNode(p::Node) = Node(p.func, copy(p.out), p.active, copy(p.inCodes))

#accessors
val(n::Node) = n.out
isActive(n::Node) = n.active
inCodes(n::Node) = n.inCodes
inNodes(n::Node) = n.inNodes
arity(n::Node) = n.func.arity

function setFunc(n::Node, f::CNF)
	n.func = f
	return nothing
end

#copies values from src to dst, preserves size of dst by filling with fval, error if src larger than dst
function copyfill!(dst::AbstractArray{T,A}, src::AbstractArray{T,B}, fval::T) where {T,A,B}
	i = 0
	for i in eachindex(src)
		dst[i] = src[i]
	end
	if checkbounds(Bool,dst,i+1)
		dst[(i+1):end] = fval
	end
	return nothing
end

#n.out is pre-allocated, constant size
function setOutput(n::Node, o::Vector{Float64})
	copyfill!(n.out, o, NaN)
end

function setActive(n::Node, a::Bool)
	n.active = a
	return nothing
end

function setInCodes(n::Node, v::Vector{Int})
	n.inCodes = v
	return nothing
end

function setInCodes(n::Node, i::Int, c::Int)
	n.inCodes[i] = c
	return nothing
end

function setInputPointers(n::Node, v::Vector{AbstractNode})
	n.inNodes = v
	return nothing
end

#evaluate node function
function execNode(n::Node)
	n.out .= n.func.f.([val(inNode) for inNode in inNodes(n)]...) ## f is broadcast, n.out must be pre-allocated
	return n.out
end


### weighted node ###
mutable struct WeightedNode <: AbstractNode
	nodeData::Node
	inWeights::Vector{Float64}
end

WeightedNode(maxSize::Int) = WeightedNode( Node(maxSize), Vector{Float64}() )

WeightedNode(f::CNF, o::Vector{Float64}, a::Bool, ic::Vector{Int}, iw::Vector{Float64}) = WeightedNode( Node(f,o,a,ic), iw )

function WeightedNode(maxSize::Int, f::CNF)
	n = WeightedNode(maxSize)
	n.nodeData.func = f
	return n
end

#copies node data, but connection pointers have to be established in new context
copyNode(p::WeightedNode) = WeightedNode(p.nodeData.func, copy(p.nodeData.out), p.nodeData.active, copy(p.nodeData.inCodes), copy(p.inWeights))

#accessors
val(n::WeightedNode) = val(n.nodeData)
isActive(n::WeightedNode) = isActive(n.nodeData)
inCodes(n::WeightedNode) = inCodes(n.nodeData)
inNodes(n::WeightedNode) = inNodes(n.nodeData)
arity(n::WeightedNode) = arity(n.nodeData)
inWeights(n::WeightedNode) = n.inWeights

setFunc(n::WeightedNode, f::CNF) = setFunc(n.nodeData, f)
setOutput(n::WeightedNode, o::Vector{Float64}) = setOutput(n.nodeData, o)
setActive(n::WeightedNode, a::Bool) = setActive(n.nodeData, a)
setInCodes(n::WeightedNode, v::Vector{Int}) = setInCodes(n.nodeData, v)
setInCodes(n::WeightedNode, i::Int, c::Int) = setInCodes(n.nodeData, i, c)
setInputPointers(n::WeightedNode, v::Vector{AbstractNode}) = setInputPointers(n.nodeData, v)

function setWeights(n::WeightedNode, v::Vector{Float64})
	n.inWeights = v
	return nothing
end

function setWeights(n::WeightedNode, i::Int, w::Float64)
	n.inWeights[i] = w
	return nothing
end

#evaluate node function
function execNode(n::WeightedNode)
	n.nodeData.out .= n.nodeData.func.f.([inWeights(n)[i] * val(inNode) for (i,inNode) in enumerate(inNodes(n))]...) ## f broadcast, n pre-alloc
	return n.nodeData.out
end


### coefficient node ###
mutable struct CoeffNode <: AbstractNode
	out::Float64 ##keep as float; node functions taking this as input will broadcast it across a vector arg
	active::Bool
end

CoeffNode() = CoeffNode(1.0,false)
CoeffNode(f::CNF) = CoeffNode(1.0,false)
CoeffNode(f::CNF, o::Float64, a::Bool, ic::Vector{Int}) = CoeffNode(o,a)
copyNode(p::CoeffNode) = CoeffNode(p.out, p.active)
#accessors
val(n::CoeffNode) = n.out
isActive(n::CoeffNode) = n.active
inCodes(n::CoeffNode) = Vector{Int}()
inNodes(n::CoeffNode) = Vector{AbstractNode}()
arity(n::CoeffNode) = 0

function setFunc(n::CoeffNode, f::CNF)
	return nothing
end

function setOutput(n::CoeffNode, o::Float64)
	n.out = o
	return nothing
end

function setActive(n::CoeffNode, a::Bool)
	n.active = a
	return nothing
end

function setInCodes(n::CoeffNode, v::Vector{Int})
	return nothing
end

function setInCodes(n::CoeffNode, i::Int, c::Int)
	return nothing
end

function setInputPointers(n::CoeffNode, v::Vector{AbstractNode})
	return nothing
end

function execNode(n::CoeffNode)
	return n.out
end

###

mutable struct Indiv
	nodes::Vector{AbstractNode}
	outputNodes::Vector{AbstractNode}
	fitness::Float64
	isNew::Bool #set this to true when the phenotype changes (change in an active node)
	id::Int
	level::Int
	score::Float64

	Indiv(fit::Float64, n::Bool, id::Int) = new(Vector{AbstractNode}(), Vector{AbstractNode}(), fit, n, id, 1, 0.0)
end

Indiv(id::Int) = Indiv(0.0, true, id)

#build a new network instead of copying nodes
Indiv(p::Indiv, id::Int) = Indiv(p.fitness, p.isNew, id)


mutable struct RunParams
	maxNodeSize::Int #to avoid re-allocating node output each time a node is evaluated
	nInputs::Int
	nOutputs::Int
	netLen::Int #number of function nodes + number of input nodes
	idxFuncNodes::UnitRange{Int} #the range of network positions corresponding to function nodes
	idxCoeff::AbstractArray{Int} #which of the function nodes should be created as coefficient nodes
	idxWeighted::AbstractArray{Int} #which of the function nodes should be created as weighted nodes
	nBack::Int #defaults to netLen
	muOut::Float64 #p of mutation in output nodes (vs function nodes)
	muF::Float64 #p of mutation in node function (vs connections)
	muR::Float64 #p of connection mutations that are recurrent (feed-back)
	muW::Float64 #p of mutation in connection weight (vs connection location)
	wtAdjRate::Float64 #mutation effect size for node weights
	selStr::Vector{Float64}

	# sim annealing
	temperature::Float64
	cooling::Float64

	# which function to use for selection & reproduction
	selFunc::Function
	
	#all individuals point to a shared input layer
	inputNodes::Vector{Node}
	targetVectors::Vector{Vector{Float64}} #one vector for each network output
	indivs::Vector{Indiv} ##Dict{Indiv,Int}
	initPopSize::Int
	nodeFuncs::Vector{CNF}
	maxArity::Int
	indCount::Int
	storedIndivs::Vector{Indiv}

	function RunParams(maxSiz::Int, nIn::Int,nOut::Int,nFuncNodes::Int,idxC::AbstractArray{Int},idxW::AbstractArray{Int},n::Int)
		r = new()
		r.maxNodeSize = maxSiz
		r.nInputs = nIn
		r.nOutputs = nOut
		r.netLen = nFuncNodes + nIn
		r.idxFuncNodes = (nIn+1 : r.netLen)
		r.idxCoeff = idxC + nIn
		r.idxWeighted = idxW + nIn
		r.nBack = r.netLen
		r.muOut = nOut / (r.netLen + nOut)
		r.muF = 0.25
		r.muR = 0.0
		r.muW = 0.5
		r.wtAdjRate = 1.0 #at 1.0, mutating a weight sets it to a new random value
		r.selStr = ones(Float64, nOut)
		r.temperature = 100.0
		r.cooling = 1.0 / r.temperature #increase for faster cooling
		r.selFunc = selectOneElite
		r.inputNodes = [Node(maxSiz) for i=1:nIn]
		r.targetVectors = [fill(NaN,maxSiz) for i=1:nOut]
		r.indivs = [Indiv(0) for i=1:n] ##empty placeholder indivs
		r.initPopSize = n
		r.nodeFuncs = Vector{CNF}()
		r.maxArity = 8
		r.indCount = 0
		r.storedIndivs = copy(r.indivs)
		return r
	end
end

###each vector in rp.targetVectors is preallocated, constant size
function setTargetVectors(rp::RunParams, tv::Vector{Vector{Float64}})
	for (d::Vector{Float64},s::Vector{Float64}) in zip(rp.targetVectors, tv)
		copyfill!(d,s,NaN)
	end
end


####################################

function connectNodes(v::Vector{AbstractNode}) ##input vector should be the whole network
	for n::AbstractNode in v
		setInputPointers(n, v[inCodes(n)])
	end
	return nothing
end

#returns valid random connections for a given network position
function randConnections(loc::Int, ncon::Int, rp::RunParams)
	minCon = max(1, loc-rp.nBack)
	return rand(minCon:(loc-1), ncon)
end

#returns valid random connections for an output node
function randConnections(ncon::Int, rp::RunParams)
	minCon = max(1,rp.netLen-rp.nBack)
	return rand(minCon:rp.netLen, ncon)
end

#adjust a weight by a random amount based on rate paramter
function randWeight(currWt::Float64, adjRate::Float64)
	#x = rand() #U(0,1)
	x = tan(pi*(rand()-0.5)) #cauchy(0,1)
	return (x*adjRate + currWt*(1.0-adjRate))
end

function randomIndiv(rp::RunParams)
	rp.indCount += 1
	x = Indiv(rp.indCount)

	#first nInput nodes point to shared input nodes
	x.nodes = copy(rp.inputNodes) #this copies the pointers into a new array

	#make function nodes
	for i in rp.idxFuncNodes
		if i in rp.idxCoeff
			n = CoeffNode()
		elseif i in rp.idxWeighted
			n = WeightedNode(rp.maxNodeSize, rand(rp.nodeFuncs))
			setWeights(n, ones(Float64, rp.maxArity))
		else
			n = Node(rp.maxNodeSize, rand(rp.nodeFuncs))
		end
		setInCodes(n, randConnections(i, arity(n), rp))
		push!(x.nodes,n)
	end

	#output nodes
	for i = 1:rp.nOutputs
		n = Node(rp.maxNodeSize, cnPass)
		setInCodes(n, randConnections(arity(n), rp))
		push!(x.outputNodes,n)
	end

	#connect all nodes in network based on inCodes
	connectNodes([x.nodes;x.outputNodes])

	return x
end

function cloneIndiv(p::Indiv, rp::RunParams)
	rp.indCount += 1
	x = Indiv(p, rp.indCount)

	#first nInput nodes point to shared input nodes
	x.nodes = copy(rp.inputNodes) #this copies the pointers into a new array

	#copy function nodes (don't copy input nodes)
	for pnode in p.nodes[rp.idxFuncNodes]
		push!(x.nodes,copyNode(pnode))
	end
	#output nodes
	for pnode in p.outputNodes
		push!(x.outputNodes,copyNode(pnode))
	end
	#connect all nodes in network based on inCodes
	connectNodes([x.nodes;x.outputNodes])

	return x
end

function mutateNodeFunction(n::AbstractNode, loc::Int, rp::RunParams)
	origArity = arity(n)
	setFunc(n, rand(rp.nodeFuncs))
	##change connections if function arity changed
	if arity(n) < origArity
		#delete connections
		setInCodes(n, inCodes(n)[1:arity(n)])
	elseif arity(n) > origArity
		#make new random connections
		setInCodes(n, [inCodes(n) ; randConnections(loc, arity(n)-origArity, rp)])
	end
end

function mutateConnection(n::Node, loc::Int, rp::RunParams)
	i = rand(1:arity(n))
	setInCodes(n, i, randConnections(loc,1,rp)[1]) ##randConnections returns an array
end

function mutateConnection(n::WeightedNode, loc::Int, rp::RunParams)
	i = rand(1:arity(n))
	if rand() < rp.muW ##change a weight
		setWeights(n, i, randWeight(inWeights(n)[i],rp.wtAdjRate))
	else #change connection
		setInCodes(n, i, randConnections(loc,1,rp)[1]) ##randConnections returns an array
	end
end

#make a random change to a node; return true if node was active
function mutateNode(n::AbstractNode, loc::Int, rp::RunParams)
	if rand() < rp.muF ##change function
		mutateNodeFunction(n, loc, rp)
	else ##change connection or connection weight
		mutateConnection(n, loc, rp)
	end
	return isActive(n)
end

#for coeff nodes, any mutation just changes the output value
function mutateNode(n::CoeffNode, loc::Int, rp::RunParams)
	setOutput(n, randWeight(val(n), rp.wtAdjRate))
	return isActive(n)
end

#should flag active nodes before creating mutants
function makeMutations(ind::Indiv, rp::RunParams)
	activeNodeChanged = false ## introduce mutations until an active node is hit

	while (!activeNodeChanged)
		if rand() < rp.muOut ##change output node connection
			n = rand(ind.outputNodes)
			setInCodes(n, randConnections(arity(n), rp))
			activeNodeChanged = true
		else
			loc = rand(rp.idxFuncNodes)
			n = ind.nodes[loc]
			activeNodeChanged = mutateNode(n, loc, rp)
		end
	end

	#update connections
	connectNodes([ind.nodes;ind.outputNodes])

	ind.isNew = true #because an active node was changed
	return nothing
end

function updateInputs(rp::RunParams, inData::Vector{Vector{Float64}}) #one vector for each input node
	for (i,n) in enumerate(rp.inputNodes)
		setOutput(n, inData[i])
	end
	return nothing
end

#recursively flag active nodes
function checkNode(n::AbstractNode)
	setActive(n, true)
	for inputN in inNodes(n)
		checkNode(inputN)
	end
	return nothing
end

function checkActiveNodes(ind::Indiv)
	#clear active flags
	for n in ind.nodes
		setActive(n, false)
	end
	#network output goes into node.out for each node in outputLayer
	for n in ind.outputNodes
		checkNode(n)
	end
	return nothing

end

#evaluate all active nodes
function executeNodes(ind::Indiv, rp::RunParams)
	for n in [ind.nodes[rp.idxFuncNodes] ; ind.outputNodes] ##don't need to execute input nodes
		if isActive(n)
			execNode(n)
		end
	end
	return nothing
end

function calcFitness(g::Indiv, rp::RunParams)
	outputFits = zeros(Float64, rp.nOutputs) #fitness is 0 (worst) to 1 (best)
	for i = 1:rp.nOutputs
		s = rp.selStr[i]
		d::Vector{Float64} = val(g.outputNodes[i]) .- rp.targetVectors[i]
		filter!(a->!isnan(a),d)
		if length(d) > 0 #if this output produced NaNs for all input values, fitness for this output says at 0
			outputFits[i] = exp(-s * dot(d,d)) #fitness is a function of distance squared, scaled by selection strength
		end
	end
	g.fitness = mean(outputFits) #combine fitnesses from multiple outputs somehow
end

###
### selection and reproduction functions 
### -- makeNewGen() calls one of these after fitness calculation
###

function selectOneElite(rp::RunParams)
	sort!(rp.indivs, by=(x->x.fitness), rev=true)
	for i = 2:rp.initPopSize
		rp.indivs[i] = cloneIndiv(rp.indivs[1],rp)
		makeMutations(rp.indivs[i],rp)
	end
	return nothing
end

function selectOneEliteReplaceParent(rp::RunParams)
	rp.indivs[1].fitness -= eps() #subtract a tiny amount from the parent's fitness
	sort!(rp.indivs, by=(x->x.fitness), rev=true)
	for i = 2:rp.initPopSize
		rp.indivs[i] = cloneIndiv(rp.indivs[1],rp)
		makeMutations(rp.indivs[i],rp)
	end
	return nothing
end

function simAnneal(rp::RunParams)
	#rp.indivs[1] is the parent
	fit1 = rp.indivs[1].fitness; fit2 = rp.indivs[2].fitness
	if (fit2 >= fit1) || (rand()<exp((fit2-fit1)/rp.temperature)) 
		#parent is replaced
		rp.indivs[1] = rp.indivs[2]
	end
	#make a new offspring
	rp.indivs[2] = cloneIndiv(rp.indivs[1],rp)
	makeMutations(rp.indivs[2],rp)
	#update temp
	rp.temperature = 1.0 / (rp.cooling + 1.0/rp.temperature)
	return nothing
end

function selectOneEliteAnneal(rp::RunParams)
	#indivs[1] (parent) is the starting point
	#sort offspring highest to lowest fitness
	#higher fitness automatically takes the parent slot, but lower fitness can replace depending on temperature
	offspringByFitDesc = sort(rp.indivs[2:rp.initPopSize], by=(x->x.fitness), rev=true)
	for ind::Indiv in offspringByFitDesc
		if ind.fitness >= rp.indivs[1].fitness || rand() < exp((ind.fitness - rp.indivs[1].fitness) / rp.temperature)
			rp.indivs[1] = ind
		end
	end
	#update temp
	rp.temperature = 1.0 / (rp.cooling + 1.0/rp.temperature)
	#make offspring from rp.indivs[1]
	for i = 2:rp.initPopSize
		rp.indivs[i] = cloneIndiv(rp.indivs[1],rp)
		makeMutations(rp.indivs[i],rp)
	end
	return nothing
end

function selectOneKeepScore(rp::RunParams)
	#rp.indivs[1] contains the parent of the current generation
	#score = recent-weighted running avergage of own fitness + fitness of offspring
	rp.indivs[1].score = 0.0*rp.indivs[1].score + 1.0*sum(x.fitness for x in rp.indivs)  #replacing the score seems to work better?
	#score = recent-weighted running avergage of own fitness
	#rp.indivs[1].score = 0.0*rp.indivs[1].score + 1.0*rp.indivs[1].fitness  #replacing the score seems to work better?
	#replace stored indiv if score is higher
	if rp.indivs[1].score > rp.storedIndivs[1].score
		rp.storedIndivs[1] = rp.indivs[1]
	end
	#make next gen from highest fitness indiv
	sort!(rp.indivs, by=(x->x.fitness), rev=true)
	for i = 2:rp.initPopSize
		rp.indivs[i] = cloneIndiv(rp.indivs[1],rp)
		makeMutations(rp.indivs[i],rp)
	end
	#drop previous stored indiv
	rp.indivs = rp.indivs[1:rp.initPopSize]
	#include current stored indiv as contender for next gen
	push!(rp.indivs, rp.storedIndivs[1])
	return nothing
end

function selectOneKeepScoreAnneal(rp::RunParams)
	#rp.indivs[1] contains the parent of the current generation
	#score = recent-weighted running avergage of own fitness + fitness of offspring
	rp.indivs[1].score = 0.0*rp.indivs[1].score + 1.0*sum(x.fitness for x in rp.indivs)  #replacing the score seems to work better?
	#score = recent-weighted running avergage of own fitness
	#rp.indivs[1].score = 0.0*rp.indivs[1].score + 1.0*rp.indivs[1].fitness  #replacing the score seems to work better?
	#replace stored indiv if score is higher
	if rp.indivs[1].score > rp.storedIndivs[1].score
		rp.storedIndivs[1] = rp.indivs[1]
	end

	#the stored indivdiual is the starting point
	p = rp.storedIndivs[1]
	#the parent of the current gen, or one of the offspring, will become the parent of the next gen if its fitness is higher
	#lower fitness might become the parent instead, depending on temperature
	candidates = sort(rp.indivs[1:rp.initPopSize], by=(x->x.fitness), rev=true)
	for ind::Indiv in candidates
		if ind.fitness >= p.fitness || rand() < exp((ind.fitness - p.fitness) / rp.temperature)
			p = ind
		end
	end
	#update temp
	rp.temperature = 1.0 / (rp.cooling + 1.0/rp.temperature)

	#p becomes parent of next gen
	rp.indivs[1] = p
	#make next gen from indivs[1]
	for i = 2:rp.initPopSize
		rp.indivs[i] = cloneIndiv(rp.indivs[1],rp)
		makeMutations(rp.indivs[i],rp)
	end
	#drop previous stored indiv
	rp.indivs = rp.indivs[1:rp.initPopSize]
	#include current stored indiv as contender for next gen
	push!(rp.indivs, rp.storedIndivs[1])
	return nothing
end

function selectNplusN(rp::RunParams)
	sort!(rp.indivs, by=(x->x.fitness), rev=true)
	resize!(rp.indivs, rp.initPopSize) #shrink back down to n, dropping lowest fitness indivs
	for i = 1:rp.initPopSize
		newInd = cloneIndiv(rp.indivs[i],rp) #each surviving individual reproduces once
		makeMutations(newInd,rp)
		push!(rp.indivs, newInd)
	end
	#then change env, calc fitness, select n again
end

function selectionZiggurat(rp::RunParams)
	#anyone who didn't win goes down a level
	for ind in rp.indivs 
		ind.level -= 1
	end
	#winner goes up a level
	sort!(rp.indivs, by=(x->x.fitness), rev=true)
	rp.indivs[1].level += 2
	#everyone below level 1 is eliminated
	filter!((x::Indiv->x.level > 0), rp.indivs)
	#new individuals created at level 1
	for i = 2:rp.initPopSize
		newInd = cloneIndiv(rp.indivs[1],rp)
		makeMutations(newInd, rp)
		push!(rp.indivs, newInd)
	end
	return nothing
end

function makeNewGen(rp::RunParams)
	rp.selFunc(rp)
end

###

#population size n with random genomes
function populateWithRandoms(rp::RunParams)
	for i = 1:rp.initPopSize
		g = randomIndiv(rp)
		rp.indivs[i] = g
		rp.storedIndivs[i] = g
	end
	return nothing
end

#called after input or target data have changed
function processGenomes(rp::RunParams)
	for g::Indiv in rp.indivs
		checkActiveNodes(g)
		executeNodes(g, rp)
		calcFitness(g, rp)
	end
	return nothing
end

function printReport(io, linetxt::String, rp::RunParams)
	print(io,linetxt," fits = ")
	#for ind::Indiv in sort(rp.indivs, by=(x->x.fitness), rev=true)
	for ind::Indiv in rp.indivs
		print(io,ind.id::Int,": ",round(ind.fitness,5)," ")
	end
	print(io,"\n")
end

#gen1: make pop of randos
#		load gen1 inputs & outputs (changeEnv)
#		calc fitnesses (processGenomes)
#pick best (ties should go to non-parent?)
#make mutants (cloneIndiv + makeMutations)
#load gen2 inputs & outputs (changeEnv)
#calc fitnesses (processGenomes)
#etc
function oneGen(rp::RunParams, printrep::Bool=false, io=STDOUT, linetxt::String="")
	processGenomes(rp)
	if printrep
		printReport(io, linetxt, rp)
	end
	makeNewGen(rp)
end

function changeEnv(rp::RunParams, inData::Vector{Vector{Float64}}, targData::Vector{Vector{Float64}}) #one vector for each input/output
	setTargetVectors(rp, targData)
	updateInputs(rp, inData)
	return nothing
end

function initRun(maxNodeSize::Int, nIn::Int, nOut::Int, nFuncNodes::Int, idxC::AbstractArray{Int}, idxW::AbstractArray{Int}, psize::Int; funcs=[cnPass,cnAdd,cnDiff,cnProd,cnRatio,cnSquare,cnNeg,cnOne,cnRecip,cnTanhSum,cnLSum])
	rp = RunParams(maxNodeSize, nIn,nOut,nFuncNodes,idxC,idxW,psize)
	rp.nodeFuncs = funcs
	populateWithRandoms(rp)
	return rp
end

function testOutput(rp::RunParams, inputVectors::Vector{Vector{Float64}}, targetVectors::Vector{Vector{Float64}}, o::Int)
	inds = sort(rp.indivs, by=(x->x.fitness), rev=true)
	g = inds[1]
	
	changeEnv(rp,inputVectors,targetVectors)
	checkActiveNodes(g)
	executeNodes(g, rp)

	sz = length(targetVectors[o])
	return collect(zip(val(g.outputNodes[o])[1:sz], targetVectors[o]))
end

#node functions
isnum(x) = (!isnan(x))

cnPass = CNF(x->x, 1)
cnAdd = CNF(+, 2)
cnDiff = CNF(-, 2)
cnProd = CNF(*, 2)
function wrappedDivision(x::Number,y::Number)
	if iszero(y)
		return one(x)
	else
		return x/y
	end
end
cnRatio = CNF(wrappedDivision, 2)
cnRecip = CNF(x->wrappedDivision(1,x), 1)

cnSquare = CNF(x->x^2, 1)
cnCube = CNF(x->x^3, 1)

cnOnePlus = CNF(x->x+one(x), 1)
cnOneMinus = CNF(x->one(x)-x, 1)
cnNeg = CNF(x->(-x), 1)
cnOne = CNF(x->one(x), 1)

cnSin = CNF(sin, 1)
cnCos = CNF(cos, 1)

cnTanhSum = CNF((x,y) -> tanh(x+y), 2)
cnLSum = CNF((x,y) -> 1.0/(1.0+exp(-x-y)), 2)

#module
end







##
