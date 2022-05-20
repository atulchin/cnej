include("sCGP.jl")
#import sCGP

tsteps = 1:2000;
in1 = Vector{Float64}();
out1 = Vector{Float64}();
x = 0.4;
for t in tsteps
	push!(in1,x);
	push!(out1,13.735*x);
	#3.5*(1.0-x))
	#x*(1.0-x))
	#2.0/(x*(1.0-x)))
	#3.5*x*(1.0-x))
	x = 3.5*x*(1.0-x);
end

flist = [sCGP.cnPass,sCGP.cnAdd,sCGP.cnDiff,sCGP.cnProd,sCGP.cnRatio,sCGP.cnSquare,sCGP.cnNeg,sCGP.cnOne,sCGP.cnRecip];
#flist = [sCGP.cnTanhSum,sCGP.cnLSum];


maxNodeSize = 100;
nInputs = 1;
nOutputs = 1;
#r = sCGP.initRun(length(dataSrc), length(targSrc), 500, 51:100, 1:0, 5; funcs=flist);
r = sCGP.initRun(maxNodeSize, nInputs, nOutputs, 100, 11:20, 1:0, 20; funcs=flist);
#r.muF = 0.25;
#r.muW = 0.5;
r.wtAdjRate = 0.2;
r.selStr = [0.1];
#r.selFunc = sCGP.selectOneKeepScore; 
#r.temperature = 10.0; r.cooling = 0.1;

gensPerTimeStep = 2000;
reportfreq = 1;
inputVectors = [in1[1:maxNodeSize]];
targetVectors = [out1[1:maxNodeSize]];
sCGP.changeEnv(r,inputVectors,targetVectors);

f = open("out.txt","w"); println(f,now());
for gen in 1:gensPerTimeStep
	if gen%reportfreq == 0
		printrep = true;
		linetxt = string("t ",1," g ",gen);
	else
		printrep = false;
	end
	sCGP.oneGen(r, printrep, f, linetxt);
end
println(f,now()); close(f);


sum([sCGP.isActive(n) for n in r.indivs[1].nodes])

x = 20.01
sCGP.testOutput(r,[[x]],[[13.735*x]],1)


maxNodeSize = 1;
nInputs = 1;
nOutputs = 1;
#r = sCGP.initRun(length(dataSrc), length(targSrc), 500, 51:100, 1:0, 5; funcs=flist);
r = sCGP.initRun(maxNodeSize, nInputs, nOutputs, 100, 11:20, 1:0, 20; funcs=flist);
#r.muF = 0.25;
#r.muW = 0.5;
r.wtAdjRate = 0.2;
r.selStr = [10.0];
#r.selFunc = sCGP.selectOneKeepScore; 
#r.temperature = 10.0; r.cooling = 0.1;

gensPerTimeStep = 1;
reportfreq = 1;
f = open("out.txt","w"); println(f,now());
for t in tsteps
	inputVectors = [in1[t:t]];
	targetVectors = [out1[t:t]];

	sCGP.changeEnv(r,inputVectors,targetVectors);

	for gen in 1:gensPerTimeStep
		if ((t-1)*gensPerTimeStep + gen)%reportfreq == 0
			printrep = true;
			linetxt = string("t ",t," g ",gen);
		else
			printrep = false;
		end
		sCGP.oneGen(r, printrep, f, linetxt);
	end
end
println(f,now()); close(f);



sum([sCGP.isActive(n) for n in r.indivs[1].nodes])

x = 0.1
sCGP.testOutput(r,[[x]],[[13.735*x]],1)







d1 = [0.4]
for i = 1:10000
	y::Float64 = d1[end]
	push!(d1,3.5*y*(1.0-y))
end
d1 = d1 .* 10.0;
d2 = [5.0]
for x in d1
	y = d2[end]
	push!(d2,x+y)
end
dataSrc = Vector{Float64}[[0.0;d1],d2];
targSrc = Vector{Float64}[d1,d2];



d1 = [0.4]
for i = 1:10000
	y = d1[end]
	push!(d1,3.81*y*(1.0-y))
end
d1 = d1 .- mean(d1);
d2 = [1.0]
for x in d1
	y = d2[end]
	push!(d2,16x^2-0.5y)
end
dataSrc = Vector{Float64}[[0.0;d1],d2];
targSrc = Vector{Float64}[d1,d2];





d2 = [0.4]
d3 = [0.16]
d4 = [0.4*3.78]
d5 = [0.24]
a1 = [0.0]
a2 = [0.0]
a3 = [0.0]
for i = 1:100000
	y = d2[end]
	push!(d2,3.78*y*(1.0-y))
	push!(d3,1.0-y)
	push!(d4,3.78*y)
	push!(d5,y*(1.0-y))
	push!(a1,1.5*y)
	push!(a2,2.0+y)
	push!(a3,2.5/y)
end
#d1 = ones(Float64,length(d2)) .* 3.78;
d1 = ones(Float64,length(d2));

#evolve to a3, then switch to d2
tvec = [a3[1:9000];d2[9001:end]];
targSrc = Vector{Float64}[tvec];
selTar = [100.0];

#fluctuate a1-a3, then switch to d2
tvec = [a1[1:1000];a2[1001:2000];a3[2001:3000];a1[3001:4000];a2[4001:5000];a3[5001:6000];a1[6001:7000];a2[7001:8000];a3[8001:9000];d2[9001:end]];
targSrc = Vector{Float64}[tvec];
selTar = [100.0];

#control
targSrc = Vector{Float64}[d2,a1,a2,a3];
selTar=[100.0,100.0,100.0,100.0]

#just evolve to d2
targSrc = Vector{Float64}[d2];
selTar = [100.0];

#evolve the pieces of d2
targSrc = Vector{Float64}[d2,d3,d4,d5];
selTar = [100.0,10.0,10.0,100.0];


dataSrc = Vector{Float64}[d1,d2];




#difficult to get from n parts to n+2 parts
d1 = [0.5]
d2 = [-1.25]
d2_part1 = [-1.25]
d2_part2 = [-1.25]
d2_part3 = [-1.25]
a1 = [0.0]
a2 = [0.0]
a3 = [0.0]

for i = 1:100000
	# x = a, y = b
	# d1 = [b]
	# d2 = [t1(x,y)]
	# x = b, y = t1
	# d1 = [b,t1]
	# d2 = [t1,t2]
	# x = t1, y = t2
	# d1 = [b, t1, t2...]
	# d2 = [t1, t2, t3(t1,t2)...]
	# data = [b, t1...], [t1, t2...]
	# targ = t3
	x = d1[end]
	y = d2[end]
	push!(d1,y)
	push!(d2, (-0.2*x + 2.75*y - y^3.0))
	#push!(d2_part1, (-0.2*x))
	#push!(d2_part2, (-1.0*y^3.0))
	#push!(d2_part3, (2.75*y - y^3.0))
	#push!(a1,1.8*x*x)
	#push!(a2,1.0+x)
	#push!(a3,3.3*x+y*y)
end

#evolve to a3, then switch to d2
#fluctuate a1-a3, then switch to d2

#control
targSrc = Vector{Float64}[d2[2:end]];
dataSrc = Vector{Float64}[d1,d2];


flist = [sCGP.cnPass,sCGP.cnAdd,sCGP.cnDiff,sCGP.cnProd,sCGP.cnRatio,sCGP.cnSquare,sCGP.cnNeg,sCGP.cnOne,sCGP.cnRecip,sCGP.cnTanhSum,sCGP.cnLSum];
#flist = [sCGP.cnTanhSum,sCGP.cnLSum];

r = sCGP.initRun(length(dataSrc),length(targSrc),250,5; funcs=flist);
r.wtAdjRate = 1.0;
r.selStr = [10.0];
sCGP.runGens(r, dataSrc, targSrc, 1:10000)


sum([n.active for n in r.indivs[1].nodes])

sCGP.testOutput(r,dataSrc,targSrc,1,85628)




