## [GZM] Inducible Transcription Factors 
#	Mariana Gómez-Schiavon
#	August, 2019
#		Julia v.1.1.1
#		Required libraries:
#			DifferentialEquations
#			ParameterizedFunctions
#			Statistics
#			Distributions
#			DelimitedFiles

## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, an : Chose analysis type);
include(string("ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters

# Load functions & parameters:
using DelimitedFiles
using Distributions
mm = include(string("Md_",iARG.mm,".jl"));
fn = include(string("FN_SSs.jl"));
pO = copy(p);

# Run analysis
if(iARG.an=="ExSSs")
	p = copy(pO);
	open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,".txt"), "w") do io
		writedlm(io, [vcat("HT","XT0",[string(i) for i in mm.myODE.syms])],'\t');
		HT = 10 .^ collect(-3:0.1:3);
        for i in 1:length(HT)
			XT0 = 0.01*p[:iM]/p[:gY];
            ssR = fn.SS(mm.myODE, p, [XT0 0 HT[i]], 1e-6);
            writedlm(io, [hcat(HT[i],XT0,ssR)],'\t');
			XT0 = p[:iM]/p[:gY];
            ssR = fn.SS(mm.myODE, p, [XT0 0 HT[i]], 1e-6);
            writedlm(io, [hcat(HT[i],XT0,ssR)],'\t');
        end
	end
else
	println("ERROR: Undetermined analysis. Options: ExSSs")
end