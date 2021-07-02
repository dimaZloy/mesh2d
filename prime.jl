

using Distributed;
const numThreads = 4;


if (numThreads != 1)

	if (nprocs() == 1)
		addprocs(numThreads,lazy=false); 
		display(workers());
	end
		
end

@everywhere using PyPlot;
@everywhere using WriteVTK;
@everywhere using CPUTime;
@everywhere using DelimitedFiles;
@everywhere using Printf
@everywhere using BSON: @load
@everywhere using BSON: @save
@everywhere using SharedArrays;


## TODO:
## calculateNode2cellsL2matrix - verify 
## computeCellStiffness2D - SPEEDUP


include("primeObjects.jl"); ## solver2d objects and controls 
include("computeCellsStiffnessMatrixDistributed.jl")
include("computeNode2CellsL2.jl")
include("readGambitNeuFile.jl");
include("utilsMesh2D.jl");
include("preprocessing.jl");
include("preprocessSimpleTriMesh.jl");
include("preprocessSimpleQuadMesh.jl");


preProcess("testMixedMesh2d.neu",numThreads);
preProcessSimpleTriMesh(numThreads);
preProcessSimpleQuadMesh(numThreads);

#preProcess("testMeshStep2d.neu");



