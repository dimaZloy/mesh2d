

using PyPlot;
using WriteVTK;
using CPUTime;
using Distributed;
using BSON: @load
using BSON: @save

## TODO:
## calculateNode2cellsL2matrix - verify 
## computeCellStiffness2D - SPEEDUP

include("primeObjects.jl"); ## solver2d objects and controls 
include("readGambitNeuFile.jl");
include("utilsMesh2D.jl");
include("preprocessing.jl");
include("preprocessSimpleTriMesh.jl");
include("preprocessSimpleQuadMesh.jl");


preProcess("testMixedMesh2d.neu");
preProcessSimpleTriMesh();
preProcessSimpleQuadMesh();
preProcess("testMeshStep2d.neu");



