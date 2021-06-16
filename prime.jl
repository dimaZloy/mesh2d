

using PyPlot;
using WriteVTK;
using CPUTime;
using BSON: @load
using BSON: @save

## TODO:
## calculateNode2cellsL2matrix - verify 
## computeCellStiffness2D - SPEEDUP


include("readGambitNeuFile.jl");
include("utilsMesh2D.jl");
include("preprocessing.jl");




preProcess("testMesh00.neu");
preProcess("testMesh01.neu");



