
struct mesh2d
	nCells::Int64
	nNodes::Int64
	nNeibCells::Int64						## max number of neighbors 
	nBSets::Int64							##  number of boundaries  
	xNodes::Array{Float64,1} 				##  mesh_nodes[nNodesx3]
	yNodes::Array{Float64,1} 				##	mesh_nodes[nNodesx3]
	mesh_connectivity::Array{Float64,2} 	## [nCellsx3]
	bc_data::Array{Int64,2}
	bc_indexes::Array{Int64,1}
	cell_nodes_X::Array{Float64,2} 			## [nCellsx4]
	cell_nodes_Y::Array{Float64,2} 			## [nCellsx4]
	cell_mid_points::Array{Float64,2} 		## [nCellsx2]
	cell_areas::Array{Float64,1} 			## [nCellsx1]
	Z::Array{Float64,1} 					## [nCellsx1] 1/cell_areas
	cell_edges_Nx::Array{Float64,2} 		## [nCellsx4]
	cell_edges_Ny::Array{Float64,2} 		## [nCellsx4]
	cell_edges_length::Array{Float64,2} 	## [nCellsx4]
	cell_stiffness::Array{Float64,2} 		## [nCellsx4]
	cell_clusters::Array{Float64,2} 		## [nNodesx8]
	node_stencils::Array{Float64,2} 		## [nNodesx8]
	#cell2nodes::Array{Float64,2} 			## [nCellsx8]
	# AUX:
	#node2cellsL2up::Array{Float64,2} 		## [nCellsx3]
	#node2cellsL2down::Array{Float64,2} 		## [nCellsx3]
end

# flux approximation:
#  1 - AUSM+
#  2 - exact Riemann solver
#  3 - hybrid, AUSM + RB

# flux limiter:
#  1 - minmod 
#  2 - van Leer

# spatial approximation
#  1 - First order upwind;
#  2 - Second order upwind

#  spatial approximation
#  1 - 2nd order TVD RK 
#  2 - 3rd order TVD RK (depricated)
#  3 - 4th order general compact RK


@everywhere struct SOLVER2D
	FLUXtype::Int8;
	FLUXlimiter::Int8;	
	SpatialDiscretization::Int8;
	TimeDiscretization::Int8;
end

@everywhere struct CONTROLS
	CFL::Float64;
	dt::Float64; # time step
	timeStepMethod::Int8;  # fixed or #adaptive
	startTime::Float64; # actual physical time to start simulation
	stopTime::Float64; # actual physical time to stop simulation 
	plotResidual::Int8; # flag to plot residuals
	densityConstrained::Int8; # flag to constrain density
	minDensityConstrained::Float64;
	maxDensityConstrained::Float64
end

struct plotCONTROLS
	contoursType::Int8; # 0 - filled contours , 1 - contour lines only
	nContours::Int8; # number of contours to plot 
	rhoMINcont::Float64; #min density
	rhoMAXcont::Float64; #max density 
	productionVideo::Int8 # flag to production video 
end

struct MESH
	nCells::Int64; #number of cells
	nNodes::Int64; #number of nodes
	HX::Float64; # max cell size  in domain 
	maxXcoord::Float64; # bounding box: max x
	minXcoord::Float64; # bounding box: min x
	maxYcoord::Float64; # bounding box: max y
	minYcoord::Float64;	# bounding box: max y
end

struct outputCONTROLS
	verbosity::Int8;  
	header::String; 
	saveResiduals::Int8;
	saveResults::Int8; 
	fileNameResults::String;
	fileNameResiduals::String;
end

mutable struct DYNAMICCONTROLS
	flowTime::Float64; # actual physical time
	cpuTime::Float64; # cpu time
	tau::Float64; # tau 
	verIter::Int64; # iterator for verbosity (output)
	curIter::Int64; # global iterator for time steppings
	rhoMax::Float64; #max denisty in domain
	rhoMin::Float64; #min density in domain;
	velmax::Float64; #max velocity in domain
    globalPath::String; # loacl path to the code
	localTestPath::String; #path to the test 
	isSolutionConverged::Int8; # flag to show if the convergence criteria is satisfied
	isRunSimulation::Int8;  # flag to run or stop Simulation
end
