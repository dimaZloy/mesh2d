
using PyPlot;

using Pkg;
using Dates;
using Printf;
using DelimitedFiles;
using CPUTime;
using WriteVTK;
using Distributed;
using BSON: @load
using BSON: @save


# principal mesh varaibales:: 
# mesh_nodes[nNodesx3]
# mesh_connectivity [nCellsx3]
# bc_data
# bc_indexes
# cell_nodes_X [nCellsx4]
# cell_nodes_Y [nCellsx4]
# cell_mid_points [nCellsx2]
# cell_areas [nCellsx1]
# cell_edges_Nx [nCellsx4]
# cell_edges_Ny [nCellsx4]
# cell_edges_length [nCellsx4]
# cell_stiffness [nCellsx4]
# cell_clusters [nNodesx8]
# node_stencils [nNodesx8]
# cell2nodes [nCellsx8]
# AUX 
# node2cellsL2up  [nCellsx3]
# node2cellsL2down [nCellsx3]


include("primeObjects.jl");	
include("utilsMesh2D.jl");	
include("preprocessing.jl");	

include("simpleTriMesh.jl");	


println("nNodes=", nNodes);
println("nCells=", nCells);


xNodes = zeros(Float64, nNodes);
yNodes = zeros(Float64, nNodes);

xNodes = mesh_nodes[:,2];
yNodes = mesh_nodes[:,3];

#meshNodesX = view(mesh_nodes,:,2,:);
#meshNodesY = view(mesh_nodes,:,3,:);

display("calculate bounding box for the computational domain ...");

(minX,txt) = findmin(xNodes); 
(maxX,txt) = findmax(xNodes);

(minY,txt) = findmin(yNodes);
(maxY,txt) = findmax(yNodes);


println("minX=", minX);
println("maxX=", maxX);
println("minY=", minY);
println("maxY=", maxY);

# const boundingBox = [
	# minX minY
	# maxX minY
	# maxX maxY
	# minX maxY
	# minX minY
# ];

nNeibCells = 8; 

display("compute cells-related data...");
(cell_nodes_X, cell_nodes_Y, cell_mid_points) = reconstructionCells2Nodes2D(nCells,mesh_nodes,mesh_connectivity); #ok


xCells = zeros(Float64, nCells);
yCells = zeros(Float64, nCells);
xCells = cell_mid_points[:,1];
yCells = cell_mid_points[:,2];



display("compute cell areas...");
cell_areas = computeCellsAreas2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y); #ok


display("compute edge normals...");
(cell_edges_Nx, cell_edges_Ny, cell_edges_length) = computeCellNormals2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y); #ok




display("compute cells connectivity...");
CPUtic();
cell_stiffness = computeCellStiffness2D(nCells, bc_indexes, bc_data, mesh_connectivity); #ok 
CPUtoc();

display("compute cell clusters...");
cell_clusters = computeCellClusters2D(nNodes,nCells,nNeibCells, mesh_connectivity); #ok 

display("compute node stencils...");
node_stencils = computeNodeStencilsSIMPLEX2D(nNodes, nNeibCells, mesh_nodes,cell_clusters, cell_mid_points); 
#node_stencils = computeNodeStencilsSIMPLEX2D(nNodes, xNodes,yNodes,cell_clusters, cell_mid_points); 

# xN = zeros(4,nCells);
# yN = zeros(4,nCells);


# for i=1:nCells

	# z = mesh_connectivity[i,2];
       

	# if (z ==2)

		# xN[1,i] =		  cell_nodes_X[i,1];
		# xN[2,i] =		  cell_nodes_X[i,2];
		# xN[3,i] =		  cell_nodes_X[i,3];
		# xN[4,i] =		  cell_nodes_X[i,1];


		# yN[1,i] =		  cell_nodes_Y[i,1];
		# yN[2,i] =		  cell_nodes_Y[i,2];
		# yN[3,i] =		  cell_nodes_Y[i,3];
		# yN[4,i] =		  cell_nodes_Y[i,1];


	# elseif (z == 3)

		# xN[1,i] =		  cell_nodes_X[i,1];
		# xN[2,i] =		  cell_nodes_X[i,2];
		# xN[3,i] =		  cell_nodes_X[i,3];
		# xN[4,i] =		  cell_nodes_X[i,1];


		# yN[1,i] =		  cell_nodes_Y[i,1];
		# yN[2,i] =		  cell_nodes_Y[i,2];
		# yN[3,i] =		  cell_nodes_Y[i,3];
		# yN[4,i] =		  cell_nodes_Y[i,1];

		
	# end #if


# end #for

# display("compute node2cellsL2 matrix ... ");
# include("nodes2cellsL2.jl");

# node2cellsL2up, node2cellsL2down = calculateNode2cellsL2matrix();

Z = 1.0 ./cell_areas;



testMesh = mesh2d(
		nCells,
		nNodes,
		nNeibCells,
		nBSets,
		xNodes,
		yNodes,
		mesh_connectivity,
		bc_data,
		bc_indexes,
		cell_nodes_X,
		cell_nodes_Y,
		cell_mid_points,
		cell_areas,
		Z,
		cell_edges_Nx,
		cell_edges_Ny,
		cell_edges_length,
		cell_stiffness,
		cell_clusters,
		node_stencils,
		# cell2nodes,
		# AUX:
		# node2cellsL2up,
		# node2cellsL2down
	);

	fname = "testMesh02";
	
	fnameBSON = string(fname,".bson")
	fnameVTK =  string(fname)
	
	@save fnameBSON testMesh
	saveMeshToVTK(nCells, nNodes, xNodes, yNodes, mesh_connectivity, fnameVTK);

# @everywhere const cell_areasX  = $cell_areas;
# @everywhere const cell_edges_NxX = $cell_edges_Nx;
# @everywhere const cell_edges_NyX = $cell_edges_Ny;
# @everywhere const cell_edges_lengthX = $cell_edges_length;
# @everywhere const cell_stiffnessX = $cell_stiffness;
# @everywhere const cell_clustersX = $cell_clusters;
# @everywhere const node_stencilsX = $node_stencils;
# @everywhere const node2cellsL2upX = $node2cellsL2up;
# @everywhere const node2cellsL2downX = $node2cellsL2down; 
# @everywhere const ZX = $Z; 

