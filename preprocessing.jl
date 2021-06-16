
using PyPlot;
using WriteVTK;
using CPUTime;
using BSON: @load
using BSON: @save


struct mesh2d
	nCells::Int64
	nNodes::Int64
	nBSets::Int64
	xNodes::Array{Float64,1} 				##	mesh_nodes[nNodesx3]
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


function saveMeshToVTK(
	nCells::Int64,  ##number of cells
	nNodes::Int64,  ##number of nodes
	xNodes::Array{Float64,1}, ## x-coordinate of nodes
	yNodes::Array{Float64,1}, ## x-coordinate of nodes
	mesh_connectivity::Array{Int64,2},  ## mesh connectivity
	vtkFileName::String  ## VTK file name
	)


	cells = MeshCell[];
	#celltypeTri = VTKCellTypes.VTK_TRIANGLE;
	#celltypeQuad = VTKCellTypes.VTK_QUAD;
	
	for i=1:nCells
	
		cellType = mesh_connectivity[i,2];
		if (cellType == 2) ## quads
			inds = Array{Int32}(undef, 4);
			inds[1] = mesh_connectivity[i,4];
			inds[2] = mesh_connectivity[i,5];
			inds[3] = mesh_connectivity[i,6];
			inds[4] = mesh_connectivity[i,7];
			c = MeshCell(VTKCellTypes.VTK_QUAD, inds);
			push!(cells, c);
		
		elseif (cellType == 3) ## triangle
	
			inds = Array{Int32}(undef, 3);
			inds[1] = mesh_connectivity[i,4];
			inds[2] = mesh_connectivity[i,5];
			inds[3] = mesh_connectivity[i,6];
			c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds);
			push!(cells, c);
		end
	end

	vtkfile = vtk_grid(vtkFileName, xNodes,yNodes, cells);
	densityNodes = zeros(Float64,nNodes);
	vtk_point_data(vtkfile, densityNodes, "dummy");
	outfiles = vtk_save(vtkfile);	

end


function preProcess(meshFile::String)

	#meshFile = "testMesh00.neu";
	#meshFile = "testMesh01.neu";
	
	debugPlotMesh = false;

	##nCells, nNodes,  nBSets,  mesh_nodes,  mesh_connectivity,  BCNames,  BCindexes,  bc_data = readGambitNeuFile2(meshFile);
	nCells, nNodes,  nBSets,  mesh_nodes,  mesh_connectivity,  BCNames, bc_indexes,  bc_data = readGambitNeuFile2(meshFile);
	
	##bc_indexes  = convertBCdata(nBSets, BCindexes,bc_data );
	
	if (debugPlotMesh)
		plotMesh2d(mesh_nodes, mesh_connectivity,BCindexes,bc_data);
	end

	println("Computing mesh structure ... \n", );

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



	#display("compute node2cellsL2 matrix ... ");
	#include("nodes2cellsL2.jl");
	#node2cellsL2up, node2cellsL2down = calculateNode2cellsL2matrix(nCells, nNodes, mesh_connectivity, cell_stiffness); ## 

	Z = 1.0 ./cell_areas;


	testMesh = mesh2d(
		nCells,
		nNodes,
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

	fname = split(meshFile, ".");
	
	fnameBSON = string(fname[1],".bson")
	fnameVTK =  string(fname[1])
	
	@save fnameBSON testMesh
	saveMeshToVTK(nCells, nNodes, xNodes, yNodes, mesh_connectivity, fnameVTK);
	

end
