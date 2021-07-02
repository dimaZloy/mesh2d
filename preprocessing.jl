
# using PyPlot;
# using WriteVTK;
# using CPUTime;
# using BSON: @load
# using BSON: @save



function preProcess(meshFile::String,nThreads::Int64)

	
	debugPlotMesh = false;

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


	nNeibCells::Int64 = 8; 

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
	
	display("compute cells connectivity serial ...");
	CPUtic();
	cell_stiffnessSerial = computeCellStiffnessM2D(nCells, bc_indexes, bc_data, mesh_connectivity); #ok 
	CPUtoc();
	
	display("compute cells connectivity distributed ... ");
	CPUtic();
	cell_stiffness = computeCellStiffnessDistributed(nCells, nThreads, bc_indexes,	bc_data, mesh_connectivity);
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

	(maxAreaI,id) = findmax(cell_areas);
	maxArea = sqrt(maxAreaI);
	
	(maxSide,id) = findmax(cell_edges_length);
	maxSideLength = maxSide;
	
	VTKCells = MeshCell[];
	
	for i=1:nCells
	
		cellType::Int64 = mesh_connectivity[i,2];
		
		if (cellType == 2) ## quads
			inds = Array{Int32}(undef, 4);
			inds[1] = mesh_connectivity[i,4];
			inds[2] = mesh_connectivity[i,5];
			inds[3] = mesh_connectivity[i,6];
			inds[4] = mesh_connectivity[i,7];
			c = MeshCell(VTKCellTypes.VTK_QUAD, inds);
			push!(VTKCells, c);
		
		elseif (cellType == 3) ## triangle
	
			inds = Array{Int32}(undef, 3);
			inds[1] = mesh_connectivity[i,4];
			inds[2] = mesh_connectivity[i,5];
			inds[3] = mesh_connectivity[i,6];
			c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds);
			push!(VTKCells, c);
		end
	end


	display("Compute computeNode2CellsL2 matrices ... ")

	node2cellL2up = zeros(Int64,nCells,8);
	node2cellL2down = zeros(Int64,nCells,8);
	
	computeNode2CellsL2(nCells,mesh_connectivity, cell_stiffness,node2cellL2up, node2cellL2down);
	
	display("done ... ")
	

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
		maxArea,
		maxSideLength,
		VTKCells,
		node2cellL2up,
		node2cellL2down
	);


	fname = split(meshFile, ".");
	
	fnameBSON = string(fname[1],".bson")
	fnameVTK =  string(fname[1])
	
	@save fnameBSON testMesh
	
	#saveMeshToVTK(nCells, nNodes, xNodes, yNodes, mesh_connectivity, fnameVTK);
	
	vtkfile = vtk_grid(fnameVTK, xNodes,yNodes, VTKCells);
	densityNodes = zeros(Float64,nNodes);
	vtk_point_data(vtkfile, densityNodes, "dummy");
	outfiles = vtk_save(vtkfile);	
	
	
	
end
