
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

	# display("calculate bounding box for the computational domain ...");

	# (minX,txt) = findmin(xNodes); 
	# (maxX,txt) = findmax(xNodes);

	# (minY,txt) = findmin(yNodes);
	# (maxY,txt) = findmax(yNodes);


	# println("minX=", minX);
	# println("maxX=", maxX);
	# println("minY=", minY);
	# println("maxY=", maxY);


	nNeibCells::Int64 = 8; 

	display("compute cells-related data...");
	CPUtic();
	(cell_nodes_X, cell_nodes_Y, cell_mid_points) = reconstructionCells2Nodes2D(nCells,mesh_nodes,mesh_connectivity); #ok
	CPUtoc();

	xCells = zeros(Float64, nCells);
	yCells = zeros(Float64, nCells);
	xCells = cell_mid_points[:,1];
	yCells = cell_mid_points[:,2];
	HX = zeros(Float64, nCells);


	display("compute cell areas...");
	CPUtic();
	cell_areas = computeCellsAreas2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y); #ok
	CPUtoc();


	display("compute edge normals...");
	CPUtic();
	(cell_edges_Nx, cell_edges_Ny, cell_edges_length, HX) = computeCellNormals2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y); #ok
	CPUtoc();
	
	# display("compute cells connectivity serial ...");
	# CPUtic();
	# cell_stiffnessSerial = computeCellStiffnessM2D(nCells, bc_indexes, bc_data, mesh_connectivity); #ok 
	# CPUtoc();
	
	display("compute cells connectivity distributed ... ");
	CPUtic();
	(cell_stiffness, cell_stiffnessSA, mesh_connectivitySA)  = computeCellStiffnessDistributed(nCells, nThreads, bc_indexes,	bc_data, mesh_connectivity);
	CPUtoc();	
	
	

	display("compute cell clusters distributed...");
	CPUtic();
	#cell_clusters  = computeCellClusters2D(nNodes,nCells,nNeibCells, mesh_connectivity); #ok 
	#display(cell_clusters)
	cell_clusters  = computeCellClustersDistributed(nNodes, nCells, nNeibCells, mesh_connectivitySA);
	CPUtoc();
	
	display("compute node stencils...");
	CPUtic();
	node_stencils = computeNodeStencilsSIMPLEX2D(nNodes, nNeibCells, mesh_nodes,cell_clusters, cell_mid_points); 
	CPUtoc();

	Z = 1.0 ./cell_areas;

	(maxAreaI,id) = findmax(cell_areas);
	maxArea = sqrt(maxAreaI);
	
	(maxSide,id) = findmax(cell_edges_length);
	maxSideLength = maxSide;
	


	display("Compute computeNode2CellsL2 matrices ... ")
	CPUtic();
	node2cellL2up = zeros(Int64,nCells,8);
	node2cellL2down = zeros(Int64,nCells,8);
	
	computeNode2CellsL2(nCells,mesh_connectivity, cell_stiffness,node2cellL2up, node2cellL2down);
	CPUtoc();
	display("done ... ")
	
	display("Compute compute cellsnodes matrix ... ")
	CPUtic();
	cells2nodes = zeros(Int64,nCells,8);
	computeCells2Nodes2D(nCells,mesh_connectivity, cell_stiffness, cells2nodes )
	CPUtoc();
	display("done ... ")
	
	
	###############################################################################################
	display("save mesh to VTK...  ")
	CPUtic();
	
	
	fname = split(meshFile, ".");	
	fnameBSON = string(fname[1],".bson")
	fnameVTK =  string(fname[1])

	
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
	
	#saveMeshToVTK(nCells, nNodes, xNodes, yNodes, mesh_connectivity, fnameVTK);
	
	vtkfile = vtk_grid(fnameVTK, xNodes,yNodes, VTKCells);
	densityNodes = zeros(Float64,nNodes);
	vtk_point_data(vtkfile, densityNodes, "dummy");
	outfiles = vtk_save(vtkfile);	
	
	CPUtoc();
	display("done")
	

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
		HX,
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
		node2cellL2down,
		cells2nodes
	);

	
	
	display("save mesh structure to *bson ")
	CPUtic();
	@save fnameBSON testMesh
	CPUtoc();
	display("done")
	
	
	
	
	##saveDistributedMesh2d(testMesh, fnameBSON);
	
end


# function saveDistributedMesh2d(testMesh::mesh2d, fnameBSON::String)


	
	# mesh_connectivity = SharedArray{Int64}(testMesh.nCells, 7); 
	# Z = SharedVector{Float64}(testMesh.nCells);	
	# cell_edges_Nx = SharedArray{Float64}(testMesh.nCells,4);
	# cell_edges_Ny = SharedArray{Float64}(testMesh.nCells,4);
	# cell_edges_length = SharedArray{Float64}(testMesh.nCells,4);
	# cell_stiffness = SharedArray{Int64}(testMesh.nCells,4);
	# node2cellsL2up = SharedArray{Int64}(testMesh.nCells,8);
	# node2cellsL2down = SharedArray{Int64}(testMesh.nCells,8); 


	# for i = 1:testMesh.nCells
		
		# mesh_connectivity[i,:] = testMesh.mesh_connectivity[i,:]
		# Z[i] = testMesh.Z[i]
		# cell_edges_Nx[i,:] = testMesh.cell_edges_Nx[i,:]
		# cell_edges_Ny[i,:] = testMesh.cell_edges_Ny[i,:]
		# cell_edges_length[i,:] = testMesh.cell_edges_length[i,:]
		# cell_stiffness[i,:] = testMesh.cell_stiffness[i,:]
		# node2cellsL2up[i,:] = testMesh.node2cellsL2up[i,:]
		# node2cellsL2down[i,:] = testMesh. node2cellsL2down[i,:];
		
	# end
	
	
	# testMesh_shared = mesh2d_shared(
		# mesh_connectivity,
		# Z,
		# cell_edges_Nx,
		# cell_edges_Ny,
		# cell_edges_length,
		# cell_stiffness,
		# node2cellsL2up,
		# node2cellsL2down
	# );

	# fnameBSON_shared = string(fnameBSON,"_shared");

	# @save fnameBSON_shared testMesh_shared

# end