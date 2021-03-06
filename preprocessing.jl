
# using PyPlot;
# using WriteVTK;
# using CPUTime;
# using BSON: @load
# using BSON: @save



function preProcess(meshFile::String,nThreads::Int32; scale::Float64 = 1.0)

	
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

	for i=1:nNodes
		xNodes[i] = mesh_nodes[i,2]*scale;
		yNodes[i] = mesh_nodes[i,3]*scale;
	end

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


	nNeibCells::Int32 = 8; 

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
	

	display("compute cell wall distances ...");
	
	CPUtic();
	cell_wall_distances = zeros(Float64,nCells,4);
	calcCellWallDistances(nCells, cell_nodes_X, cell_nodes_Y, mesh_connectivity, cell_mid_points, cell_wall_distances);
	
	#display(cell_wall_distances)
	
	CPUtoc();
	
	
	# display("compute cells connectivity serial ...");
	# CPUtic();
	# cell_stiffnessSerial = computeCellStiffnessM2D(nCells, bc_indexes, bc_data, mesh_connectivity); #ok 
	# CPUtoc();
	
	display("compute cells connectivity distributed ... ");
	CPUtic();
	(cell_stiffness, mesh_connectivitySA)  = computeCellStiffnessDistributed(nCells, nThreads, bc_indexes, bc_data, mesh_connectivity);
	CPUtoc();	
	
	
	#display(mesh_connectivitySA)

	display("compute cell clusters distributed...");
	CPUtic();
	#cell_clusters  = computeCellClusters2D(nNodes,nCells,nNeibCells, mesh_connectivity); #ok 
	#display(cell_clusters)
	cell_clusters  = computeCellClustersDistributed(nNodes, nCells, nNeibCells, nThreads, mesh_connectivitySA);
	CPUtoc();
	
	#display(cell_clusters)
	
	display("compute node stencils...");
	CPUtic();
	node_stencils = computeNodeStencilsSIMPLEX2D(nNodes, nNeibCells, mesh_nodes,cell_clusters, cell_mid_points); 
	CPUtoc();
	
	#display(node_stencils)

	# Z = 1.0 ./cell_areas;

	# (maxAreaI,id) = findmax(cell_areas);
	# maxArea = sqrt(maxAreaI);
	
	# (maxSide,id) = findmax(cell_edges_length);
	# maxSideLength = maxSide;
	


	display("Compute computeNode2CellsL2 matrices ... ")
	CPUtic();
	node2cellL2up = zeros(Int32,nCells,8);
	node2cellL2down = zeros(Int32,nCells,8);
	
	computeNode2CellsL2(nCells,mesh_connectivity, cell_stiffness,node2cellL2up, node2cellL2down);
	CPUtoc();
	display("done ... ")
	
	display("Compute compute cellsnodes matrix ... ")
	CPUtic();
	cells2nodes = zeros(Int32,nCells,8);
	computeCells2Nodes2D(nCells,mesh_connectivity, cell_stiffness, cells2nodes )
	#computeCells2Nodes2Dhybrid(nCells,mesh_connectivity, cell_stiffness, cells2nodes )
	CPUtoc();
	display("done ... ")
	
	
	###############################################################################################
	# display("save mesh to VTK...  ")
	# CPUtic();
	
	
	# fname = split(meshFile, ".");	
	# fnameBSON = string(fname[1],".bson")
	# fnameVTK =  string(fname[1])

	
	# VTKCells = MeshCell[];
	
	# for i=1:nCells
	
		# cellType::Int64 = mesh_connectivity[i,2];
		
		# if (cellType == 2) ## quads
			# inds = Array{Int32}(undef, 4);
			# inds[1] = mesh_connectivity[i,4];
			# inds[2] = mesh_connectivity[i,5];
			# inds[3] = mesh_connectivity[i,6];
			# inds[4] = mesh_connectivity[i,7];
			# c = MeshCell(VTKCellTypes.VTK_QUAD, inds);
			# push!(VTKCells, c);
		
		# elseif (cellType == 3) ## triangle
	
			# inds = Array{Int32}(undef, 3);
			# inds[1] = mesh_connectivity[i,4];
			# inds[2] = mesh_connectivity[i,5];
			# inds[3] = mesh_connectivity[i,6];
			# c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds);
			# push!(VTKCells, c);
		# end
	# end
	
	# #saveMeshToVTK(nCells, nNodes, xNodes, yNodes, mesh_connectivity, fnameVTK);
	
	# vtkfile = vtk_grid(fnameVTK, xNodes,yNodes, VTKCells);
	# densityNodes = zeros(Float64,nNodes);
	# vtk_point_data(vtkfile, densityNodes, "dummy");
	# outfiles = vtk_save(vtkfile);	
	
	# CPUtoc();
	# display("done")

	
	# display("save mesh structure to *bson ")
	# CPUtic();
	# @save fnameBSON VTKCells
	# CPUtoc();
	# display("done")
	
	
	display("save mesh structure to *hdf5 ")
	CPUtic();
	
	fname = split(meshFile, ".");		
	fn = string(fname[1],".h5");
	
	h5open(fn,"w") do file
	
		write(file,"nCells", nCells);
		write(file,"nNodes", nNodes);
		write(file,"nNeibCells", nNeibCells);
		write(file,"nBSets",nBSets );
		write(file,"xNodes",xNodes );
		write(file,"yNodes",yNodes);
		write(file,"mesh_connectivity",mesh_connectivity );
		write(file,"bc_data", bc_data);
		write(file,"bc_indexes",bc_indexes );
		write(file,"cell_nodes_X", cell_nodes_X);
		write(file,"cell_nodes_Y", cell_nodes_Y);
		write(file,"cell_mid_points", cell_mid_points);
		write(file,"cell_areas", cell_areas);
		write(file,"cell_wall_distances", cell_wall_distances);
		write(file,"HX", HX);
		write(file,"cell_edges_Nx",cell_edges_Nx );
		write(file,"cell_edges_Ny", cell_edges_Ny);
		write(file,"cell_edges_length",cell_edges_length );
		write(file,"cell_stiffness",cell_stiffness );
		write(file,"cell_clusters", cell_clusters);
		write(file,"node_stencils", node_stencils);
		#write(file,"VTKCells", VTKCells);
		write(file,"node2cellL2up",node2cellL2up );
		write(file,"node2cellL2down", node2cellL2down);
		write(file,"cells2nodes",cells2nodes );
	
	end
	
	CPUtoc();
	display("done")
	
	
	
	##saveDistributedMesh2d(testMesh, fnameBSON);
	
end

