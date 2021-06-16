
# mesh_connectivity has the following format (nCells x 7):
# 1 element - global id,
# 2 element - cell type (2 for quads /3 for triangles)
# 3 element - number of nodes
# 4-7 - nodes ids
# 3 nodes for tri mesh 
# 4 nodes for quad mesh 

#bc_data (nBCesll x3)
# 1 element - global cell id,
# 2 element - cell type (2 for quads /3 for triangles)
# 3 element - global node index 

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
# AUX:
# node2cellsL2up  [nCellsx3]
# node2cellsL2down [nCellsx3]



function convertBCdata(nBSets, BCindexes, BCData)
	
	## BCData::Array{Int64,2}
	## BCIndexes::Array{Any,1}
	
	nBC::Int64 = size(BCData,1);
	#sumBC::Int64 = 0; 
	
	#for z=1:size(BCData,1)
	#	sumBC += BCindexes[z];
	#end

	#if ( (sumBC != size(BCData,1))  ||  (nBSets != size(BCindexes,1) ) )
	#	println("Warning! number of elements in boundary arrays is not correct ...") ;
	#end

	bc_indexes = ones(Int, nBC)*-nBSets;

	counter::Int64 = 0;	
	for k=1:nBSets-1
		
		for v= 1:BCindexes[k]
			counter +=1;	
			bc_indexes[counter] = -1*k	
		end
		
	end

	
	return bc_indexes;
end



function readGambitNeuFile2(fileName)

	if ( isfile(fileName) != true ) ### ~= true )
		println("Error to open file: ", fileName);
		return;
	end
		
	io = open(fileName,"r");
	
	
	
		println("reading gambit neu file: ", fileName);
	

		head1 = readline(io);
		head2 = readline(io);
		head3 = readline(io);
		head4 = readline(io);
		head5 = readline(io);
		head6 = readline(io);
		head7 = readline(io);
		head8 = readline(io);
	
	
		z = split(head7);
		MESH_DATA = zeros(Int64,6);
	
		for i=1:size(z,1)
			MESH_DATA[i] = parse(Int,z[i]);
		end
	
		println("NUMNP:\t",MESH_DATA[1] );
		println("NELEM:\t",MESH_DATA[2] );
		println("NGRPS:\t",MESH_DATA[3] );
		println("NBSETS:\t",MESH_DATA[4] );	
		println("NDFCD:\t",MESH_DATA[5] );
		println("NDFVL:\t",MESH_DATA[6] );
		
		nNodes = MESH_DATA[1];
		nCells = MESH_DATA[2];
		nBSets = MESH_DATA[4];
		
		
		line = readline(io);
	
		mesh_nodes = zeros(Float64,nNodes,3);
	
		for N=1:nNodes
			line = readline(io);
			z = split(line);
			for i=1:size(z,1)
				mesh_nodes[N,i] = parse(Float64,z[i]);
			end
		end
		line = readline(io);
	
		
		line = readline(io);
	
		mesh_connectivity = zeros(Int64,nCells,7);
		for N=1:nCells
			line = readline(io);
			z = split(line);
			for i=1:size(z,1)
				mesh_connectivity[N,i] = parse(Int64,z[i]);
			end
		end
		
		line = readline(io);
		
		# skip element groups
		#println( readline(io) );
		#println( readline(io) );
		#println( readline(io) );
	
		readline(io);
		readline(io);
		readline(io);
		
	
		while(true)
			line = readline(io);
			#println(line);
			if line == "ENDOFSECTION"
				break;
			end
		end
	
		######################################################################
	
	
		BCNames = [];
		BCData = zeros(Int,0,3);
		BCindexes = [];
	
		for B=1:nBSets

			head1 = readline(io);
			head2 = readline(io);
			#rintln(head1);
			#rintln(head2);
	
		
			z = split(head2);
			#isplay(z)
	
			BCNames = [BCNames; z[1] ];	
	
			id1 =parse(Int, z[2]);
			id2 =parse(Int, z[3]);
			id3 =parse(Int, z[4]);
			id4 =parse(Int, z[5]);	
	
			#println(z[1], '\t', id1,'\t', id2, '\t',  id3, '\t', id4);
	
			BCX = zeros(Int,id2,3);
			for C=1:id2
				line = readline(io);
				#display(line)
				z = split(line);
				#println(z);
				for i =1:size(z,1)
					BCX[C,i] = parse(Int,z[i]);	
				end
			end
			#BCData = [BCData, BCX];
			BCData = vcat( BCData, BCX );
			
			#display(BCX)
			BCindexes = [BCindexes; id2];
		    
			if (B != nBSets)
				readline(io);
			end
		
		end
	
		######################################################################
	
	close(io);	
	
	##nCells, nNodes,  nBSets,  mesh_nodes,  mesh_connectivity,  BCNames,  BCindexes,  bc_data = readGambitNeuFile2(meshFile);
	bc_indexes  = convertBCdata(nBSets, BCindexes, BCData);
	
	return nCells, nNodes, 	nBSets, mesh_nodes, mesh_connectivity, BCNames, bc_indexes, BCData;	
	

end



function plotMesh2d(mesh_nodes, mesh_connectivity,BCindexes,bc_data)

	debugPlotNodes = false;
	debugPlotInterior = true;
	debugPlotBoundaries = true;

	nbc = size(bc_data,1);
	xcb = [];
	ycb = [];

	for i =1:nbc
	
		idb = bc_data[i,1]; ## cell id 
		idv1 = bc_data[i,2]; ##  cell type
		idv2 = bc_data[i,3]; ## node id
	
		#p1 = mesh_connectivity[idb,3+idv1];
		p2 = mesh_connectivity[idb,3+idv2];
	
		xc2 =mesh_nodes[p2,2];
		push!(xcb, xc2);

		yc2 =mesh_nodes[p2,3];
		push!(ycb, yc2);
	
	end

	

	figure(1)
	clf();

	if (debugPlotNodes)
		plot(mesh_nodes[:,2],mesh_nodes[:,3],"or",markersize = 1.0);
	end
	
	if (debugPlotInterior)

		nCells = size(mesh_connectivity,1);
		for i=1:nCells

			if (mesh_connectivity[i,3] == 3)
				id1 = mesh_connectivity[i,4];
				id2 = mesh_connectivity[i,5];
				id3 = mesh_connectivity[i,6];
				
				xc = zeros(Float64,3);
				xc[1] =mesh_nodes[id1,2];
				xc[2] =mesh_nodes[id2,2];
				xc[3] =mesh_nodes[id3,2];

				yc = zeros(Float64,3);
				yc[1] =mesh_nodes[id1,3];
				yc[2] =mesh_nodes[id2,3];
				yc[3] =mesh_nodes[id3,3];
				
				plot( xc,yc,"-g",linewidth = 1.0);
			
			elseif (mesh_connectivity[i,3] == 4)
			
				id1 = mesh_connectivity[i,4];
				id2 = mesh_connectivity[i,5];
				id3 = mesh_connectivity[i,6];
				id4 = mesh_connectivity[i,7];

				xc = zeros(Float64,4);
				xc[1] =mesh_nodes[id1,2];
				xc[2] =mesh_nodes[id2,2];
				xc[3] =mesh_nodes[id3,2];
				xc[4] =mesh_nodes[id4,2];

				yc = zeros(Float64,4);
				yc[1] =mesh_nodes[id1,3];
				yc[2] =mesh_nodes[id2,3];
				yc[3] =mesh_nodes[id3,3];
				yc[4] =mesh_nodes[id3,3];
				
				plot( xc,yc,"-b",linewidth = 1.0);
				
			end
		end

	end	

	if (debugPlotBoundaries)
	
		bcid0 = 1;
		bcidI = BCindexes[1];
		plot( xcb[bcid0:bcidI],ycb[bcid0:bcidI],"sk",linewidth = 2.0,markersize = 1.0);

		for i = 2: size(BCindexes,1)

			bcid0 = bcidI+1;
			bcidI = bcidI + BCindexes[i];
			plot( xcb[bcid0:bcidI],ycb[bcid0:bcidI],"sk",linewidth = 2.0,markersize = 1.0);
		end
		
	end

	xlabel("x")
	ylabel("y")

end




