
#reconstruction_2D_cells_nodes(); 
#compute_2D_cells_stiffmatrix(); 
#compute_2D_cells_normals();
#compute_2D_cells_clusters(); 
#compute_2D_nodes_stensils_SIMPLEX();
#compute_2D_cells2nodes();
#compute2DMeshBoundingBox();
#distibuteCellsInThreads(); 

function computeQuadArea2D(x1::Float64, x2::Float64, x3::Float64, x4::Float64, y1::Float64, y2::Float64, y3::Float64, y4::Float64)::Float64
   return 0.5*( (x3-x1)*(y4-y2)-(x4-x2)*(y3-y1));
end

function computeTrgArea2D(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64, y3::Float64)::Float64
   return 0.5*((y3-y1)*(x2-x1)-(x3-x1)*(y2-y1));
end


function computeCellsAreas2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y)

	cell_areas = zeros(Float64, nCells);

	for i =1:nCells

		
		ct = mesh_connectivity[i,2];

                if (ct == 2)  
		   # quad
		   cell_areas[i] = computeQuadArea2D(
			cell_nodes_X[i,1],cell_nodes_X[i,2], cell_nodes_X[i,3],cell_nodes_X[i,4],
			cell_nodes_Y[i,1],cell_nodes_Y[i,2], cell_nodes_Y[i,3],cell_nodes_Y[i,4]);
		elseif (ct ==3)

		   # triangle
		    cell_areas[i] = computeTrgArea2D(
			cell_nodes_X[i,1],cell_nodes_X[i,2],cell_nodes_X[i,3],
			cell_nodes_Y[i,1],cell_nodes_Y[i,2],cell_nodes_Y[i,3]);

		end

	end # for 
	return cell_areas;
end


function computeNormal2Edge2D(x1::Float64, y1::Float64, x2::Float64, y2::Float64)


#   data 27.04.2010
#   Функция возвращает компоненты (внешней) нормали для прямой, заданной двумя точками,  а также длину прямой
#   Обход осуществляется против часовой сирелки...
#   Точка с координатами (mx,my) - середина данной прямой,
#   Точка с координатами (nx,ny) - вершина нормали. 

	N = zeros(Float64, 3);

	sp = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
	

	N[1] =   -(y2-y1)/sp;
	N[2] =    (x2-x1)/sp; 
	N[3] =   sp; 

	return N;
end







function computeCellNormals2D(nCells,mesh_connectivity,cell_nodes_X,cell_nodes_Y)

cell_edjes_Nx = zeros(Float64, nCells,4);
cell_edjes_Ny = zeros(Float64, nCells,4);
cell_edjes_lengths = zeros(Float64, nCells,4);

N1 = zeros(Float64, 3);
N2 = zeros(Float64, 3);
N3 = zeros(Float64, 3);
N4 = zeros(Float64, 3);


for i=1:nCells

	z = mesh_connectivity[i,2];

	if (z==2)

		N1 = computeNormal2Edge2D(cell_nodes_X[i,1], cell_nodes_Y[i,1],cell_nodes_X[i,2], cell_nodes_Y[i,2]);
		N2 = computeNormal2Edge2D(cell_nodes_X[i,2], cell_nodes_Y[i,2],cell_nodes_X[i,3], cell_nodes_Y[i,3]);
		N3 = computeNormal2Edge2D(cell_nodes_X[i,3], cell_nodes_Y[i,3],cell_nodes_X[i,4], cell_nodes_Y[i,4]);
		N4 = computeNormal2Edge2D(cell_nodes_X[i,4], cell_nodes_Y[i,4],cell_nodes_X[i,1], cell_nodes_Y[i,1]);
		
		cell_edjes_Nx[i,1] = N1[1];
		cell_edjes_Nx[i,2] = N2[1];
		cell_edjes_Nx[i,3] = N3[1];
		cell_edjes_Nx[i,4] = N4[1];

		cell_edjes_Ny[i,1] = N1[2];
		cell_edjes_Ny[i,2] = N2[2];
		cell_edjes_Ny[i,3] = N3[2];
		cell_edjes_Ny[i,4] = N4[2];

		cell_edjes_lengths[i,1] = N1[3];
		cell_edjes_lengths[i,2] = N2[3];
		cell_edjes_lengths[i,3] = N3[3];
		cell_edjes_lengths[i,4] = N4[3];


	elseif (z==3)

		N1 = computeNormal2Edge2D(cell_nodes_X[i,1], cell_nodes_Y[i,1],cell_nodes_X[i,2], cell_nodes_Y[i,2]);
		N2 = computeNormal2Edge2D(cell_nodes_X[i,2], cell_nodes_Y[i,2],cell_nodes_X[i,3], cell_nodes_Y[i,3]);
		N3 = computeNormal2Edge2D(cell_nodes_X[i,3], cell_nodes_Y[i,3],cell_nodes_X[i,1], cell_nodes_Y[i,1]);
		
		cell_edjes_Nx[i,1] = N1[1];
		cell_edjes_Nx[i,2] = N2[1];
		cell_edjes_Nx[i,3] = N3[1];

		cell_edjes_Ny[i,1] = N1[2];
		cell_edjes_Ny[i,2] = N2[2];
		cell_edjes_Ny[i,3] = N3[2];

		cell_edjes_lengths[i,1] = N1[3];
		cell_edjes_lengths[i,2] = N2[3];
		cell_edjes_lengths[i,3] = N3[3];

	end # if

end # for

   return cell_edjes_Nx, cell_edjes_Ny,cell_edjes_lengths; 

end


function reconstructionCells2Nodes2D(nCells::Int64,mesh_nodes,mesh_connectivity)

	#cout << endl << "cell's nodes reconstruction ... ";
	
	p1 = 0;
	p2 = 0;
	p3 = 0;
	p4 = 0;


	cell_nodes_X = zeros(Float64,nCells,4);
	cell_nodes_Y = zeros(Float64,nCells,4);
	cells_mid_points = zeros(Float64,nCells,2);


	for i=1:nCells
	

		z = mesh_connectivity[i,2]; 

		if (z == 2)

		#cell type - 2 for quads /3 for triangles


			p1 = mesh_connectivity[i,4];
			p2 = mesh_connectivity[i,5];
			p3 = mesh_connectivity[i,6];
			p4 = mesh_connectivity[i,7];


			cell_nodes_X[i,1] = mesh_nodes[p1,2];
			cell_nodes_X[i,2] = mesh_nodes[p2,2];
			cell_nodes_X[i,3] = mesh_nodes[p3,2];
			cell_nodes_X[i,4] = mesh_nodes[p4,2];

			cell_nodes_Y[i,1] = mesh_nodes[p1,3];
			cell_nodes_Y[i,2] = mesh_nodes[p2,3];
			cell_nodes_Y[i,3] = mesh_nodes[p3,3];
			cell_nodes_Y[i,4] = mesh_nodes[p4,3];
	

			cells_mid_points[i,1] = (cell_nodes_X[i,1]+cell_nodes_X[i,2]+cell_nodes_X[i,3]+cell_nodes_X[i,4])/4.0;
			cells_mid_points[i,2] = (cell_nodes_Y[i,1]+cell_nodes_Y[i,2]+cell_nodes_Y[i,3]+cell_nodes_Y[i,4])/4.0;
				

		elseif (z==3)
			
		
			p1 = mesh_connectivity[i,4];
			p2 = mesh_connectivity[i,5];
			p3 = mesh_connectivity[i,6];



			cell_nodes_X[i,1] = mesh_nodes[p1,2];
			cell_nodes_X[i,2] = mesh_nodes[p2,2];
			cell_nodes_X[i,3] = mesh_nodes[p3,2];

			cell_nodes_Y[i,1] = mesh_nodes[p1,3];
			cell_nodes_Y[i,2] = mesh_nodes[p2,3];
			cell_nodes_Y[i,3] = mesh_nodes[p3,3];


			cells_mid_points[i,1] = (cell_nodes_X[i,1]+cell_nodes_X[i,2]+cell_nodes_X[i,3])/3.0;
			cells_mid_points[i,2] = (cell_nodes_Y[i,1]+cell_nodes_Y[i,2]+cell_nodes_Y[i,3])/3.0;



		end #if 

		
	end #for 

	#cout << "   done" << endl;
	#cout << "cell's midpoints calculation ... done" << endl;

	return cell_nodes_X, cell_nodes_Y, cells_mid_points;

end




function computeCellClusters2D(nNodes::Int64, nCells::Int64, nNeibCells::Int64, mesh_connectivity)

#cluster_size = 7;
cell_clusters =zeros(Int64, nNodes, nNeibCells);



for i=1:nNodes

	c_node = i;
	counter = 0; 

	#cluster = zeros(1,cluster_size);

	cluster = zeros(nNeibCells);

	for j=1:nCells

		num_nodes::Int64 = mesh_connectivity[j,3];

		if (num_nodes == 3 )

			if (c_node == mesh_connectivity[j,4] || c_node == mesh_connectivity[j,5] || c_node == mesh_connectivity[j,6] )
				counter = counter +1;
				#cluster[1,counter] = j;
				cluster[counter] = j;
			end	

		elseif (num_nodes == 4)

			if (c_node == mesh_connectivity[j,4] || c_node == mesh_connectivity[j,5] || c_node == mesh_connectivity[j,6] || c_node == mesh_connectivity[j,7])
				counter = counter +1;
				#cluster[1,counter] = j;
				cluster[counter] = j;
			end	


		end

	end #end of j
	cell_clusters[i,:] = cluster;		

end

return cell_clusters;
end


function computeNodeStencilsSIMPLEX2D(nNodes::Int64, nNeibCells::Int64, mesh_nodes,cell_clusters, cell_mid_points)

#RX = [];
#RY = [];
#IXX = [];
#IYY = [];
#IXY = [];

rx = 0.0;
ry = 0.0;
ixx = 0.0;
iyy = 0.0;
ixy = 0.0;
x0 = 0.0;
y0 = 0.0;

#n_neib_cells = 7; 
#J = 1;

node_stencils = zeros(Float64, nNodes, nNeibCells);
RX = zeros(Float64, nNodes);
RY = zeros(Float64, nNodes);
IXX = zeros(Float64, nNodes);
IYY = zeros(Float64, nNodes);
IXY = zeros(Float64, nNodes);

for p =1:nNodes

	x0 = mesh_nodes[p,2];
	y0 = mesh_nodes[p,3];

	#x0 = mesh_nodes[J,1];
	#y0 = mesh_nodes[J,2];

	rx = 0.0;
	ry = 0.0;
	ixx = 0.0;
	iyy = 0.0;
	ixy = 0.0;
	
	
	for j =1:nNeibCells

		neib_cell = cell_clusters[p,j];

		if (neib_cell !=0)

			xI = cell_mid_points[neib_cell,1];
			yI = cell_mid_points[neib_cell,2];

			a = xI-x0;
			b = yI-y0;
			rx = rx +a;
			ry = ry +b;

			ixx = ixx + a*a;
			iyy = iyy + b*b;
			ixy = ixy + a*b;

		end #if

			

	end #j
	
	#RX = [RX; rx];
	#RY = [RY; ry];

	#IXX = [IXX; ixx];
	#IYY = [IYY; iyy];
	#IXY = [IXY; ixy];

	#J = J +1;

	RX[p] =  rx;
	RY[p] =  ry;

	IXX[p] = ixx;
	IYY[p] = iyy;
	IXY[p] = ixy;



end # p 


J = 1;
lymdaX = 0.0;
lymdaY = 0.0;



for p=1:nNodes

	x0 = mesh_nodes[p,2];
	y0 = mesh_nodes[p,3];

	#x0 = mesh_nodes[J,1];
	#y0 = mesh_nodes[J,2];


	stencils = zeros(nNeibCells);
	
	for j=1:nNeibCells

		neib_cell = cell_clusters[J,j];

		wi = 1.0;
		#wi = Float64(j);
		

		if (neib_cell != 0)

			xI = cell_mid_points[neib_cell,1];
			yI = cell_mid_points[neib_cell,2];
			a = xI-x0;
			b = yI-y0;			
			c = IXX[p] *IYY[p] - IXY[p] * IXY[p];

			if (c!=0.0)
				lymdaX = (IXY[p]*RY[p]-IYY[p]*RX[p])/c;
				lymdaY = (IXY[p]*RX[p]-IXX[p]*RY[p])/c;
			else
				lymdaX = 1.0;
				lymdaY = 1.0;
			end 
			
			dwi = lymdaX*a + lymdaY*b;
			wi = 1.0 + dwi;

			if (wi<=1e-2)
				wi =1.0;
			end

			
			

		end #if
		
		stencils[j] = wi;
		
	end #j

	#J = J+1;
	node_stencils[p,:] = stencils;	
end #p

return node_stencils;
end

function cells2nodesConnectivity2D(nCells, stiff_matrix, mesh_connectivity)
cells2nodes = zeros(Int64,nCells,7);

for C=1:nCells

	nodesC = mesh_connectivity[C,4:8];
	num_nodes=mesh_connectivity[C,3];
	
	if (num_nodes ==3)
		for T=1:3
			neib_cell = stiff_matrix[C,T];
			node1 = 0;
			node2 = 0;
			if (neib_cell >0) #internal cell
				nodesT = mesh_connectivity[neib_cell,4:8];
				neibCellType = mesh_connectivity[neib_cell, 2];
				NodesP =[]; 
				if (neibCellType ==3)
					if (nodesC[1] == nodesT[1] || nodesC[1] == nodesT[2] || nodesC[1] == nodesT[3])
						NodesP = [NodesP; nodesC[1]];
					end
					if (nodesC[2] == nodesT[1] || nodesC[2] == nodesT[2] || nodesC[2] == nodesT[3])
						NodesP = [NodesP; nodesC[2]];
					end
					if (nodesC[3] == nodesT[1] || nodesC[3] == nodesT[2] || nodesC[3] == nodesT[3])
						NodesP = [NodesP; nodesC[3]];
					end
				elseif (neibCellType ==2)
					if (nodesC[1] == nodesT[1] || nodesC[1] == nodesT[2] || nodesC[1] == nodesT[3] || nodesC[1] == nodesT[4])
						NodesP = [NodesP; nodesC[1]];
					end
					if (nodesC[2] == nodesT[1] || nodesC[2] == nodesT[2] || nodesC[2] == nodesT[3] || nodesC[2] == nodesT[4])
						NodesP = [NodesP; nodesC[2]];
					end
					if (nodesC[3] == nodesT[1] || nodesC[3] == nodesT[2] || nodesC[3] == nodesT[3] || nodesC[3] == nodesT[4])
						NodesP = [NodesP; nodesC[3]];
					end
					if (nodesC[4] == nodesT[1] || nodesC[4] == nodesT[2] || nodesC[4] == nodesT[3] || nodesC[4] == nodesT[4])
						NodesP = [NodesP; nodesC[4]];
					end
				

				end	

				if (size(NodesP,1) == 2)
					node1 = NodesP[1];
					node2 = NodesP[2];
				end
	
			else #boundary cell
				
				if (T ==1)
					node1 = NodesC[1];
					node2 = NodesC[2];
				elseif (T ==2)
					node1 = NodesC[2];
					node2 = NodesC[3];
				else 
					node1 = NodesC[3];
					node2 = NodesC[1];

				end

			end

			if (T==1)
				cells2nodes[C,1] = node1;
				cells2nodes[C,2] = node2;
			elseif (T==2)
				cells2nodes[C,3] = node1;
				cells2nodes[C,4] = node2;
			elseif (T==3)
				cells2nodes[C,5] = node1;
				cells2nodes[C,6] = node2;
			else
				println("somewthing wrong in cell2nodes fun ... ");
			end

		end

	elseif (num_nodes ==4)


		for T =1:4
			neib_cell = stiff_matrix[C,T];
			node1 =0;
			node2 =0;
			if (neib_cell >0)
				nodesT = mesh_connectivity[C,4:8];
				NodesP = [];

				if (nodesC[1] == nodesT[1] || nodesC[1] == nodesT[2] || nodesC[1] == nodesT[3] || nodesC[1] == nodesT[4])
					NodesP = [NodesP; nodesC[1]];
				end
				if (nodesC[2] == nodesT[1] || nodesC[2] == nodesT[2] || nodesC[2] == nodesT[3] || nodesC[2] == nodesT[4])
					NodesP = [NodesP; nodesC[2]];
				end
				if (nodesC[3] == nodesT[1] || nodesC[3] == nodesT[2] || nodesC[3] == nodesT[3] || nodesC[3] == nodesT[4])
					NodesP = [NodesP; nodesC[3]];
				end
				if (nodesC[4] == nodesT[1] || nodesC[4] == nodesT[2] || nodesC[4] == nodesT[3] || nodesC[4] == nodesT[4])
					NodesP = [NodesP; nodesC[4]];
				end

				if (size(NodesP,1) == 2)
					node1 = NodesP[1];
					node2 = NodesP[2];
				end
			else
				
				if (T ==1)
					node1 = NodesC[1];
					node2 = NodesC[2];
				elseif (T ==2)
					node1 = NodesC[2];
					node2 = NodesC[3];
				elseif (T ==3)
					node1 = NodesC[3];
					node2 = NodesC[4];
				else
					node1 = NodesC[4];
					node2 = NodesC[1];
				end

				if (T==1)
					cells2nodes[C,1] = node1;
					cells2nodes[C,2] = node2;
				elseif (T==2)
					cells2nodes[C,3] = node1;
					cells2nodes[C,4] = node2;
				elseif (T==3)
					cells2nodes[C,5] = node1;
					cells2nodes[C,6] = node2;
				elseif (T==4)
					cells2nodes[C,7] = node1;
					cells2nodes[C,8] = node2;
				else
					println("somewthing wrong in cell2nodes fun ... ");
				end

											
			end

		end




	end

end #C 

end #function 


# DEPRICATED!!!
# function computeCellStiffness2D(nCells,bc_indexes,bc_data,mesh_connectivity)
# cell_stiffness = zeros(Int64, nCells,4);

# #nbc = size(bc_indexes,1);
# nbc = length(bc_indexes);
# for m=1:nbc
    # i = bc_data[m,1];
    # k = bc_data[m,3];
    # cell_stiffness[i,k] = bc_indexes[m];				
# end

# p1 = 0;
# p2 = 0;
# p3 = 0;
# p4 = 0;

# k1 = 0;
# k2 = 0;
# k3 = 0;
# k4 = 0;

# #fCells = [];


# for i = 1:nCells

	# et::Int64 = mesh_connectivity[i,2];
	# if (et == 2)
		# p1 = mesh_connectivity[i,4];
		# p2 = mesh_connectivity[i,5];
		# p3 = mesh_connectivity[i,6];
		# p4 = mesh_connectivity[i,7];
	# elseif (et == 3)

		# p1 = mesh_connectivity[i,4];
		# p2 = mesh_connectivity[i,5];
		# p3 = mesh_connectivity[i,6];
	# end

	# for k=1:nCells

		# if (et == 2)
		
		# k1 = mesh_connectivity[k,4];
		# k2 = mesh_connectivity[k,5];
		# k3 = mesh_connectivity[k,6];
		# k4 = mesh_connectivity[k,7];

		# if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (i!=k) )
			# cell_stiffness[i,1] = k;
		# end
		# if ( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (i!=k) )
			# cell_stiffness[i,2] = k;
		# end
		# if ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (i!=k) )
			# cell_stiffness[i,3] = k;
		# end
		# if ( (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (i!=k) )
			# cell_stiffness[i,4] = k;
		# end

		# elseif (et == 3)

		# k1 = mesh_connectivity[k,4];
		# k2 = mesh_connectivity[k,5];
		# k3 = mesh_connectivity[k,6];

		# if ( (p1 == k1 || p1 == k2 || p1 == k3 ) && (p2 == k1 || p2 == k2 || p2 == k3 ) && (i!=k) )
			# cell_stiffness[i,1] = k;
		# end
		# if ( (p2 == k1 || p2 == k2 || p2 == k3 ) && (p3 == k1 || p3 == k2 || p3 == k3 ) && (i!=k) )
			# cell_stiffness[i,2] = k;
		# end
		# if ( (p3 == k1 || p3 == k2 || p3 == k3 ) && (p1 == k1 || p1 == k2 || p1 == k3 ) && (i!=k) )
			# cell_stiffness[i,3] = k;
		# end


		# end #if

	# end #k

# end #i


# return  cell_stiffness;
# end


function computeCellStiffnessM2D(nCells,bc_indexes,bc_data,mesh_connectivity)

	cell_stiffness = zeros(Int64, nCells,4);

	#nbc = size(bc_indexes,1);
	nbc = length(bc_indexes);
	for m=1:nbc
		i = bc_data[m,1];
		k = bc_data[m,3];
		cell_stiffness[i,k] = bc_indexes[m];				
	end

	p1::Int64 = 0;
	p2::Int64 = 0;
	p3::Int64 = 0;
	p4::Int64 = 0;

	k1::Int64 = 0;
	k2::Int64 = 0;
	k3::Int64 = 0;
	k4::Int64 = 0;

	fCells = [];


	for i = 1:nCells

		et::Int64 = mesh_connectivity[i,2];
		if (et == 2)
			p1 = mesh_connectivity[i,4];
			p2 = mesh_connectivity[i,5];
			p3 = mesh_connectivity[i,6];
			p4 = mesh_connectivity[i,7];
		elseif (et == 3)

			p1 = mesh_connectivity[i,4];
			p2 = mesh_connectivity[i,5];
			p3 = mesh_connectivity[i,6];
		end

		for k=1:nCells

			if (et == 2)
			
				k1 = mesh_connectivity[k,4];
				k2 = mesh_connectivity[k,5];
				k3 = mesh_connectivity[k,6];
				k4 = mesh_connectivity[k,7];

				if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (i!=k) )
					cell_stiffness[i,1] = k;
				end
				if ( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (i!=k) )
					cell_stiffness[i,2] = k;
				end
				if ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (i!=k) )
					cell_stiffness[i,3] = k;
				end
				if ( (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (i!=k) )
					cell_stiffness[i,4] = k;
				end

			elseif (et == 3)

				k1 = mesh_connectivity[k,4];
				k2 = mesh_connectivity[k,5];
				k3 = mesh_connectivity[k,6];

				if ( (p1 == k1 || p1 == k2 || p1 == k3 ) && (p2 == k1 || p2 == k2 || p2 == k3 ) && (i!=k) )
					cell_stiffness[i,1] = k;
				end
				if ( (p2 == k1 || p2 == k2 || p2 == k3 ) && (p3 == k1 || p3 == k2 || p3 == k3 ) && (i!=k) )
					cell_stiffness[i,2] = k;
				end
				if ( (p3 == k1 || p3 == k2 || p3 == k3 ) && (p1 == k1 || p1 == k2 || p1 == k3 ) && (i!=k) )
					cell_stiffness[i,3] = k;
				end


			end #if

		end #k
		
		 ## check for those cell who has zero index in connectivity        
		 if ( (et == 3) && (cell_stiffness[i,1] == 0 || cell_stiffness[i,2] == 0 || cell_stiffness[i,3] == 0)) 
			 ##display("Somethnig wrong in stiffMatrix2d fun ...");
			 #fCells = [fCells; i];
			 push!(fCells,i);
		 elseif ( (et == 2) && (cell_stiffness[i,1] == 0 || cell_stiffness[i,2] == 0 || cell_stiffness[i,3] == 0 || cell_stiffness[i,4] == 0)) 
			 ##disp("Somethnig wrong in stiffMatrix2d fun ...#);
			 #fCells = [fCells; i];
			 push!(fCells,i);
			 
		 end

		

	end #i


	
	n::Int64 = size(fCells,1);
 

	for z = 1:n
		
		i::Int64 = fCells[z];
		  
		element_type::Int64 = mesh_connectivity[i,2];
		
		p1 = mesh_connectivity[i,4];
		p2 = mesh_connectivity[i,5];
		p3 = mesh_connectivity[i,6];
		p4 = mesh_connectivity[i,7];
			
		zIndex::Int64 = -1000;
		
		if (cell_stiffness[i,1] == 0)
			zIndex = 1;
		elseif (cell_stiffness[i,2] == 0)
			zIndex = 2;
		elseif (cell_stiffness[i,3] == 0)
			zIndex = 3;
		else
			display("something wrong in cell stiffness calc ... ");
			break;
		end
				
		for k = 1:nCells
						
			cell_type::Int64 = mesh_connectivity[k,2];
			num_nodes::Int64 = mesh_connectivity[k,3];
										
			k1 = mesh_connectivity[k,4];
			k2 = mesh_connectivity[k,5];
			k3 = mesh_connectivity[k,6];
			k4 = mesh_connectivity[k,7];

						
						
			if ( cell_type == 2 && element_type == 3 )

			    ##  find first edge
                if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (i!=k) )
                    cell_stiffness[i,zIndex] = k;
                end
     			##  finding second edge
                if ( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (i!=k) )
                    cell_stiffness[i,zIndex] = k;
                end
 				##  finding third edge
                if ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (i!=k) )
                    cell_stiffness[i,zIndex] = k;
                end                   
						
			end

		end # end for k 
				
				
	 end #%end for i 



	return  cell_stiffness;
end



function distibuteCellsInThreads(nThreads::Int64, nCells::Int64 )

 
 	if (nThreads>1)
		
  		cellsThreads = zeros(Int64,nThreads,2);
		#cellsThreads = SharedArray{Int64}(nThreads,2);

    	#cout << "nThreads: " <<  nThreads << endl;
    	#cout << "nCells: " <<  get_num_cells() << endl;
    	nParts = floor(nCells/nThreads);
    	#cout << "nParts: " <<  nParts << endl;


	    for i=1:nThreads
    		cellsThreads[i,2] =  nCells - nParts*(nThreads-i );
		end
	

    	for i=1:nThreads
      		cellsThreads[i,1] =  cellsThreads[i,2] - nParts + 1;
		end

    	cellsThreads[1,1] = 1;

		#display(cellsThreads);	
		
		return cellsThreads;

    	#cout << "Partitioning mesh via threads ... " << endl;
    	#  for (int i=0;i<nThreads;i++)
    	#    cout << "cellsThreads[i][]: " << cellsThreads1[i] << '\t' << cellsThreads2[i] << endl;
    	#cout << "done" << endl;
		
	else				
		return 0;			
	end

end

