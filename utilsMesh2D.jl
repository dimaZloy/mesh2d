
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


function computeCellsAreas2D(nCells::Int32,mesh_connectivity::Array{Int32,2},cell_nodes_X::Array{Float64,2},cell_nodes_Y::Array{Float64,2})::Array{Float64,1}

	cell_areas = zeros(Float64, nCells);

	for i =1:nCells

		
		ct::Int32 = mesh_connectivity[i,2];

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


function computeNormal2Edge2D(x1::Float64, y1::Float64, x2::Float64, y2::Float64)::Array{Float64,1}


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







function computeCellNormals2D(nCells::Int32,mesh_connectivity::Array{Int32,2},cell_nodes_X::Array{Float64,2},cell_nodes_Y::Array{Float64,2})

	cell_edjes_Nx = zeros(Float64, nCells,4);
	cell_edjes_Ny = zeros(Float64, nCells,4);
	cell_edjes_lengths = zeros(Float64, nCells,4);
	HX = zeros(Float64, nCells);

	N1 = zeros(Float64, 3);
	N2 = zeros(Float64, 3);
	N3 = zeros(Float64, 3);
	N4 = zeros(Float64, 3);


	for i=1:nCells

		z::Int32 = mesh_connectivity[i,2];

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
		
		
		
		max_edge_length, id = findmax( cell_edjes_lengths[i,:]);
		HX[i] = max_edge_length;
		
		

	end # for



   return cell_edjes_Nx, cell_edjes_Ny,cell_edjes_lengths, HX;

end


function calcPoint2LineDistance(x1::Float64,y1::Float64,x2::Float64,y2::Float64,px::Float64,py::Float64)::Float64

	m::Float64 = 0.0;
	if abs(x2-x1) > eps(Float64)
		m = (y2-y1)/(x2-x1);
	end
	
	A::Float64 = m;
	B::Float64 = -1.0;
	C = -y1 + m*x1;
	
	return abs(A*px + B*py+C)/sqrt(A*A + B*B);
	
end


function calcCellWallDistances(nCells::Int32, cell_nodes_X::Array{Float64,2}, cell_nodes_Y::Array{Float64,2}, 
		mesh_connectivity::Array{Int32,2}, cell_mid_points::Array{Float64,2}, cell_wall_distances::Array{Float64,2})

	## compute abs length of perpendecular from the cell mid point to the edges 

	for i=1:nCells

		z::Int32 = mesh_connectivity[i,2];

		if (z==2)

			cell_wall_distances[i,1] = calcPoint2LineDistance(cell_nodes_X[i,1], cell_nodes_Y[i,1],cell_nodes_X[i,2], cell_nodes_Y[i,2], cell_mid_points[i,1], cell_mid_points[i,2]);
			cell_wall_distances[i,2] = calcPoint2LineDistance(cell_nodes_X[i,2], cell_nodes_Y[i,2],cell_nodes_X[i,3], cell_nodes_Y[i,3], cell_mid_points[i,1], cell_mid_points[i,2]);
			cell_wall_distances[i,3] = calcPoint2LineDistance(cell_nodes_X[i,3], cell_nodes_Y[i,3],cell_nodes_X[i,4], cell_nodes_Y[i,4], cell_mid_points[i,1], cell_mid_points[i,2]);
			cell_wall_distances[i,4] = calcPoint2LineDistance(cell_nodes_X[i,4], cell_nodes_Y[i,4],cell_nodes_X[i,1], cell_nodes_Y[i,1], cell_mid_points[i,1], cell_mid_points[i,2]);
			


		elseif (z==3)
		
			cell_wall_distances[i,1] = calcPoint2LineDistance(cell_nodes_X[i,1], cell_nodes_Y[i,1],cell_nodes_X[i,2], cell_nodes_Y[i,2], cell_mid_points[i,1], cell_mid_points[i,2]);
			cell_wall_distances[i,2] = calcPoint2LineDistance(cell_nodes_X[i,2], cell_nodes_Y[i,2],cell_nodes_X[i,3], cell_nodes_Y[i,3], cell_mid_points[i,1], cell_mid_points[i,2]);
			cell_wall_distances[i,3] = calcPoint2LineDistance(cell_nodes_X[i,3], cell_nodes_Y[i,3],cell_nodes_X[i,4], cell_nodes_Y[i,4], cell_mid_points[i,1], cell_mid_points[i,2]);


		end # if
		
	end # for


end

function reconstructionCells2Nodes2D(nCells::Int32,mesh_nodes::Array{Float64,2},mesh_connectivity::Array{Int32,2})

	#cout << endl << "cell's nodes reconstruction ... ";
	
	p1::Int32 = 0;
	p2::Int32 = 0;
	p3::Int32 = 0;
	p4::Int32 = 0;


	cell_nodes_X = zeros(Float64,nCells,4);
	cell_nodes_Y = zeros(Float64,nCells,4);
	cells_mid_points = zeros(Float64,nCells,2);


	for i=1:nCells
	

		z::Int32 = mesh_connectivity[i,2]; 

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




function computeNodeStencilsSIMPLEX2D(nNodes::Int32, nNeibCells::Int32,
	mesh_nodes::Array{Float64,2}, cell_clusters::Array{Int32,2}, cell_mid_points::Array{Float64,2})::Array{Float64,2}

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

function cells2nodesConnectivity2D(nCells::Int32, stiff_matrix::Array{Int32,2}, mesh_connectivity::Array{Int32,2})::Array{Int32,2}

	cells2nodes = zeros(Int32,nCells,7);

for C=1:nCells

	nodesC = mesh_connectivity[C,4:8];
	num_nodes = mesh_connectivity[C,3];
	
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



# DO NOT USE, USE DISTRIBUTED VERSION 
# function computeCellStiffnessM2D(nCells,bc_indexes,bc_data,mesh_connectivity)::Array{Int64,2}

	# cell_stiffness = zeros(Int64, nCells,4);

	# #nbc = size(bc_indexes,1);
	# nbc = length(bc_indexes);
	# for m=1:nbc
		# i = bc_data[m,1];
		# k = bc_data[m,3];
		# cell_stiffness[i,k] = bc_indexes[m];				
	# end

	# p1::Int64 = 0;
	# p2::Int64 = 0;
	# p3::Int64 = 0;
	# p4::Int64 = 0;

	# k1::Int64 = 0;
	# k2::Int64 = 0;
	# k3::Int64 = 0;
	# k4::Int64 = 0;

	# fCells = [];


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
		
		 # ## check for those cell who has zero index in connectivity        
		 # if ( (et == 3) && (cell_stiffness[i,1] == 0 || cell_stiffness[i,2] == 0 || cell_stiffness[i,3] == 0)) 
			 # ##display("Somethnig wrong in stiffMatrix2d fun ...");
			 # #fCells = [fCells; i];
			 # push!(fCells,i);
		 # elseif ( (et == 2) && (cell_stiffness[i,1] == 0 || cell_stiffness[i,2] == 0 || cell_stiffness[i,3] == 0 || cell_stiffness[i,4] == 0)) 
			 # ##disp("Somethnig wrong in stiffMatrix2d fun ...#);
			 # #fCells = [fCells; i];
			 # push!(fCells,i);
			 
		 # end

		

	# end #i


	
	# n::Int64 = size(fCells,1);
 

	# for z = 1:n
		
		# i::Int64 = fCells[z];
		  
		# element_type::Int64 = mesh_connectivity[i,2];
		
		# p1 = mesh_connectivity[i,4];
		# p2 = mesh_connectivity[i,5];
		# p3 = mesh_connectivity[i,6];
		# p4 = mesh_connectivity[i,7];
			
		# zIndex::Int64 = -1000;
		
		# if (cell_stiffness[i,1] == 0)
			# zIndex = 1;
		# elseif (cell_stiffness[i,2] == 0)
			# zIndex = 2;
		# elseif (cell_stiffness[i,3] == 0)
			# zIndex = 3;
		# else
			# display("something wrong in cell stiffness calc ... ");
			# break;
		# end
				
		# for k = 1:nCells
						
			# cell_type::Int64 = mesh_connectivity[k,2];
			# num_nodes::Int64 = mesh_connectivity[k,3];
										
			# k1 = mesh_connectivity[k,4];
			# k2 = mesh_connectivity[k,5];
			# k3 = mesh_connectivity[k,6];
			# k4 = mesh_connectivity[k,7];

						
						
			# if ( cell_type == 2 && element_type == 3 )

			    # ##  find first edge
                # if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (i!=k) )
                    # cell_stiffness[i,zIndex] = k;
                # end
     			# ##  finding second edge
                # if ( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (i!=k) )
                    # cell_stiffness[i,zIndex] = k;
                # end
 				# ##  finding third edge
                # if ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (i!=k) )
                    # cell_stiffness[i,zIndex] = k;
                # end                   
						
			# end

		# end # end for k 
				
				
	 # end #%end for i 



	# return  cell_stiffness;
# end







function computeCells2Nodes2D(nCells::Int32, mesh_connectivity::Array{Int32,2}, cell_stiffness::Array{Int32,2}, cells2nodes::Array{Int32,2} )


  ##cells2nodes = zeros{Int64, nCells,8};
  

  for C = 1: nCells
  

     
      num_nodes::Int64 = mesh_connectivity[C,3];
      nodesC1::Int64 =   mesh_connectivity[C,4];
      nodesC2::Int64 =   mesh_connectivity[C,5];
      nodesC3::Int64 =   mesh_connectivity[C,6];
      nodesC4::Int64 =  mesh_connectivity[C,7];


      if (num_nodes == 3) ## triangle
      
        for T = 1:num_nodes
          

              neib_cell::Int64 = cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;

			  
			  nodesP = [];

              if (neib_cell >0) ## internal cell
              
				  nib_cell_num_nodes  = 	mesh_connectivity[neib_cell,3];
				  nodesT1::Int64 = mesh_connectivity[neib_cell,4];
				  nodesT2::Int64 = mesh_connectivity[neib_cell,5];
				  nodesT3::Int64 = mesh_connectivity[neib_cell,6];
				  nodesT4::Int64 = mesh_connectivity[neib_cell,7];
				  
				  if (nib_cell_num_nodes == 3) ## neib cell is triangle
				  
					  

					  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3)
						push!(nodesP, nodesC1);
						##nodesP.push_back(nodesC1);
					  end
					  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3)
						push!(nodesP, nodesC2);
						##nodesP.push_back(nodesC2);
					  end
					  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3)
						push!(nodesP, nodesC3);
						##nodesP.push_back(nodesC3);
					  end

					  if (size(nodesP,1) == 2)
					  
						node1 = nodesP[1];
						node2 = nodesP[2];
					  
					  else
					  
						  display("something wrong in creating cells2node matrix T1 ... ");
						  ##throw(-1);
					  end
				
					elseif (nib_cell_num_nodes == 4) ## neib cell is quad
					
					
							  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3 || nodesC1 == nodesT4)
								push!(nodesP,nodesC1);
								#nodesP.push_back(nodesC1);
							  end
							  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3 || nodesC2 == nodesT4)
								push!(nodesP,nodesC2);
								##nodesP.push_back(nodesC2);
							  end
							  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3 || nodesC3 == nodesT4)
								push!(nodesP,nodesC3); 
								##nodesP.push_back(nodesC3);
							  end
							  if (nodesC4 == nodesT1 ||  nodesC4 == nodesT2 || nodesC4 == nodesT3 || nodesC4 == nodesT4)
								push!(nodesP,nodesC4);
								##nodesP.push_back(nodesC4);
							  end


							  if (size(nodesP,1)  == 2)
							  
								node1 = nodesP[1];
								node2 = nodesP[2];
							  
							  else
							  
								  display("something wrong in creating cells2node matrix Q1... ");
								  ##throw(-1);
							  end
					
					
					end
				
              
              else ## boundary cell
              

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC1;
                  end


              end ##  end if 


              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
              
              elseif (T == 3)
              
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;
              
              else
              
                  display("something wrong in creating cells2node matrix T2... ");
                  ##throw(-1);
              end




		end ## end for particular triangular cell

		## end for triangle cells
	  
	  
      elseif (num_nodes == 4) ## quad element
      

          for T = 1:num_nodes
          

              neib_cell::Int64 = cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;


              nodesP = [];


              if (neib_cell >0) ## internal cell
              
                  nodesT1::Int64 = mesh_connectivity[neib_cell,4];
                  nodesT2::Int64 = mesh_connectivity[neib_cell,5];
                  nodesT3::Int64 = mesh_connectivity[neib_cell,6];
                  nodesT4::Int64 = mesh_connectivity[neib_cell,7];

                  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3 || nodesC1 == nodesT4)
				    push!(nodesP,nodesC1);
                    #nodesP.push_back(nodesC1);
				  end
                  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3 || nodesC2 == nodesT4)
				    push!(nodesP,nodesC2);
                    ##nodesP.push_back(nodesC2);
				  end
                  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3 || nodesC3 == nodesT4)
				    push!(nodesP,nodesC3); 
                    ##nodesP.push_back(nodesC3);
				  end
                  if (nodesC4 == nodesT1 ||  nodesC4 == nodesT2 || nodesC4 == nodesT3 || nodesC4 == nodesT4)
					push!(nodesP,nodesC4);
                    ##nodesP.push_back(nodesC4);
				  end


                  if (size(nodesP,1)  == 2)
                  
                    node1 = nodesP[1];
                    node2 = nodesP[2];
                  
                  else
                  
                      display("something wrong in creating cells2node matrix Q1... ");
                      ##throw(-1);
                  end

              
              else

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC4;
                  
                  elseif (T == 4)
                  
                       node1 = nodesC4;
                       node2 = nodesC1;
                  end


              end  ##  end boundary cell




              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
            
              elseif (T == 3)
			  
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;

              elseif (T == 4)
              
                  cells2nodes[C,7] = node1;
                  cells2nodes[C,8] = node2;

              else
              
                  display("something wrong in creating cells2node matrix Q2... ");
                  throw(-1);
              end

            end ## 
			
			## end for quad element

		else
      
          display("something wrong in creating cells2node matrix ... ");
          display("unknown element type: must be triangle or quad ... ");
		  
          ##throw(-1);
      end




  end ## end global loop for cells
  
  return cells2nodes;

  ## cout << "done " << endl;

end ## function



function computeCells2Nodes2Dhybrid(nCells::Int32, mesh_connectivity::Array{Int32,2}, cell_stiffness::Array{Int32,2}, cells2nodes::Array{Int32,2} )


  ##cells2nodes = zeros{Int64, nCells,8};
  

  for C = 1: nCells
  

     
      num_nodes::Int64 = mesh_connectivity[C,3];
      nodesC1::Int64 =   mesh_connectivity[C,4];
      nodesC2::Int64 =   mesh_connectivity[C,5];
      nodesC3::Int64 =   mesh_connectivity[C,6];
      nodesC4::Int64 =  mesh_connectivity[C,7];


      if (num_nodes == 3) ## triangle
      
        for T = 1:num_nodes
          

              neib_cell::Int64 = cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;

			  
			  nodesP = [];

              if (neib_cell >0) ## internal cell
              
				  neib_cell_num_nodes = mesh_connectivity[ neib_cell, 3];	
				  
				  nodesT1::Int64 = mesh_connectivity[neib_cell,4];
				  nodesT2::Int64 = mesh_connectivity[neib_cell,5];
				  nodesT3::Int64 = mesh_connectivity[neib_cell,6];
				  nodesT4::Int64 = mesh_connectivity[neib_cell,7];
			  
				  if (neib_cell_num_nodes == 3) 	
			  
					 

					  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3)
						push!(nodesP, nodesC1);
						##nodesP.push_back(nodesC1);
					  end
					  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3)
						push!(nodesP, nodesC2);
						##nodesP.push_back(nodesC2);
					  end
					  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3)
						push!(nodesP, nodesC3);
						##nodesP.push_back(nodesC3);
					  end

					  if (size(nodesP,1) == 2)
					  
						node1 = nodesP[1];
						node2 = nodesP[2];
					  
					  else
					  
						  display("something wrong in creating cells2node matrix T1 ... ");
						  ##throw(-1);
					  end
				  
				  elseif(neib_cell_num_nodes == 4)
				  
					  

					  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3 || nodesC1 == nodesT4)
						push!(nodesP,nodesC1);
						#nodesP.push_back(nodesC1);
					  end
					  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3 || nodesC2 == nodesT4)
						push!(nodesP,nodesC2);
						##nodesP.push_back(nodesC2);
					  end
					  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3 || nodesC3 == nodesT4)
						push!(nodesP,nodesC3); 
						##nodesP.push_back(nodesC3);
					  end
					  if (nodesC4 == nodesT1 ||  nodesC4 == nodesT2 || nodesC4 == nodesT3 || nodesC4 == nodesT4)
						push!(nodesP,nodesC4);
						##nodesP.push_back(nodesC4);
					  end


					  if (size(nodesP,1)  == 2)
					  
						node1 = nodesP[1];
						node2 = nodesP[2];
					  
					  else
					  
						display("something wrong in creating cells2node matrix Q1... ");
						##throw(-1);
					  end
				  
				  end

              
              else ## boundary cell
              

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC1;
                  end


              end ##  end if 


              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
              
              elseif (T == 3)
              
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;
              
              else
              
                  display("something wrong in creating cells2node matrix T2... ");
                  ##throw(-1);
              end




		end ## end for particular triangular cell

		## end for triangle cells
	  
	  
      elseif (num_nodes == 4) ## quad element
      

          for T = 1:num_nodes
          

              neib_cell::Int64 = cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;


              nodesP = [];


              if (neib_cell >0) ## internal cell
              
                  nodesT1::Int64 = mesh_connectivity[neib_cell,4];
                  nodesT2::Int64 = mesh_connectivity[neib_cell,5];
                  nodesT3::Int64 = mesh_connectivity[neib_cell,6];
                  nodesT4::Int64 = mesh_connectivity[neib_cell,7];

                  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3 || nodesC1 == nodesT4)
				    push!(nodesP,nodesC1);
                    #nodesP.push_back(nodesC1);
				  end
                  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3 || nodesC2 == nodesT4)
				    push!(nodesP,nodesC2);
                    ##nodesP.push_back(nodesC2);
				  end
                  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3 || nodesC3 == nodesT4)
				    push!(nodesP,nodesC3); 
                    ##nodesP.push_back(nodesC3);
				  end
                  if (nodesC4 == nodesT1 ||  nodesC4 == nodesT2 || nodesC4 == nodesT3 || nodesC4 == nodesT4)
					push!(nodesP,nodesC4);
                    ##nodesP.push_back(nodesC4);
				  end


                  if (size(nodesP,1)  == 2)
                  
                    node1 = nodesP[1];
                    node2 = nodesP[2];
                  
                  else
                  
                      display("something wrong in creating cells2node matrix Q1... ");
                      ##throw(-1);
                  end

              
              else

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC4;
                  
                  elseif (T == 4)
                  
                       node1 = nodesC4;
                       node2 = nodesC1;
                  end


              end  ##  end boundary cell




              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
            
              elseif (T == 3)
			  
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;

              elseif (T == 4)
              
                  cells2nodes[C,7] = node1;
                  cells2nodes[C,8] = node2;

              else
              
                  display("something wrong in creating cells2node matrix Q2... ");
                  throw(-1);
              end

            end ## 
			
			## end for quad element

		else
      
          display("something wrong in creating cells2node matrix ... ");
          display("unknown element type: must be triangle or quad ... ");
		  
          ##throw(-1);
      end




  end ## end global loop for cells
  
  return cells2nodes;

  ## cout << "done " << endl;

end ## function

