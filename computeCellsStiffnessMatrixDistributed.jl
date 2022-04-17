
function distibuteCellsInThreadsSA(nThreads::Int32, nCells::Int32) ::SharedArray{Int32,2}

	cellsThreads = SharedArray{Int32}(Int64(nThreads),2);
 
 	if (nThreads>1)
		
    	nParts = floor(nCells/nThreads);

	    for i=1:nThreads
    		cellsThreads[i,2] =  nCells - nParts*(nThreads-i );
		end
	

    	for i=1:nThreads
      		cellsThreads[i,1] =  cellsThreads[i,2] - nParts + 1;
		end

    	cellsThreads[1,1] = 1;

		#display(cellsThreads);	
		
	end
	
	return cellsThreads;

end





## nCells, nThreads, bc_indexes, bc_data, mesh_connectivity

function computeCellStiffnessDistributed(nCells::Int32, nThreads::Int32,
		bc_indexes::Array{Int32}, bc_data::Array{Int32,2}, mesh_connectivity::Array{Int32,2})
	
	debug::Bool = false;
	
	if debug
		display("1");
	end
	
	cell_stiffness = zeros(Int32,nCells,4);
	cell_stiffnessSA = SharedArray{Int32}(Int64(nCells),4);
	mesh_connectivitySA = SharedArray{Int32}(Int64(nCells),7);
	
	for i = 1:nCells
		
		mesh_connectivitySA[i,1] = mesh_connectivity[i,1];
		mesh_connectivitySA[i,2] = mesh_connectivity[i,2];
		mesh_connectivitySA[i,3] = mesh_connectivity[i,3];
		mesh_connectivitySA[i,4] = mesh_connectivity[i,4];
		mesh_connectivitySA[i,5] = mesh_connectivity[i,5];
		mesh_connectivitySA[i,6] = mesh_connectivity[i,6];
		
		if size(mesh_connectivity,2) == 7
			mesh_connectivitySA[i,7] = mesh_connectivity[i,7];
		end
		
	end
	
	if debug
		display("2");
	end


	#nbc = size(bc_indexes,1);
	nbc = length(bc_indexes);
	for m=1:nbc
		i = bc_data[m,1];
		k = bc_data[m,3];
		cell_stiffnessSA[i,k] = bc_indexes[m];				
	end

	p1::Int32 = 0;
	p2::Int32 = 0;
	p3::Int32 = 0;
	p4::Int32 = 0;

	k1::Int32 = 0;
	k2::Int32 = 0;
	k3::Int32 = 0;
	k4::Int32 = 0;

	cellsThreads = distibuteCellsInThreadsSA(nThreads, nCells); ## partition mesh 
	
	if debug
		display("3");
	end

	
	@everywhere nCellsX = $nCells;
	
	mixedCellsX = SharedVector{Int32}(Int64(nCells));
	
	@sync @distributed for p in workers()	
	
		beginCell::Int32 = cellsThreads[p-1,1];
		endCell::Int32 = cellsThreads[p-1,2];	
		kernelStiffness(nCellsX, beginCell, endCell, mesh_connectivitySA, cell_stiffnessSA, mixedCellsX)	
						
	end
		
	@everywhere finalize(kernelStiffness);		
	
	if debug
		display("4");
	end


	mixedCells = [];

	 for i = 1:nCellsX
	 
		if (mixedCellsX[i] != 0)
			push!(mixedCells, mixedCellsX[i]);
		end
		
	 end
	
	
	n::Int32 = size(mixedCells,1);
 

	for z = 1:n
		
		i::Int32 = mixedCells[z];
		  
		element_type::Int32 = mesh_connectivitySA[i,2];
		
		p1 = mesh_connectivitySA[i,4];
		p2 = mesh_connectivitySA[i,5];
		p3 = mesh_connectivitySA[i,6];
		p4 = mesh_connectivitySA[i,7];
			
		zIndex::Int32 = -1000;
		
		if (cell_stiffnessSA[i,1] == 0)
			zIndex = 1;
		elseif (cell_stiffnessSA[i,2] == 0)
			zIndex = 2;
		elseif (cell_stiffnessSA[i,3] == 0)
			zIndex = 3;
		else
			display("something wrong in cell stiffness calc ... ");
			break;
		end
				
		for k = 1:nCells
						
			cell_type::Int32 = mesh_connectivitySA[k,2];
			num_nodes::Int32 = mesh_connectivitySA[k,3];
										
			k1 = mesh_connectivitySA[k,4];
			k2 = mesh_connectivitySA[k,5];
			k3 = mesh_connectivitySA[k,6];
			k4 = mesh_connectivitySA[k,7];

						
						
			if ( cell_type == 2 && element_type == 3 )

			    ##  find first edge
                if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (i!=k) )
                    cell_stiffnessSA[i,zIndex] = k;
                end
     			##  finding second edge
                if ( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (i!=k) )
                    cell_stiffnessSA[i,zIndex] = k;
                end
 				##  finding third edge
                if ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (i!=k) )
                    cell_stiffnessSA[i,zIndex] = k;
                end                   
						
			end

		end # end for k 
				
				
	 end #%end for i 


	if debug
		display("5");
	end
	
	
	for i = 1:nCellsX
		cell_stiffness[i,1] = cell_stiffnessSA[i,1];
		cell_stiffness[i,2] = cell_stiffnessSA[i,2];
		cell_stiffness[i,3] = cell_stiffnessSA[i,3];
		cell_stiffness[i,4] = cell_stiffnessSA[i,4];
	end



	return  cell_stiffness, mesh_connectivitySA; 
	
end


@everywhere function kernelStiffness(nCells::Int32, beginCell::Int32, endCell::Int32,
	mesh_connectivity::SharedArray{Int32,2}, cell_stiffness::SharedArray{Int32,2}, fCells::SharedVector{Int32})
	
	#fCells = [];
	
	for i = beginCell:endCell

		et::Int32 = mesh_connectivity[i,2];
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
			 #push!(fCells,i);
			 fCells[i] = i;
		 elseif ( (et == 2) && (cell_stiffness[i,1] == 0 || cell_stiffness[i,2] == 0 || cell_stiffness[i,3] == 0 || cell_stiffness[i,4] == 0)) 
			 ##disp("Somethnig wrong in stiffMatrix2d fun ...#);
			 #fCells = [fCells; i];
			 #push!(fCells,i);
			 
			 fCells[i] = i;
			 
		 end

		

	end #i


	#println(fCells)
	
	
end