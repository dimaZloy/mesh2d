

function distibuteNodesInThreadsSA(nThreads::Int64, nNodes::Int64 )::SharedArray{Int64,2}

	nodesThreads = SharedArray{Int64}(nThreads,2);
 
 	if (nThreads>1)
		
    	nParts = floor(nNodes/nThreads);

	    for i=1:nThreads
    		nodesThreads[i,2] =  nodes - nParts*(nThreads-i );
		end
	

    	for i=1:nThreads
      		nodesThreads[i,1] =  nodesThreads[i,2] - nParts + 1;
		end

    	nodesThreads[1,1] = 1;

		#display(cellsThreads);	
		
	end
	
	return nodesThreads;

end

function computeCellClustersDistributed(nNodes::Int64, nCells::Int64, nNeibCells::Int64, 
		mesh_connectivitySA::SharedArray{Int64,2})::Array{Int64,2}


	nodesThreads = distibuteCellsInThreadsSA(numThreads, nNodes); ## partition mesh 
	#display(nodesThreads)
	
	cell_clustersSA = SharedArray{Int64}(nNodes,nNeibCells);
	
	@everywhere nCellsX = $nCells;
	@everywhere nNeibCellsX = $nNeibCells;
	
	
	@sync @distributed for p in workers()	
	
		beginNode::Int64 = nodesThreads[p-1,1];
		endNode::Int64 = nodesThreads[p-1,2];		
		# println("worker: ",p, "\tbegin node: ",beginNode,"\tend node: ", endNode,"\tnCells: ", nCellsX,"\tnNeibCells: ", nNeibCellsX  );
		
		kernelCellClusters(beginNode, endNode, nCellsX, nNeibCellsX,  mesh_connectivitySA, cell_clustersSA)				
	end
		
	@everywhere finalize(kernelCellClusters);		
	
	
	cell_clusters = zeros(Int64, nNodes, nNeibCells );
	
	for i = 1:nNodes
		cell_clusters[i,:] = cell_clustersSA[i,:];
	end
	
	#display(cell_clusters)
	
	return cell_clusters;

end



@everywhere function kernelCellClusters(beginNode::Int64, endNode::Int64, nCells::Int64, nNeibCells::Int64, 
	mesh_connectivity::SharedArray{Int64,2}, cell_clusters::SharedArray{Int64,2})

	#cluster_size = 7;
	#cell_clusters =zeros(Int64, nNodes, nNeibCells);


	for i=beginNode:endNode

		c_node = i;
		counter = 0; 

		#cluster = zeros(1,cluster_size);

		cluster = zeros(Int64,nNeibCells);

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

	

end
