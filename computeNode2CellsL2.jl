

function computeNode2CellsL2(nCells::Int32, mesh_connectivity::Array{Int32,2}, cell_stiffness::Array{Int32,2},
	node2cellL2up::Array{Int32,2},node2cellL2down::Array{Int32,2})


	##node2cellL2up = size(4x2)
	##node2cellL2up = size(4x2)
	## for triangular control volume: use idexes 1,3,6
	## for quad control volume: use idexes 1,2; 3,4; 5;6 - then take average value 
	
	
	node2cellL2up1= zeros(Int32,nCells,4);
	node2cellL2up2= zeros(Int32,nCells,4);
	
	node2cellL2down1= zeros(Int32,nCells,4);
	node2cellL2down2= zeros(Int32,nCells,4);
	

for i=1:nCells, 
    
    ck::Int32 = mesh_connectivity[i,3]; ##  number of nodes (edges) in the i-cell
		
	p1::Int32 = 0;
	p2::Int32 = 0;
	p3::Int32 = 0;
	p4::Int32 = 0;
	
	k1::Int32  = 0;
	k2::Int32  = 0;
	k3::Int32  = 0;
	k4::Int32  = 0;
				
	pDown1::Int32 = 0;
	pDown2::Int32 = 0;
				
	pUp1::Int32 = 0;
	pUp2::Int32 = 0;
	
       
    if (ck == 3)  ##  triangle cell with 3 nodes        
	
        p1 = mesh_connectivity[i,4];
        p2 = mesh_connectivity[i,5];
        p3 = mesh_connectivity[i,6];

    elseif (ck == 4) ## quad cell with 4 nodes
        
        p1 = mesh_connectivity[i,4];
        p2 = mesh_connectivity[i,5];
        p3 = mesh_connectivity[i,6];
        p4 = mesh_connectivity[i,7];
        
    else
        println("something wrong in inviscid flux ... ");
        break;
    end

    
    for k = 1:ck ## for each node 
        
			ek::Int32 = cell_stiffness[i,k]; ## get right cell 
			
			
            
            if (ek >=1 && ek<=nCells)
                
				# k1::Int32  = 0;
				# k2::Int32  = 0;
				# k3::Int32  = 0;
				# k4::Int32  = 0;
				
				# pDown1::Int32 = 0;
				# pDown2::Int32 = 0;
				
				# pUp1::Int32 = 0;
				# pUp2::Int32 = 0;
				
                 ekNodes::Int32 = mesh_connectivity[ek,3];  ## get num nodes in neib cell
                
                if (ck == 3 )  ##  tri-tri
                     
                     k1 = mesh_connectivity[ek,4];
                     k2 = mesh_connectivity[ek,5];
                     k3 = mesh_connectivity[ek,6];
                     k4 = mesh_connectivity[ek,7];
                     
                     if (ekNodes == 3) ## tri
                                             
						 ## find first edge
						 if ( (p1 == k1 || p1 == k2 || p1 == k3 ) && (p2 == k1 || p2 == k2 || p2 == k3 ) )
							 pDown1 = p3;
						 end
						 ## finding second edge
						 if ( (p2 == k1 || p2 == k2 || p2 == k3 ) && (p3 == k1 || p3 == k2 || p3 == k3 )  )
							  pDown1 = p1;
						 end
						 ##  finding third edge
						 if ( (p3 == k1 || p3 == k2 || p3 == k3 ) && (p1 == k1 || p1 == k2 || p1 == k3 )  )
							  pDown1 = p2;
						 end                   
						 
						  ## find first edge
						 if ( (k1 == p1 || k1 == p2 || k1 == p3 ) && (k2 == p1 || k2 == p2 || k2 == p3 ) )
							 pUp1 = k3;
						 end
						 ## finding second edge
						 if ( (k2 == p1 || k2 == p2 || k2 == p3 ) && (k3 == p1 || k3 == p2 || k3 == p3 )  )
							  pUp1 = k1;
						 end
						##  finding third edge
						 if ( (k3 == p1 || k3 == p2 || k3 == p3 ) && (k1 == p1 || k1 == p2 || k1 == p3 )  )
							  pUp1 = k2;
						 end                   

                    
                     elseif (ekNodes == 4) ## quad
                         
                        if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) ) # && (i!=k) )
                   
							pDown1 = p3;
							pDown2 = p4;
						elseif( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) ) # && (i!=k) )
                     
							pDown1 = p1;
							pDown2 = p4;
 
						elseif ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) ) #&& (i!=k) )
                     
							pDown1 = p1;
							pDown2 = p2;
 
						elseif ( (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) ) #&& (i!=k) )
                     
							pDown1 = p2;
							pDown2 = p3;
 
						end                   
                                          
						if ( (k1 == p1 || k1 == p2 || k1 == p3 || k1 == p4) && (k1 == p1 || k2 == p2 || k2 == p3 || k2 == p4) )# && (i!=k) )
							pUp1 = k3;
							pUp2 = k4;
						elseif ( (k2 == p1 || k2 == p2 || k2 == p3 || k2 == p4) && (k3 == p1 || k3 == p2 || k3 == p3 || k3 == p4) ) # && (i!=k) )
                  
							pUp1 = k1;
							pUp2 = k4;
                     
						elseif ( (k3 == p1 || k3 == p2 || k3 == p3 || k3 == p4) && (k4 == p1 || k4 == p2 || k4 == p3 || k4 == p4) )# && (i!=k) )
							pUp1 = k1;
							pUp2 = k2;
 
						elseif ( (k4 == p1 || k4 == p2 || k4 == p3 || k4 == p4) && (k1 == p1 || k1 == p2 || k1 == p3 || k1 == p4) )#&& (i!=k) )
                               
							pUp1 = k2;
							pUp2 = k3;
 
						end                   

                     end
                     
                                       
                   #node2cellL2up[i,k] = pUp1;
                   #node2cellL2down[i,k] = pDown1;
				   
				   node2cellL2up1[i,k] = pUp1;
				   node2cellL2up2[i,k] = pUp2;
				   node2cellL2down1[i,k] = pDown1;
				   node2cellL2down2[i,k] = pDown2;
				                       
                    
                elseif ( ck == 4 )  ##  quad-quad
                    
                    k1 = mesh_connectivity[ek,4];
                    k2 = mesh_connectivity[ek,5];
                    k3 = mesh_connectivity[ek,6];
                    k4 = mesh_connectivity[ek,7];
         
                    if (ekNodes == 4) ## quad
                     
						if ( (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) && (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) ) # && (i!=k) )
					   
							 pDown1 = p3;
							 pDown2 = p4;
						elseif( (p2 == k1 || p2 == k2 || p2 == k3 || p2 == k4) && (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) ) # && (i!=k) )
						 
							 pDown1 = p1;
							 pDown2 = p4;
	 
						elseif ( (p3 == k1 || p3 == k2 || p3 == k3 || p3 == k4) && (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) ) #&& (i!=k) )
						 
							 pDown1 = p1;
							 pDown2 = p2;
	 
						elseif ( (p4 == k1 || p4 == k2 || p4 == k3 || p4 == k4) && (p1 == k1 || p1 == k2 || p1 == k3 || p1 == k4) ) #&& (i!=k) )
						 
							 pDown1 = p2;
							 pDown2 = p3;
	 
						end                   
                                          
						if ( (k1 == p1 || k1 == p2 || k1 == p3 || k1 == p4) && (k1 == p1 || k2 == p2 || k2 == p3 || k2 == p4) )# && (i!=k) )
							 pUp1 = k3;
							 pUp2 = k4;
						elseif ( (k2 == p1 || k2 == p2 || k2 == p3 || k2 == p4) && (k3 == p1 || k3 == p2 || k3 == p3 || k3 == p4) ) # && (i!=k) )
					  
							 pUp1 = k1;
							 pUp2 = k4;
						 
						elseif ( (k3 == p1 || k3 == p2 || k3 == p3 || k3 == p4) && (k4 == p1 || k4 == p2 || k4 == p3 || k4 == p4) )# && (i!=k) )
							 pUp1 = k1;
							 pUp2 = k2;
	 
						elseif ( (k4 == p1 || k4 == p2 || k4 == p3 || k4 == p4) && (k1 == p1 || k1 == p2 || k1 == p3 || k1 == p4) )#&& (i!=k) )
								   
							 pUp1 = k2;
							 pUp2 = k3;
	 
						end                   

                    elseif (ekNodes == 3) ## tri
                         
                        #% [pDown1, pUp1] = findEdge3x(p1,p2,p3,k1,k2,k3);
                         
                        #  % find first edge
						if ( (p1 == k1 || p1 == k2 || p1 == k3 ) && (p2 == k1 || p2 == k2 || p2 == k3 ) )
							pDown1 = p3;
						end
						#% finding second edge
						if ( (p2 == k1 || p2 == k2 || p2 == k3 ) && (p3 == k1 || p3 == k2 || p3 == k3 )  )
							pDown1 = p1;
						end
						#  finding third edge
						if ( (p3 == k1 || p3 == k2 || p3 == k3 ) && (p1 == k1 || p1 == k2 || p1 == k3 )  )
                          pDown1 = p2;
						end                   
                     
						#find first edge
						if ( (k1 == p1 || k1 == p2 || k1 == p3 ) && (k2 == p1 || k2 == p2 || k2 == p3 ) )
							pUp1 = k3;
						end
						# finding second edge
						if ( (k2 == p1 || k2 == p2 || k2 == p3 ) && (k3 == p1 || k3 == p2 || k3 == p3 )  )
                          pUp1 = k1;
						end
						# finding third edge
						if ( (k3 == p1 || k3 == p2 || k3 == p3 ) && (k1 == p1 || k1 == p2 || k1 == p3 )  )
                          pUp1 = k2;
						end                   
					
					end
                  
					node2cellL2up1[i,k] = pUp1;
					node2cellL2up2[i,k] = pUp2;
					node2cellL2down1[i,k] = pDown1;
					node2cellL2down2[i,k] = pDown2;


                else
                    println("something wrong in computeNode2CellsL2 matrixes ... ");
                    break;
                end
                                         
                   
            else
			
            end
            
         
            
		end ## k - for each cell edge   
		
	end ## i - loop for all cells
	
	
	
	for i = 1:nCells
	
		node2cellL2up[i,1] = node2cellL2up1[i,1];
		node2cellL2up[i,2] = node2cellL2up2[i,1];
		
		node2cellL2up[i,3] = node2cellL2up1[i,2];
		node2cellL2up[i,4] = node2cellL2up2[i,2];
				
		node2cellL2up[i,5] = node2cellL2up1[i,3];
		node2cellL2up[i,6] = node2cellL2up2[i,3];

		node2cellL2up[i,7] = node2cellL2up1[i,4];
		node2cellL2up[i,8] = node2cellL2up2[i,4];
		

		node2cellL2down[i,1] = node2cellL2down1[i,1];
		node2cellL2down[i,2] = node2cellL2down2[i,1];
		
		node2cellL2down[i,3] = node2cellL2down1[i,2];
		node2cellL2down[i,4] = node2cellL2down2[i,2];
				
		node2cellL2down[i,5] = node2cellL2down1[i,3];
		node2cellL2down[i,6] = node2cellL2down2[i,3];

		node2cellL2down[i,7] = node2cellL2down1[i,4];
		node2cellL2down[i,8] = node2cellL2down2[i,4];
		
		
		
		
	end
	
	
	
	

end ## end function