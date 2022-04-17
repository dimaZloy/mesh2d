

using Distributed;
using PyPlot;

using WriteVTK;
using CPUTime;
using DelimitedFiles;
using Printf
using BSON: @load
using BSON: @save
using SharedArrays;

using HDF5;
using ProfileView;


@everywhere struct mesh2d_Int32
	nCells::Int64
	nNodes::Int64
	nNeibCells::Int32						## max number of neighbors 
	nBSets::Int32							##  number of boundaries  
	xNodes::Array{Float64,1} 				##  mesh_nodes[nNodesx3]
	yNodes::Array{Float64,1} 				##	mesh_nodes[nNodesx3]
	mesh_connectivity::Array{Int32,2} 		## [nCellsx7]
	bc_data::Array{Int32,2}
	bc_indexes::Array{Int32,1}
	cell_nodes_X::Array{Float64,2} 			## [nCellsx4]
	cell_nodes_Y::Array{Float64,2} 			## [nCellsx4]
	cell_mid_points::Array{Float64,2} 		## [nCellsx2]
	cell_areas::Array{Float64,1} 			## [nCellsx1]
	HX::Array{Float64,1}
	cell_edges_Nx::Array{Float64,2} 		## [nCellsx4]
	cell_edges_Ny::Array{Float64,2} 		## [nCellsx4]
	cell_edges_length::Array{Float64,2} 	## [nCellsx4]
	cell_stiffness::Array{Int32,2} 			## [nCellsx4]
	cell_clusters::Array{Int32,2} 			## [nNodesx8]
	node_stencils::Array{Float64,2} 		## [nNodesx8]
	node2cellsL2up::Array{Int32,2} 			## [nCellsx8]
	node2cellsL2down::Array{Int32,2} 		## [nCellsx8]
	cells2nodes::Array{Int32,2}
	VTKCells::Array{MeshCell,1}
	triangles::Array{Int32,2}
end


include("utilsIO.jl");


function plot_cell(cell, tm, lw, color )


	cell_type = tm.mesh_connectivity[cell,2];	
	cell_num_nodes = tm.mesh_connectivity[cell,3];	
		
	p1 = tm.mesh_connectivity[cell,4];	
	p2 = tm.mesh_connectivity[cell,5];	
	p3 = tm.mesh_connectivity[cell,6];	
	p4 = tm.mesh_connectivity[cell,7];	
		
	if (cell_num_nodes == 3) 
		plot([tm.xNodes[p1], tm.xNodes[p2], tm.xNodes[p3], tm.xNodes[p1] ], [tm.yNodes[p1], tm.yNodes[p2], tm.yNodes[p3], tm.yNodes[p1]], color,linewidth = lw )
	elseif (cell_num_nodes == 4) 
		plot([tm.xNodes[p1], tm.xNodes[p2], tm.xNodes[p3], tm.xNodes[p4], tm.xNodes[p1] ], [tm.yNodes[p1], tm.yNodes[p2], tm.yNodes[p3], tm.yNodes[p4], tm.yNodes[p1]], color,linewidth = lw )
	end
	
end


function test()

	tm = readMesh2dHDF5("cyl2d_supersonic1_BL_test");

    display(tm.nCells)
	display(tm.nNodes)

	ms = 2.5;
	

	display(tm.cells2nodes)
	display( findall(x->x==0, tm.cells2nodes) )
	

	figure(1)
	clf()

	#cell_stiffness
	#mesh_connectivity

	plot(tm.xNodes, tm.yNodes,"ok",markersize=ms);

	axis("equal")
	
	cell = 1003;
	display(tm.cells2nodes[cell,:])
	display(tm.mesh_connectivity[cell,2])
	display(tm.mesh_connectivity[cell,3])
	
	#cell = 1800;

	plot_cell(cell,tm, 1.5, "-g");
	
	#cells2nodes 
	

	for i = 1:2: tm.nNeibCells
	
		p1 = tm.cells2nodes[cell,i]
		p2 = tm.cells2nodes[cell,i+1]

		if (p1 > 0 && p2 > 0)
			plot([tm.xNodes[p1], tm.xNodes[p2]], [tm.yNodes[p1], tm.yNodes[p2]], "sb", markersize = ms+1 )			
		end

	end
	

	
	# cell clusters is okay
	# for nodeI = 1:4
	
		# node = tm.mesh_connectivity[cell,3+ nodeI];
		
		
		# for j = 1:tm.nNeibCells
			# if node > 0
				# celli = tm.cell_clusters[node, j];
				# if celli > 0
					# plot_cell(celli,tm, 1.0, "-.c");
				# end
			# end	 
		# end
		
	# end

	
	# # cell stiffnes is okay 	
	 # for neibi = 1:4
		
		 # celli = tm.cell_stiffness[cell,neibi];
		
		 # if celli > 0
			 # plot_cell(celli,tm, 1.0, "-.r");
		 # end
		
	 # end
	
	
	
		
	
	#end
	
end

test()