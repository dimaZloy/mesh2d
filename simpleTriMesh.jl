


display("reading mesh data...");

display("read mesh nodes ... ");
include( "p.jl");
nNodes = size(mesh_nodes,1); # number of nodes;

display("read mech connectivity ... ");
include( "c.jl" );
nCells = size(mesh_connectivityTri,1); #number of cells


mesh_connectivity = zeros(Int64, nCells, 7);
for f=1:nCells
    mesh_connectivity[f,1:6] = 	mesh_connectivityTri[f,1:6];	
end

# mesh_connectivity has the following format (nCells x 7):
# 1 element - global id,
# 2 element - cell type (2 for quads /3 for triangles)
# 3 element - number of nodes
# 4-7 - nodes ids
# 3 nodes for tri mesh 
# 4 nodes for quad mesh 

display("read mesh boundaries ... ");
include( "b.jl" );


#bc_indexes[1:size(bctop,1)] = -1;
#bc_indexes[nBCtop+1:nBCtop+nBCright] = -2;
#bc_indexes[nBCtop+nBCright+1:nBCtop+nBCright+nBCbottom] = -3;


bc_data = [bctop; bcright; bcbottom; bcleft;];

nBC = size(bc_data,1);
nBSets = nBC;
nBCtop = size(bctop,1);
nBCright = size(bcright,1);
nBCbottom = size(bcbottom,1);


bc_indexes = ones(Int64, nBC)*-4;
for i=1:nBCtop
	bc_indexes[i] = -1;
end
for i=nBCtop+1 : nBCtop+nBCright
	bc_indexes[i] = -2;
end
for i=nBCtop+nBCright+1 : nBCtop+nBCright+nBCbottom
	bc_indexes[i] = -3;
end

