bcleft  =[#     1      20       0       6
	6    3    2
        20    3    2
        32    3    2
        42    3    2
        68    3    2
        82    3    2
       100    3    2
       118    3    2
       122    3    2
       150    3    2
       142    3    3
       130    3    3
       110    3    3
        94    3    3
        74    3    3
        62    3    3
        52    3    3
        28    3    3
        16    3    3
         1    3    1
];

bctop=[#       1      82       0       6
	8    3    2
        10    3    2
        38    3    2
        46    3    2
        58    3    2
        86    3    2
        90    3    2
       106    3    2
       128    3    2
       140    3    2
       160    3    2
       164    3    2
       172    3    2
       184    3    2
       188    3    2
       306    3    2
       308    3    2
       310    3    2
       312    3    2
       314    3    2
       326    3    2
       328    3    2
       330    3    2
       332    3    2
       334    3    2
       376    3    2
       378    3    2
       380    3    2
       382    3    2
       384    3    2
       400    3    3
       398    3    3
       396    3    3
       394    3    3
       364    3    3
       362    3    3
       360    3    3
       358    3    3
       356    3    3
       304    3    3
       302    3    3
       300    3    3
       298    3    3
       288    3    3
       286    3    3
       284    3    3
       282    3    3
       280    3    3
       268    3    3
       266    3    3
       264    3    3
       262    3    3
       260    3    3
       248    3    3
       246    3    3
       244    3    3
       242    3    3
       240    3    3
       228    3    3
       226    3    3
       220    3    3
       218    3    3
       216    3    3
       214    3    3
       212    3    3
       200    3    3
       196    3    3
       192    3    3
       180    3    3
       178    3    3
       168    3    3
       156    3    3
       146    3    3
       134    3    3
       114    3    3
        98    3    3
        78    3    3
        66    3    3
        50    3    3
        26    3    3
        14    3    3
         5    3    1
];
bcright=[#       1      20       0       6
	4    3    2
        24    3    2
        34    3    2
        56    3    2
        70    3    2
        80    3    2
       102    3    2
       116    3    2
       136    3    2
       148    3    2
       152    3    3
       126    3    3
       120    3    3
       104    3    3
        84    3    3
        72    3    3
        44    3    3
        36    3    3
        22    3    3
         7    3    1
];
bcbottom=[#       1      82       0       6
	2    3    2
        18    3    2
        30    3    2
        54    3    2
        64    3    2
        76    3    2
        96    3    2
       112    3    2
       132    3    2
       144    3    2
       154    3    2
       166    3    2
       174    3    2
       176    3    2
       190    3    2
       194    3    2
       198    3    2
       202    3    2
       204    3    2
       206    3    2
       208    3    2
       210    3    2
       222    3    2
       224    3    2
       230    3    2
       232    3    2
       234    3    2
       236    3    2
       238    3    2
       250    3    2
       252    3    2
       254    3    2
       256    3    2
       258    3    2
       270    3    2
       272    3    2
       274    3    2
       276    3    2
       278    3    2
       290    3    2
       292    3    2
       294    3    2
       296    3    2
       346    3    2
       348    3    2
       350    3    2
       352    3    2
       354    3    2
       386    3    2
       388    3    2
       390    3    2
       392    3    2
       374    3    3
       372    3    3
       370    3    3
       368    3    3
       366    3    3
       344    3    3
       342    3    3
       340    3    3
       338    3    3
       336    3    3
       324    3    3
       322    3    3
       320    3    3
       318    3    3
       316    3    3
       186    3    3
       182    3    3
       170    3    3
       162    3    3
       158    3    3
       138    3    3
       124    3    3
       108    3    3
        92    3    3
        88    3    3
        60    3    3
        48    3    3
        40    3    3
        12    3    3
         3    3    1
];


bc_data = [bctop; bcright; bcbottom; bcleft;];

nBC = size(bc_data,1);
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
