%%% collects the position vectors for each array element, considering HFFS
%%% model's lumped ports annotation 

clc;
nbrPorts = 40;
%x
orderX{1} = [23 24 21 22];
orderX{2} = [27 31 35 39 7];
orderX{3} = [11 12 9 10];
orderX{4} = [28 32 36 40 8];
orderX{5} = [15 16 13 14];
orderX{6} = [25 29 33 37 5];
orderX{7} = [19 20 17 18];
orderX{8} = [26 30 34 38 6];
orderX{9} = [3 4 1 2];
%y
orderY{1} = [26 25 28 27 ];
orderY{2} = [3 19 15 11 23];
orderY{3} = [30 29 32 31];
orderY{4} = [4 20 16 12 24];
orderY{5} = [34 33 36 35];
orderY{6} = [1 17 13 9 21];
orderY{7} = [38 37 40 39];
orderY{8} = [2 18 14 10 22];
orderY{9} = [6 5 8 7];

equiRowX = zeros(1,nbrPorts);
equiRowY = zeros(1,nbrPorts);
for k=1:nbrPorts
    for i=1:9
        tmpX = orderX{i};
        tmpY = orderY{i};
        for j=1:size(tmpX,2)
          if(k == tmpX(j))
            equiRowX(k) = i;
          end
        end
        for j=1:size(tmpY,2)
          if(k == tmpY(j))
            equiRowY(k) = i;
          end
        end  
    end    
end
fprintf('\n');
fprintf('%d ', equiRowX);
fprintf('\n');
fprintf('%d ', equiRowY);
fprintf('\n');