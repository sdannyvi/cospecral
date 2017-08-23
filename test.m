
function output = test(n,d)

%B is the adjacency matrix of a d-regular 2n-vertex random graph G=(V1 U
%V2,E). Rows are indexed by V1, columns by V2. 
B = creatRndG(n,d);

% Returns true if the generated graph is indeed d-regular, and then we can
% proceed.
while(B == zeros(n,n))   
    B = creatRndG(n,d);
end

% Now we use B and the Bipartite Graph Product to create to graphs which we
% know are co-spectral. We use two templates: XXX, and XYX. Finally we
% check whether the resulting graph are also non-isomporhic. 


% create XXX vs YYY matrix
rowIndx = 0;
C1 = zeros(n^3,n^3);
for ix1=1:n
    for ix2=1:n
        for ix3=1:n
        rowIndx=rowIndx+1;
        for iy1=1:n
            for iy2=1:n              
                for iy3=1:n              
                if (B(ix1,iy1) == 1) || (B(ix2,iy2) == 1) || ( B(ix3,iy3) == 1) 
                    colIndx = n^2*(iy1-1)+n*(iy2-1)+iy3;    
                    C1(rowIndx,colIndx)=1;
                end  
                end                    
            end
        end
        end
    end
end



A1 = horzcat(zeros(n^3,n^3),C1);
A2 = horzcat(C1',zeros(n^3,n^3));
Cxxyy = vertcat(A1,A2);


% create XYX vs YXY matrix
rowIndx1 = 0;
C2 = zeros(n^3,n^3);
for ix1=1:n
    for iy2=1:n
        for ix3=1:n
        rowIndx1=rowIndx1+1;
        for iy1=1:n
            for ix2=1:n              
                for iy3=1:n              
                if (B(ix1,iy1) == 1) || (B(ix2,iy2) == 1) || ( B(ix3,iy3) == 1) 
                    colIndx1 = n^2*(iy1-1)+n*(ix2-1)+iy3;    
                    C2(rowIndx1,colIndx1)=1;
                end  
                    
            end
        end
        end
        end
    end
end




A1 = horzcat(zeros(n^3,n^3),C2);
A2 = horzcat(C2',zeros(n^3,n^3));
Cxyyx = vertcat(A1,A2);

%check if the graphs are not isomoprhic by looking at the sequence of
%number of closed paths of length 4.

degSeq2 = sort(diag(Cxyyx^4));
degSeq1 = sort(diag(Cxxyy^4));


if (any(abs(degSeq2-degSeq1)>1))
   ' not isomorphic ...'
    B
end





% [V, Lambda1] = eig(Cxxyy); 
% Lambda1 = diag(Lambda1); 
% sort(Lambda1,'descend');
% 
% 
% [V, Lambda2] = eig(Cxyyx); 
% Lambda2 = diag(Lambda2); 
% sort(Lambda2,'descend');
% 
% (Lambda1-Lambda2)'*(Lambda1-Lambda2)
% 
% 
% if(Lambda1==Lambda2)
%     1
% else 0
% end


end

function output = creatRndG(vertNum, deg)

A=zeros(vertNum,vertNum);

ydeg = deg*ones(1,vertNum);

for x=1:vertNum
        for cur=1:deg      
 
            y = ceil(rand*vertNum);
            max_rep = 10;        
            while( (A(x,y)==1 || ydeg(y)==0) && max_rep>0)
                y = ceil(rand*vertNum);
                max_rep = max_rep-1;
            end
            
            if(max_rep == 0)
                y=-1;
                for z=1:vertNum
                    if(A(x,z)==0 && ydeg(z)>0)   
                        y=z;                
                    end
                end
            end
            
            if(y==-1)
                output = zeros(vertNum,vertNum);
                return;
            end
                                      
            A(x,y)=1;
            ydeg(y)=ydeg(y)-1;          
        end
        
        
end    
        output = A;

end 
                       


