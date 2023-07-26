%% non-uniform mesh for the fractional Laplacian
%% Based on the Oberman's algorithm
%% Introduced non-uniform mesh and revise the original scheme
%% Two scale: H and h

function A = fractional1Dmat_adap(para,node)

%% primitive function F
if para.s == 0.5
    F = @(x) -para.Cs*log(abs(x)+eps);
    DF =@(x) -sign(x).*para.Cs./(abs(x) + eps);
else
    F = @(x) para.Cs*(abs(x)+eps).^(1-2*para.s)/2/para.s/(2*para.s - 1);
    DF =@(x) -para.Cs*sign(x).*((abs(x)+eps).^(-2*para.s))/para.s/2;
end


%% node & mesh

NV =  size(node,2);             %% vertex number
NT = NV -1;                     %% element number

node1 = -3: 2/NT: -1;
node2 = 1:2/NT:3;

NV2 = size(node1,2);
elem = [1:NV-1;2:NV]';          %% element 

NT2 = NV2 - 1;
node_total = [node1(1:end-1),node,node2(2:end)];

left_H = node(2:NV-1)-para.H(2:NV-1);                        % interior point
right_H = node(2:NV-1) + para.H(2:NV-1);
left_elem_num = whereis(left_H',node);                       % if left_H at the grid point x_k, then the elem_num = k-1
right_elem_num = whereis(right_H',node);



%% stiffnes matrix 
startNode = NV2+1; %% first interior point
endNode = NV2+NV-2; %% last interior point
freenode = startNode: endNode;


%% singular part
i1 = [];
j1 = [];
v1 = [];
A1 =sparse(i1,j1,v1,NT+1, 2*NT2+NT+1);
for k = 2:NV-1                                                             %% loop for nodes, can be implemented vectorizedly
    i1 = [k,k,k,k,k];
    j1 = [NT2+k, NT2+left_elem_num(k-1), NT2+left_elem_num(k-1)+1,...
                  NT2+right_elem_num(k-1), NT2+right_elem_num(k-1)+1];
    %% Compute the linear interpolation weight
    left_weight1 = left_H(k-1) - node_total(left_elem_num(k-1)+NT2);
    left_weight2 = node_total(left_elem_num(k-1)+1+NT2) - left_H(k-1);
    left_total = left_weight1+left_weight2;  
    right_weight1 = right_H(k-1) - node_total(right_elem_num(k-1)+NT2);
    right_weight2 = node_total(right_elem_num(k-1)+1+NT2) - right_H(k-1);
    right_total = right_weight1 + right_weight2;                           
    v1 = [2*para.H(k)^(-2*para.s),...
         -(left_weight2/left_total)*para.H(k)^(-2*para.s),-(left_weight1/left_total)*para.H(k)^(-2*para.s),...
         -(right_weight2/right_total)*para.H(k)^(-2*para.s),-(right_weight1/right_total)*para.H(k)^(-2*para.s)];
    tempA = sparse(i1,j1,v1, NT+1,2*NT2+NT+1);
    A1 = A1+tempA;
end
A1 = para.Cs*A1/(2-2*para.s);
clear i1 j1 v1


%% tail part 
% integral on the x_i - L, x_i + L, where L = 2;
% compute the left_h and right_h for x_i, x_i is interior node
i2 = [];
j2 = [];
v2 = [];
p = node_total;
A2 = sparse(i2,j2,v2,NT+1,2*NT2+NT+1);
for k = 2:NV-1
    % step 1: find the elem number xi-H, xi+H
    nl = left_elem_num(k-1);
    nr = right_elem_num(k-1);
    if isempty(find(node_total == left_H(k-1),1))
        nl = nl+1;
    else
        nl = nl+2;
    end
    % step 2: for unirorm mesh, only need compute the weight for y > x_i+H.
    %         for non-uniform, need to be modified. 
    wr = uniform_weight(node(k),nr, node_total,para,NV,1); %% note w only has 3NT-1 components
    wl = uniform_weight(node(k),nr, node_total,para,NV,0); 
    % x+H belongs to (node_total(NT+nr),node_total(NT+nr+1)]
    w2 = -(p(NT2+nr+1)-right_H(k-1))*DF(right_H(k-1)-node(k))/(p(NT2+nr+1)-p(NT2+nr))...
         +(F(p(NT2+nr+1)-node(k)) - F(right_H(k-1)-node(k)))/(p(NT2+nr+1)-p(NT2+nr)); % wk
    w3 = (F(p(NT2+nr+2)-node(k)) - F(p(NT2+nr+1)-node(k)))/(p(NT2+nr+2) - p(NT2+nr+1))...
         - (right_H(k-1) - p(NT2+nr))*DF(right_H(k-1)-node(k))/(p(NT2+nr+1) - p(NT2+nr)) ...
         - (F(p(NT2+nr+1)-node(k)) - F(right_H(k-1) - node(k)))/(p(NT2+nr+1)-p(NT2+nr));  % wk+1
     
    w4 = +(left_H(k-1)-p(NT2+nl-1))*DF(left_H(k-1)-node(k))/(p(NT2+nl)-p(NT2+nl-1))...
         -(F(left_H(k-1)-node(k)) - F(p(NT2+nl-1)-node(k)))/(p(NT2+nl)-p(NT2+nl-1)); % wk
     
    
    w5 = -(F(p(NT2+nl-2)-node(k)) - F(p(NT2+nl-1)-node(k)))/(p(NT2+nl-2) - p(NT2+nl-1))...
         + (left_H(k-1) - p(NT2+nl))*DF(left_H(k-1)-node(k))/(p(NT2+nl-1) - p(NT2+nl)) ...
         + (F(p(NT2+nl-1)-node(k)) - F(left_H(k-1) - node(k)))/(p(NT2+nl-1)-p(NT2+nl));  % wk+1

    n_right = find(node<right_H(k-1), 1, 'last' )-k+1;
    n_left = k-find(node>left_H(k-1), 1 )+1;  
    
    if n_right == 1
        wr(n_right) = w3;
        temp_v= [w2,wr];
    else
        wr(n_right-1) = w2;
        wr(n_right) = w3;
        wr(1:n_right-2) = 0;
        temp_v = [0,wr];
    end
    
    if n_left == 1
        wl(n_left) = w5;
        temp_u = [w4,wl];
    else
        wl(n_left-1) = w4;
        wl(n_left) = w5;
        wl(1:n_left-2) = 0;
        temp_u = [0,wl];
    end
    
    j2 = [(NT2+k):-1:k, NT2+k:2*NT2+k];
    i2 = [k*ones(1,NV2), k*ones(1,NV2)];
    v2 = [-temp_u, -temp_v];
    temA = sparse(i2,j2,v2,NT+1,2*NT2+NT+1);
    A2 = A2+temA;   
end
clear i2 j2 v2


%% truncated part
i = 2:NV-1;
j = freenode;
v = para.Cs*para.H(2:end-1).^(-2*para.s)/para.s;
A3 = sparse(i,j,v,NT2+1,2*NT2+NT+1);

A = A1 + A2 + A3;
% B1 = full(A1(2:end-1,freenode));
% B2 = full(A2(2:end-1,freenode));
A = A(2:end-1,freenode);

    
