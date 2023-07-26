function [node,hi,mesh] = generate_adaptmesh(para)
%% generate_adaptmesh: generate adaptive mesh for fractional laplacian
%% To improve the order s to order 1
%% Input: para, structure, parameter
%% Outpu: node, row vector, node coordinate
%% Latest Version: 7/24/2021
%% Coder: Rubing Han
N = 1/para.h;

size =1;
a(1) = 1;
sum(1) = 0;
if sum(end)*para.h^(1/para.mu) < 0.1
    a =[a,zeros(1,N)];
    sum = [sum,zeros(1,N)];
    for k = size+1:size+N
        sum(k) = sum(k-1)+a(k-1);
        a(k) = sum(k)^(1-para.mu);
    end
    size = size+N;
end
sum = sum*para.h^(1/para.mu);
h_left = a*para.h^(1/para.mu);
h_right = a*para.h^(1/para.mu);
h_right = fliplr(h_right);

node_left = -1+sum;
node_right = fliplr(1-sum);
node_inter = node_left(end):para.h:node_right(1);
h_inter = node_inter(2:end) - node_inter(1:end-1);

% keyboard;
if node_inter(end)-node_inter(end-1) < 0.2*para.h
    node_inter(end-1) = node_inter(end);
    node_inter = node_inter(1:end-1);
end
node = [node_left(1:end-1),node_inter,node_right(2:end)];
hi = [h_left(1:end-1),h_inter,h_right(2:end)];
mesh = hi;
hi = [hi(1),min(hi(1:end-1),hi(2:end)),hi(end)];

end