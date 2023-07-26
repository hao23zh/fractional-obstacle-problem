function [ele_num,flag] = whereis(p,node)
%% whereis: find the number of the element which contains point p
%% Input: p, scalar or column vector,   position of the point p
%%        node,   row vector,        vertex of the mesh
%%
%% Output: ele_num, scalar or column vector, element number contains the point p
%%
%% Date: 7/19
flag = zeros(size(p));
ele_num = sum(node<p,2);
ele_num2= sum(node<=p,2);
flag(ele_num==ele_num2) = 1;
end