function [err,uII] = error1(uh,node,uI,nodeI)
uII = uh;

% 找node在nodeI中的位置
flag = zeros(size(node));
ele_num = sum(nodeI<node');
ele_num2= sum(nodeI<=node');
flag(ele_num==ele_num2) = 1;

% uh和uI重合的顶点
uII(flag==0) = uI(ele_num(flag==0)+1);

% uh在uI单元内的顶点
ln = ele_num(flag==1);rn = ln+1;
p = node(flag==1);
uII(flag==1) = (uI(ln).*(nodeI(rn)-p)+uI(rn).*(p-nodeI(ln)))./(nodeI(rn)-nodeI(ln));

err = uh-uII;
end