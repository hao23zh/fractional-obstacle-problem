function [err,uII] = error1(uh,node,uI,nodeI)
uII = uh;

% ��node��nodeI�е�λ��
flag = zeros(size(node));
ele_num = sum(nodeI<node');
ele_num2= sum(nodeI<=node');
flag(ele_num==ele_num2) = 1;

% uh��uI�غϵĶ���
uII(flag==0) = uI(ele_num(flag==0)+1);

% uh��uI��Ԫ�ڵĶ���
ln = ele_num(flag==1);rn = ln+1;
p = node(flag==1);
uII(flag==1) = (uI(ln).*(nodeI(rn)-p)+uI(rn).*(p-nodeI(ln)))./(nodeI(rn)-nodeI(ln));

err = uh-uII;
end