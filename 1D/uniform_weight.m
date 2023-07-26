function [w] = uniform_weight(t,nr,node_total,para,NV,flag)
%% uniform_weight: generate the tail weight for uniform mesh
%% Input: p,    scalar,      coordinate of point
%%        node, vector,      computaion domain
%%
%%        para, structure,   parameter
%%        NV,   scalar,      degree of freedom
%% Output:w,    vector,      weight for node.
%% primitive function F
if para.s == 0.5
    F = @(x) -para.Cs*log(abs(x)+eps);
    DF =@(x) -sign(x).*para.Cs./(abs(x) + eps);
else
    F = @(x) para.Cs*(abs(x)+eps).^(1-2*para.s)/2/para.s/(2*para.s - 1);
    DF =@(x) -para.Cs*sign(x).*((abs(x)+eps).^(-2*para.s))/para.s/2;
end
G = @(x) para.Cs./(abs(x)+eps).^(1+2*para.s);
% node_total = [node_total,node_total(end)+para.h];
    
if flag == 1
    w = zeros(1,NV-1);
%     w2 = zeros(1,NV-1);
    p1 = node_total - t;
    p1 = p1(p1>0);
    p1 = p1(1:NV-1);
%     [lambda,weight] = quadpts1(para.order);
%     m = size(lambda,1);
%     for i = 1:m
%         point1 = lambda(i,1)*p1(1:end-1) + lambda(i,2)*p1(2:end);
%         w2(2:end) = w2(2:end) + weight(i)*G(point1).*(point1-p1(1:end-1));
%         point2 = lambda(i,1)*p1(2:end-1) + lambda(i,2)*p1(3:end);
%        w2(2:end-1) = w2(2:end-1) + weight(i)*G(point2).*(p1(3:end)-point2);
%     end
    w(2:end-1) = (F(p1(3:end)) - F(p1(2:end-1)))./(p1(3:end) - p1(2:end-1)) -...
                 (F(p1(2:end-1)) - F(p1(1:end-2)))./(p1(2:end-1) - p1(1:end-2));
    w(NV-1) = DF(p1(end)) - (F(p1(end)) - F(p1(end-1)))/(p1(end) - p1(end-1));
else
    w = zeros(1,NV-1);
    p1 = node_total - t;
    p1 = p1(p1<0);
    p1 = p1(end-NV+2:end);
%     [lambda,weight] = quadpts1(para.order);
%     m = size(lambda,1);
%     keyboard;
%     for i = 1:m
%         point1 = lambda(i,1)*p1(1:end-1) + lambda(i,2)*p1(2:end);
%         w2(2:end) = w2(2:end) + weight(i)*G(point1).*(point1-p1(1:end-1));
%         point2 = lambda(i,1)*p1(1:end-2) + lambda(i,2)*p1(2:end-1);
%        w2(2:end-1) = w2(2:end-1) + weight(i)*G(point2).*(p1(2:end-1)-point2);
%     end
    w(2:end-1) = (F(p1(3:end)) - F(p1(2:end-1)))./(p1(3:end) - p1(2:end-1)) -...
                 (F(p1(2:end-1)) - F(p1(1:end-2)))./(p1(2:end-1) - p1(1:end-2));
    w(1) = -DF(p1(1)) + (F(p1(1)) - F(p1(2)))/(p1(1) - p1(2));
    w = fliplr(w); 
end





end