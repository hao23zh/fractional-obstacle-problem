%% 1d obstacle problem
%% adapt mesh for the fractional Laplacian

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.N = 512;
para.h = 2/para.N;
para.rho = 1;                 %% boundary scale
para.s = 0.8;           
para.mu = para.s/(2-para.s);
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);

eg_num = 3;    
problem = eg_data(para,eg_num);
% u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 10^(-12);               %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 10;                      %% maximum iterating times 


%% mesh & stiffnes matrix 
node = (-1:para.h:1)';
A = fractional1Dmat(para);
NV = size(node,1);  %% vertex number
freeNode = 2:NV-1;

% uI = u_exact(node);
psih = psi(node);         
fh = f(node);                      
uh = 0*node;

%% Iteration
wn = 1; k = 0;
while(wn>0 && k<nmax)
k = k+1;    
e1 = A*uh(freeNode)-fh(freeNode);
e2 = uh(freeNode)-psih(freeNode);
wn = sum(abs(min(e1,e2))>error);    % number of points did not reach the abort condition
% disp(find(abs(min(e1,e2))>error))
fprintf("iter:%d   notpass :%d \n",k,wn);

num1 = freeNode(e1<e2);
num2 = setdiff(freeNode,num1);

uh(num2) = psih(num2);
uh(num1) = A(num1-1,num1-1)\(fh(num1)-A(num1-1,num2-1)*uh(num2));

%if wn ~= 0
subplot(3,3,k)
plot(node,psih,'k','linewidth',2);
hold on 
plot(node,uh,'k')
plot(node(num1),uh(num1),'b.') % untouch
plot(node(num2),uh(num2),'r.') % touch
title(['iter ',int2str(k),'  notpass ',int2str(wn)])
hold off
%end

end



%% end of file