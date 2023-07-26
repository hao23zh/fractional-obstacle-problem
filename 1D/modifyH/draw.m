%% 1d obstacle problem
%% adapt mesh for the fractional Laplacian

%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.rho = 1;                 %% boundary scale
para.s = 0.6;           
para.mu = para.s/(2-para.s);
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);
alpha = 1/2;

eg_num = 5;    
problem = eg_data(para,eg_num);
%u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 10^(-11);               %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 600;                      %% maximum iterating times 


%% mesh & stiffnes matrix 
para.N = 256;
para.h = 1/para.N;
[node,para.hi,~] = generate_adaptmesh(para);

delta = dist2bd(node);     %% distance to the boundary 
H = (para.hi.^alpha).*min(delta.^(1-alpha), para.rho^(1-alpha)); %% scale for sigular part

node = node';
NV = size(node,1);  %% vertex number

uI = u_exact(node);
psih = psi(node);         
fh = f(node);                      
uh = max(psih,0);
rerr = [];

%% Iteration
wn = 1; k = 0;
para.H = H;
A = fractional1Dmat_adap(para,node');
while(wn>0 && k<nmax)
k = k+1;   
% contact
e1 = A*uh(2:end-1)-fh(2:end-1);
e2 = uh(2:end-1)-psih(2:end-1);
rerr(k) = max(abs(min(e1,e2)));
wn = sum(abs(min(e1,e2))>error);    % number of points did not reach the abort condition
num1 = find(e1<e2)+1;
num2 = setdiff(2:NV-1,num1);
fprintf("iter:%d  notpass:%d  contact:%d \n",k,wn,size(num2,2));

% Hi A
para.H = modifyH1D(H,para.hi,node,num1,num2);
A = fractional1Dmat_adap(para,node');

% uh
uh(num2) = psih(num2);
uh(num1) = A(num1-1,num1-1)\(fh(num1)-A(num1-1,num2-1)*uh(num2));

% draw
subplot(3,3,k)
plot(node,psih,'k','linewidth',2);
hold on 
plot(node,uh,'k')
plot(node(num1),uh(num1),'b.') % untouch
plot(node(num2),uh(num2),'r.') % touch
ylim([0,3])
title(['k=',int2str(k),'  |C_h^{(k)}|=',int2str(size(num2,2))])
hold off
end

e1 = A*uh(2:end-1)-fh(2:end-1);
e2 = uh(2:end-1)-psih(2:end-1);
rerr(k+1) = max(abs(min(e1,e2)));
%% 
subplot(3,3,9)
plot(1:k,log10(rerr(2:end)),'o-','linewidth',1.5);
xlabel('iteration');
ylabel('log(||G_h[u_h^{(k)}]||_{L_\infty})');
title('residual')

% plot(node,psih,'k','linewidth',2);
% hold on 
% plot(node,uI)
% plot(node,uh,'k')
% plot(node(num1),uh(num1),'b.') % untouch
% plot(node(num2),uh(num2),'r.') % touch
% title(['iter ',int2str(k),'  notpass ',int2str(wn)])
% hold off

%% end of file