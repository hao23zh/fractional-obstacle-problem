%% 1d obstacle problem
%% adapt mesh for the fractional Laplacian

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.N = 128;
para.h = 2/para.N;
para.rho = 1;                 %% boundary scale
para.s = 0.6;           
para.mu = para.s/(2-para.s);
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);
alpha = 1/2;

eg_num = 5;    
problem = eg_data(para,eg_num);
% u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 10^(-11);               %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 20;                      %% maximum iterating times 

%% mesh 
node = -1:para.h:1;
NV = size(node,2);  %% vertex number
freeNode = 2:NV-1;

para.hi = [node(2)-node(1),min(node(2:end-1)-node(1:end-2),...
           node(3:end)-node(2:end-1)),node(end)-node(end-1)]; %% compute the local hi
delta = dist2bd(node);     %% distance to the boundary 
H = (para.hi.^alpha).*min(delta.^(1-alpha), para.rho^(1-alpha)); %% scale for sigular part

% uI = u_exact(node);
node = node';
psih = psi(node);         
fh = f(node);                      
uh = max(0*node,psih);

%% Iteration
wn = 1; k = 0; % num1 = freeNode; num2 = [];
para.H = H;
A = fractional1Dmat(para,node');
while(wn>0 && k<nmax)
k = k+1;
% contact
e1 = A*uh(freeNode)-fh(freeNode);
e2 = uh(freeNode)-psih(freeNode);
wn = sum(abs(min(e1,e2))>error);    % number of points did not reach the abort condition
num1 = freeNode(e1<e2); % increase
num2 = setdiff(freeNode,num1); % contact
fprintf("iter:%d  notpass:%d  contact:%d \n",k,wn,size(num2,2));

% if wn == 0
%     break
% end

% Hi, A
para.H = modifyH1D(H,para.hi,node,num1,num2);
A = fractional1Dmat(para,node');

% uh
uh(num2) = psih(num2);
uh(num1) = A(num1-1,num1-1)\(fh(num1)-A(num1-1,num2-1)*uh(num2));

% draw
subplot(2,3,k)
plot(node,psih,'k','linewidth',2);
hold on 
plot(node,uh,'k')
plot(node(num1),uh(num1),'b.') % untouch
plot(node(num2),uh(num2),'r.') % touch
ylim([0,5])
title(['k=',int2str(k),'  |C_h^k|=',int2str(size(num2,2))])
hold off

end
% plot(node,psih,'k','linewidth',2);
% hold on 
% plot(node,uh,'k')
% plot(node(num1),uh(num1),'b.') % untouch
% plot(node(num2),uh(num2),'r.') % touch
% title(['iter ',int2str(k),'  notpass ',int2str(wn)])
% hold off
%fprintf("left:%f  right:%f  dist:%f \n",node(num2(1)),node(num2(end)),node(num2(end))-node(num2(1)));


%% end of file