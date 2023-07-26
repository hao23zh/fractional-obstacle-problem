%% 1d obstacle problem
%% adapt mesh for the fractional Laplacian

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.rho = 1;                 %% boundary scale
para.s = 0.6;           
para.mu = para.s/(2-para.s);
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);

eg_num = 1;    
problem = eg_data(para,eg_num);
u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 5*10^(-11);               %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 600;                      %% maximum iterating times 

para.level = 5;
err_Linf = zeros(para.level,1);
num_mesh = [32,64,128,256,512]';%128]';
err_order = zeros(para.level,1);
num_dof = zeros(para.level,1);

for iter = 1:para.level
%% mesh & stiffnes matrix 
para.N = num_mesh(iter);
para.h = 2/para.N;
node = -1:para.h:1;
NV = size(node,2);  %% vertex number
freeNode = 2:NV-1;

para.hi = [node(2)-node(1),min(node(2:end-1)-node(1:end-2),...
           node(3:end)-node(2:end-1)),node(end)-node(end-1)]; %% compute the local hi
delta = dist2bd(node);     %% distance to the boundary 
H = (para.hi.^alpha).*min(delta.^(1-alpha), para.rho^(1-alpha)); %% scale for sigular part

node = node';
uI = u_exact(node);
psih = psi(node);         
fh = f(node);                      
uh = max(0*node,psih);


%% Iteration
wn = 1; k = -1;
para.H = H;
A = fractional1Dmat(para,node');
while(wn>0 && k<nmax)
k = k+1;    
e1 = A*uh(2:end-1)-fh(2:end-1);
e2 = uh(2:end-1)-psih(2:end-1);
wn = sum(abs(min(e1,e2))>error);    % number of points did not reach the abort condition
% disp(find(abs(min(e1,e2))>error))
% fprintf("iter:%d   notpass :%d \n",k,wn);

num1 = find(e1<e2)+1;
num2 = setdiff(2:NV-1,num1);
%fprintf("iter:%d  notpass:%d  contact:%d \n",k,wn,size(num2,2));

para.H = modifyH1D(H,para.hi,node,num1,num2);
A = fractional1Dmat(para,node');
uh(num2) = psih(num2);
uh(num1) = A(num1-1,num1-1)\(fh(num1)-A(num1-1,num2-1)*uh(num2));

end

%% error and update
max_err = max(abs(uh-uI));
err_Linf(iter) = max_err;
%num_mesh(iter) = 2*para.N;
fprintf("times:%d  Nmesh:%d  e:%f\n",k-1,NV-1,max_err);
num_dof(iter) = NV-1; 
% figure 
% plot(node,uII-uh);

if iter>=2
   err_order(iter) =  -log(err_Linf(iter-1)/err_Linf(iter))/log(num_dof(iter-1)/num_dof(iter)); 
end

end

%% Display error and time
   disp(['Example ',int2str(eg_num)])
   disp('Table: Error')
   colname = {'#Dof','||u-u_h||','Linf order'};
   disptable(colname,num_mesh,[],err_Linf,'%0.5e',err_order, '%0.2f');


%% end of file