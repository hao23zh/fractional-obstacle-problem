%% 1d obstacle problem

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.h = 1/16;
para.rho = 1;                 %% boundary scale
para.s = 0.8;           
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);

eg_num = 3;    
problem = eg_data(para,eg_num);
u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 10^(-12);                       %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 10;                             %% maximum iterating times 


para.level = 5;
err_Linf = zeros(para.level,1);
num_mesh = zeros(para.level,1);
err_order = zeros(para.level,1);


for iter = 1:para.level
%% generate mesh 
node = (-1: para.h : 1)';               %% unifrom mesh
NV =  size(node,1);                     %% vertex number
A = fractional1Dmat(para);

uI = u_exact(node);
psih = psi(node);         
fh = f(node);                      
uh = 0*node;


%% Iteration
wn = 1; k = 0;
while(wn>0 && k<nmax)
k = k+1;    
e1 = A*uh(2:end-1)-fh(2:end-1);
e2 = uh(2:end-1)-psih(2:end-1);
wn = sum(abs(min(e1,e2))>error);    % number of points did not reach the abort condition
% disp(find(abs(min(e1,e2))>error))
% fprintf("iter:%d   notpass :%d \n",k,wn);

num1 = find(e1<e2)+1;
num2 = setdiff(2:NV-1,num1);

uh(num2) = psih(num2);
uh(num1) = A(num1-1,num1-1)\(fh(num1)-A(num1-1,num2-1)*uh(num2));

end


%% error and update
max_err = max(abs(uh-uI));
err_Linf(iter) = max_err;
num_mesh(iter) = NV-1;
% fprintf("times:%d  Nmesh:%d  e:%f\n",k-1,NV-1,max_err);
% figure 
% plot(node,u_exact-uh);

if iter>=2
   err_order(iter) =  log(err_Linf(iter-1)/err_Linf(iter))/log(2); 
end

para.h = para.h/2;

end


%% Display error and time
    disp(['Example ',int2str(eg_num)])
    disp('Table: Error')
    colname = {'#Dof','||u-u_h||','Linf order'};
    disptable(colname,num_mesh,[],err_Linf,'%0.5e',err_order, '%0.2f');


%% end of file


