%% 1d obstacle problem
%% adapt mesh for the fractional Laplacian

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
para.rho = 1;                 %% boundary scale
para.s = 0.8;           
para.mu = para.s/(2-para.s);
para.Cs = (para.s*2^(2*para.s)*gamma(para.s + 0.5))/(sqrt(pi)*gamma(1-para.s));
para.Ks = 2^(-2*para.s)*gamma(0.5)/gamma(1+para.s)/gamma(0.5+para.s);

eg_num = 3;    
problem = eg_data(para,eg_num);
u_exact = problem.u_exact;
psi = problem.psi;              %% obstacle function
f = problem.f;

error = 10^(-12);               %% termination condition  min(uh-psih,Lu-fh) < error
nmax = 600;                      %% maximum iterating times 


NI = 512;
para.level = 5;
err_Linf = zeros(para.level,1);
num_mesh = [NI,16,32,64,128]';
err_order = zeros(para.level,1);


for iter = 1:para.level
%% mesh & stiffnes matrix 
para.N = num_mesh(iter);
[A,node] = fractional1Dmat_adap(para);
node = node';
NV = size(node,1);  %% vertex number

% uI = u_exact(node);
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

%%
if iter == 1
    uI = uh;
    nodeI= node;
end

%% error and update
if iter >= 2
[err,uII] = error1(uh,node,uI,nodeI);
max_err = max(abs(err));
err_Linf(iter) = max_err;
fprintf("times:%d  Nmesh:%d  e:%f\n",k-1,NV-1,max_err);
% figure 
% plot(node,uII-uh);

if iter>=3
   err_order(iter) =  log(err_Linf(iter-1)/err_Linf(iter))/log(2); 
end

end
end

%% Display error and time
   disp(['Example ',int2str(eg_num)])
   disp('Table: Error')
   colname = {'#Dof','||u-u_h||','Linf order'};
   disptable(colname,num_mesh,[],err_Linf,'%0.5e',err_order, '%0.2f');


%% end of file