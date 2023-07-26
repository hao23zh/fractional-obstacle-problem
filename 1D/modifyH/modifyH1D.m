function H1 = modifyH1D(H,h,node,num1,num2)
% num2 is contact: uh(num2) = psih(num2);
% h <= H1 <= H 
H1 = H;
Nc = size(num2,2);
Nnc = size(num1,2);

if Nc > 0
%    H1(num2) = h(num2).^0.5;
    nct = node(num1)';
    ct = node(num2);
    dist = abs(repmat(nct,Nc,1) - repmat(ct,1,Nnc));
    distnc = min(dist);
    H1(num1) = min(H1(num1),distnc);
    %
    distc = min(dist,[],2);
    H1(num2) = min(H1(num2),distc');
end

end