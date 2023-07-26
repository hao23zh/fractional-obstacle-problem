function problem = eg_data(para,num)
% Data for obstacle

problem = struct('psi',@psi,...
             'f',@f,...
             'u_exact',@u_exact);

         
    function psi = psi(x)
        if num == 1
            psi = para.Ks*(1-x.^2).^(para.s)-(abs(x)>1/3);
        elseif num == 2
            psi = para.Ks*(1-x.^2).^(para.s)-(abs(x)<1/6)/4-(1/2>abs(x))/4-(abs(x)>3/4)/4;
        elseif num == 3
            psi = max(0*x,1-4*abs(x-1/4));
        end
    end


    function f = f(x)
        if num == 1
            f = 0*x+1;
        elseif num == 2
            f = 0*x;
        elseif num == 3
            f = 0*x+1/2;
        end
    end


    function u = u_exact(x)
        % Lu = 1
        % u may not be the exact solution to the obstacle problem 123
        u = para.Ks*(1-x.^2).^(para.s);
    end

end