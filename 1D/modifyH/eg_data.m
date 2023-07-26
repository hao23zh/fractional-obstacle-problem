function problem = eg_data(para,num)
% Data for obstacle

problem = struct('psi',@psi,...
             'f',@f,...
             'u_exact',@u_exact);

    s = para.s;
    function psi = psi(x)
        if num == 1
        %    y = 1/2-abs(x);
        %    psi = para.Ks*(1-x.^2).^s.*(y>=0)+...
        %        para.Ks*((0.75)^s+s*(0.75)^(s-1)*y+s*(0.75)^(s-2)*(2*s-5)/4*y.^2).*(y<0);
            psi = para.Ks*(1-x.^2).^(para.s)+0.5*min(1/4-x.^2,0);
        elseif num == 2
            psi = para.Ks*(1-x.^2).^(para.s)-(abs(x)<1/6)/4-(1/2>abs(x))/4-(abs(x)>3/4)/4;
        elseif num == 3
            psi = max(0*x,1-4*abs(x-1/4));
        elseif num == 4
            psi = max(1-3*(x-1/4).^2,-0.1);%min(max(0*x,sin(2*pi*x)),2*(x+0.5));
        elseif num == 5
            psi = 3*(1-2*abs(x-1/4));
        end
    end


    function f = f(x)
        if num == 1
            f = 1-5*max(1/2-abs(x),0);
        elseif num == 2
            f = 0*x;
        elseif num == 3
            f = 0*x;
        elseif num == 4
            f = 0*x;
        elseif num == 5
            f = 0*x+1;
        end
    end


    function u = u_exact(x)
        % Lu = 1
        % u may not be the exact solution to the obstacle problem 123
        u = para.Ks*(1-x.^2).^(para.s);
    end

end