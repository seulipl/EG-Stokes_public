function pde = Stokes2

ProbType = 1;
nu = 1e+0; lambda = 1/(2*nu) - sqrt(1/(4*nu^2)+4*pi^2); Ra = 1000;

pde = struct('exact_u',@ExaU,'exact_p',@ExaP,'rhs',@Rhs,'nu',nu,...
    'streamline',@streamline,'exact_gu',@GradU);

    function u = ExaU(node)
        x = node(:,1); y = node(:,2);
        switch ProbType
            case 1
                u(:,1) = 10*x.^2.*y.*(x-1).^2.*(2*y-1).*(y-1);
                u(:,2) = -10*x.*y.^2.*(2*x-1).*(x-1).*(y-1).^2;
            case 2
                u(:,1) = x-x;
                u(:,2) = y-y;
            case 3
                u(:,1) = 1-exp(lambda*x).*cos(2*pi*y);
                u(:,2) = lambda*exp(lambda*x).*sin(2*pi*y)/(2*pi);
            case 4
                u(:,1) = sin(pi*x).*sin(pi*y);
                u(:,2) = cos(pi*x).*cos(pi*y);
        end
    end


    function Du = GradU(node)
        x = node(:,1); y = node(:,2);
        switch ProbType
            case 1
                Du(:,1,1) = 20*(x-1).*x.*(2*x-1).*(y-1).*y.*(2*y-1);
                Du(:,2,1) = 10*(x-1).^2.*x.^2.*(6*y.^2-6*y+1);
                Du(:,1,2) = -10*(6*x.^2-6*x+1).*(y-1).^2.*y.^2;
                Du(:,2,2) = 20*(x-1).*x.*(2*x-1).*(1-2*y).*(y-1).*y;
            case 2
                Du(:,1,1) = x-x;
                Du(:,2,1) = x-x;
                Du(:,1,2) = y-y;
                Du(:,2,2) = y-y;
            case 4
                Du(:,1,1) = pi*cos(pi*x).*sin(pi*y);
                Du(:,2,1) = pi*sin(pi*x).*cos(pi*y);
                Du(:,1,2) = -pi*sin(pi*x).*cos(pi*y);
                Du(:,2,2) = -pi*cos(pi*x).*sin(pi*y);
        end
    end

    function p = ExaP(node)
        x = node(:,1); y = node(:,2);
        switch ProbType
            case 1
                p = 10*(2*x-1).*(2*y-1);
            case 2
                p = -Ra*y.^2/2 + Ra*y - Ra/3;
            case 3
                p = -exp(2*lambda*x)/2 + exp(x);
            case 4
                p = (y>=0).*y + (y<0).*(-y);
        end
    end


    function f = Rhs(node)
       u = ExaU(node);
        x = node(:,1); y = node(:,2);
        u1 = u(:,1); u2 = u(:,2);
        switch ProbType
            case 1
                f(:,1) = nu*(-120*(y-1/2).*(x.^4-2*x.^3+(2*y.^2-2*y+1) ...
                    .*x.^2+(-2*y.^2+2*y).*x+y.^2/3-y/3)) + 20*(2*y-1);% + (2*x-1).*(y.^2-y); % 
                f(:,2) = nu*(240*(x-1/2).*((y.^2-y+1/6).*x.^2+ ...
                    (-y.^2+y-1/6).*x+y.^2.*(y-1).^2/2)) + 20*(2*x-1);% + (x.^2-x).*(2*y-1); %
            case 2
                f(:,1) = x-x;
                f(:,2) = -Ra*y + Ra;

            case 3
                f(:,1) = lambda*(1-u1)-lambda*exp(2*lambda*x) + exp(x);
                f(:,2) = -lambda*u2;
            case 4
                f(:,1) = nu*2*pi^2*u1...
                    +0;
                f(:,2) = nu*2*pi^2*u2...
                    +(y>=0).*1 + (y<0).*(-1);
        end
    end

    function u = streamline(x,y)
        u=5*x.^2.*(x-1).^2.*y.^2.*(y-1).^2;
    end

end

