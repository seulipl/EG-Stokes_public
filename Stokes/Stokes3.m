function pde = Stokes3

ProbType = 1;

nu = 1e-4; k = 1;

pde = struct('exact_u',@ExaU,'exact_p',@ExaP,'rhs',@Rhs,'nu',nu,...
    'streamline',@streamline,'exact_gu',@GradU);

    function u = ExaU(node)
        x = node(:,1); y = node(:,2); z = node(:,3);
        switch ProbType
            case 1
                u(:,1) = sin(pi*x).*cos(pi*y) - sin(pi*x).*cos(pi*z);
                u(:,2) = sin(pi*y).*cos(pi*z) - sin(pi*y).*cos(pi*x);
                u(:,3) = sin(pi*z).*cos(pi*x) - sin(pi*z).*cos(pi*y);
            case 2
                u(:,1) = -y./(x.^2+y.^2+1);
                u(:,2) = x./(x.^2+y.^2+1);
                u(:,3) = z-z;
            case 3
                u(:,1) = sin(k*pi*x).*cos(k*pi*y) - sin(k*pi*x).*cos(k*pi*z);
                u(:,2) = sin(k*pi*y).*cos(k*pi*z) - sin(k*pi*y).*cos(k*pi*x);
                u(:,3) = sin(k*pi*z).*cos(k*pi*x) - sin(k*pi*z).*cos(k*pi*y);
        end
    end


    function Du = GradU(node)
        x = node(:,1); y = node(:,2); z = node(:,3);
        switch ProbType
            case 1
                Du(:,1,1) = pi*cos(pi*x).*cos(pi*y) - pi*cos(pi*x).*cos(pi*z);
                Du(:,2,1) = -pi*sin(pi*x).*sin(pi*y);
                Du(:,3,1) = pi*sin(pi*x).*sin(pi*z);
                Du(:,1,2) = pi*sin(pi*y).*sin(pi*x);
                Du(:,2,2) = pi*cos(pi*y).*cos(pi*z) - pi*cos(pi*y).*cos(pi*x);
                Du(:,3,2) = -pi*sin(pi*y).*sin(pi*z);
                Du(:,1,3) = -pi*sin(pi*z).*sin(pi*x);
                Du(:,2,3) = pi*sin(pi*z).*sin(pi*y);
                Du(:,3,3) = pi*cos(pi*z).*cos(pi*x) - pi*cos(pi*z).*cos(pi*y);
            case 2
                Du(:,1,1) = 2*x.*y./(x.^2+y.^2+1).^2;
                Du(:,2,1) = (-x.^2+y.^2-1)./(x.^2+y.^2+1).^2;
                Du(:,3,1) = z-z;
                Du(:,1,2) = (-x.^2+y.^2+1)./(x.^2+y.^2+1).^2;
                Du(:,2,2) = -2*x.*y./(x.^2+y.^2+1).^2;
                Du(:,3,2) = z-z;
                Du(:,1,3) = z-z;
                Du(:,2,3) = z-z;
                Du(:,3,3) = z-z;
            case 3
                Du(:,1,1) = k*(pi*cos(k*pi*x).*cos(k*pi*y) - pi*cos(k*pi*x).*cos(k*pi*z));
                Du(:,2,1) = -k*pi*sin(k*pi*x).*sin(k*pi*y);
                Du(:,3,1) = k*pi*sin(k*pi*x).*sin(k*pi*z);
                Du(:,1,2) = k*pi*sin(k*pi*y).*sin(k*pi*x);
                Du(:,2,2) = k*(pi*cos(k*pi*y).*cos(k*pi*z) - pi*cos(k*pi*y).*cos(k*pi*x));
                Du(:,3,2) = -k*pi*sin(k*pi*y).*sin(k*pi*z);
                Du(:,1,3) = -k*pi*sin(k*pi*z).*sin(k*pi*x);
                Du(:,2,3) = k*pi*sin(k*pi*z).*sin(k*pi*y);
                Du(:,3,3) = k*(pi*cos(k*pi*z).*cos(k*pi*x) - pi*cos(k*pi*z).*cos(k*pi*y));
        end
    end

    function p = ExaP(node)
        x = node(:,1); y = node(:,2); z = node(:,3);
        switch ProbType
            case 1
                p = pi*sin(pi*x).*sin(pi*y).*sin(pi*z)-1/pi^2;
            case 2
                p = (x>=0.5).*(2*x-1) + (x<0.5).*(-2*x+1);
            case 3
                p = sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
        end
    end


    function f = Rhs(node)
       u = ExaU(node);
        x = node(:,1); y = node(:,2); z = node(:,3);
        u1 = u(:,1); u2 = u(:,2);
        switch ProbType
            case 1
                f(:,1) = pi*(nu*2*pi*sin(pi*x).*cos(pi*y) - nu*2*pi*sin(pi*x).*cos(pi*z))...
                    + pi*pi*cos(pi*x).*sin(pi*y).*sin(pi*z);
                f(:,2) = pi*(nu*2*pi*sin(pi*y).*cos(pi*z) - nu*2*pi*sin(pi*y).*cos(pi*x))...
                    + pi*pi*sin(pi*x).*cos(pi*y).*sin(pi*z);
                f(:,3) = pi*(nu*2*pi*sin(pi*z).*cos(pi*x) - nu*2*pi*sin(pi*z).*cos(pi*y))...
                    + pi*pi*sin(pi*x).*sin(pi*y).*cos(pi*z);
            case 2
                f(:,1) = -8*nu*y./(x.^2+y.^2+1).^3 + (x>=0.5).*2 + (x<0.5).*(-2);
                f(:,2) = 8*nu*x./(x.^2+y.^2+1).^3;
                f(:,3) = z-z;
            case 3
                f(:,1) = k*pi*(nu*2*k*pi*sin(k*pi*x).*cos(k*pi*y) - nu*2*k*pi*sin(k*pi*x).*cos(k*pi*z))...
                    + k*pi*cos(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
                f(:,2) = k*pi*(nu*2*k*pi*sin(k*pi*y).*cos(k*pi*z) - nu*2*k*pi*sin(k*pi*y).*cos(k*pi*x))...
                    + k*pi*sin(k*pi*x).*cos(k*pi*y).*sin(k*pi*z);
                f(:,3) = k*pi*(nu*2*k*pi*sin(k*pi*z).*cos(k*pi*x) - nu*2*k*pi*sin(k*pi*z).*cos(k*pi*y))...
                    + k*pi*sin(k*pi*x).*sin(k*pi*y).*cos(k*pi*z);
        end


    end

    function u = streamline(x,y)
        u=5*x.^2.*(x-1).^2.*y.^2.*(y-1).^2;
    end

end

