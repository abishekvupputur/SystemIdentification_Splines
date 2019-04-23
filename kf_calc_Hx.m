%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = kf_calcDHx(x) Calculates the Jacobian of the output dynamics equation f(x,u,t) 
%   
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = kf_calcDHx(t, x, u)

%     zpred(1) = atan ( x(3) / x(1) ) * (1 + x(4));
%     zpred(2) = atan ( x(2) / sqrt(x(1)^2 + x(3)^2));
%     zpred(3) = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    x1=x(1);
    x2=x(2);
    x3=x(3);
    x4=x(4);
    Hx(1,1) = -(x3*(x4 + 1))/(x1^2*(x3^2/x1^2 + 1));
    Hx(1,2) = 0;
    Hx(1,3) = (x4 + 1)/(x1*(x3^2/x1^2 + 1));
    Hx(1,4) = atan(x3/x1);
    Hx(2,1) = -(x1*x2)/((x2^2/(x1^2 + x3^2) + 1)*(x1^2 + x3^2)^(3/2));
    Hx(2,2) = 1/((x2^2/(x1^2 + x3^2) + 1)*(x1^2 + x3^2)^(1/2));
    Hx(2,3) = -(x2*x3)/((x2^2/(x1^2 + x3^2) + 1)*(x1^2 + x3^2)^(3/2));
    Hx(2,4) = 0;
    Hx(3,1) = x1/(x1^2 + x2^2 + x3^2)^(1/2);
    Hx(3,2) = x2/(x1^2 + x2^2 + x3^2)^(1/2);
    Hx(3,3) = x3/(x1^2 + x2^2 + x3^2)^(1/2);
    Hx(3,4) = 0;
end
