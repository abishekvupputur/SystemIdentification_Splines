%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%   
%   Author: C.C. de Visser, Delft University of Technology, 2013
%   email: c.c.devisser@tudelft.nl
%   Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpred = kf_calcHx(t, x, u)
    
    zpred(1) = atan ( x(3) / x(1) ) * (1 + x(4));
    zpred(2) = atan ( x(2) / sqrt(x(1)^2 + x(3)^2));
    zpred(3) = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    end
    