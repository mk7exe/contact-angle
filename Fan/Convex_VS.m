%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function calculates the volume and surface area of a data set
% through convex hull triangulation
%
% Author: Mohammad Khalkhali, Sep 2016
% Reference: Khalkhali et al. J. Chem. Phys. (2017)
% 
% This program is a free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details: http://www.gnu.org/licenses/ 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,S] = VS_Fan(x_final,y_final,z_final);

%apply convex hull algorithm to define the boundary

[K,V] = convhull(x_final,y_final,z_final,'simplify',true);

S = sum(sqrt(sum(( ...
    [y_final(K(:,1)).*z_final(K(:,2)) - z_final(K(:,1)).*y_final(K(:,2)) ...
    z_final(K(:,1)).*x_final(K(:,2)) - x_final(K(:,1)).*z_final(K(:,2))  ...
    x_final(K(:,1)).*y_final(K(:,2)) - y_final(K(:,1)).*x_final(K(:,2))] + ...
    [y_final(K(:,2)).*z_final(K(:,3)) - z_final(K(:,2)).*y_final(K(:,3)) ...
    z_final(K(:,2)).*x_final(K(:,3)) - x_final(K(:,2)).*z_final(K(:,3))  ...
    x_final(K(:,2)).*y_final(K(:,3)) - y_final(K(:,2)).*x_final(K(:,3))] + ...
    [y_final(K(:,3)).*z_final(K(:,1)) - z_final(K(:,3)).*y_final(K(:,1)) ...
    z_final(K(:,3)).*x_final(K(:,1)) - x_final(K(:,3)).*z_final(K(:,1))  ...
    x_final(K(:,3)).*y_final(K(:,1)) - y_final(K(:,3)).*x_final(K(:,1))]) ...
    .^2,2)))/2;

end
