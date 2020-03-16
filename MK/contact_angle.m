%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the contact angle of a liquid droplet of a solid
% surface through convex hull triangulation.
% 
% Parameters: 
% z_min: Triangles with at least one vertex position lower that z_min*max(z)
% are contributing in the contact angle calculation.
% precision = This parameter is used to recognize triangles contributing
% to the area of the base of the droplet.
%
% Authors: 
% Mohammad Khalkhali
% Nasser Kazemi 
% (Sep 2016)
%
% References:
% Khalkhali et al. J. Chem. Phys. (2017)
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
function [K,V,angle,BaseA_temp] = contact_angle(x_final,y_final,z_final,...
    z_min,precision)

%apply convex hull algorithm to define the boundary

[K,V] = convhull(x_final,y_final,z_final,'simplify',true);

angle=zeros(2,length(K(:,1)));
BaseA_temp = 0;
%calculating angles
for jj=1:length(K(:,1))
    % applying z_limit
    if (any(z_final(K(jj,:)) < min(z_final) + z_min*max(z_final)))
        x=x_final(K(jj,:));
        y=y_final(K(jj,:));
        z=z_final(K(jj,:));
        
        p1=[x(1) y(1) z(1)];
        p2=[x(2) y(2) z(2)];
        p3=[x(3) y(3) z(3)];
        
        n = cross(p1-p2, p1-p3);%normal to the triangular facet
        angle(1,jj)=acos(n(3)./norm(n))*180/pi;
        angle(2,jj)=norm(n);
        % applying precision
        if ((abs(z(1)-z(2))<precision && ...
             abs(z(1)-z(3))<precision && ...
             abs(z(3)-z(2))<precision))
            
            BaseA_temp = BaseA_temp + norm(n); % area of the base
            angle(2,jj)=0;
        end
        
        if (all(z_final(K(jj,:)) < min(z_final) + z_min*max(z_final)))
            angle(2,jj)=0;
        end
        
    else
        angle(1,jj)=0;
        angle(2,jj)=0;
    end
end
end
