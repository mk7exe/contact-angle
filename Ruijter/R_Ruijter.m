%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the contact angle of a liquid droplet of a solid
% surface through a circular fitting to the boundary points. Boundary of
% the liquid droplet at different heights from the solid surface is
% calculated through fitting the radial density profile to a sigmoidal
% function (see the reference for details)
% 
% Parameters: 
% z_min: minimun z of data points
% z_max: maximum z of data points
% dz and dr: size of the radial density profile
%
% Author: Mohammad Khalkhali, Sep 2016
% Reference: 
% Khalkhali et al. J. Chem. Phys. (2017)
% M. J. de Ruijter, T. D. Blake, and J. De Coninck, Langmuir 15,7836(1999).
% 
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
function [theta,x,y] = R_Ruijter(x_final,y_final,z_final,zmin,zmax,dz,dr)

PI = 3.14159265359;
R_temp = zeros(fix((zmax-zmin)/dz),1);
ro_l = 0.033;
ro_v = 0.0;
ii = 0;
jj = 0;

Rmax = 0.5*sqrt((max(x_final)-min(x_final))^2+(max(y_final)-min(y_final))^2);

opts = optimset('Display','off');

for i=fix(zmin/dz)*dz:dz:fix(zmax/dz)*dz
    ii = ii + 1;
    j = find(z_final>= i & z_final < i + dz);
    if (size(j,1) ~= 0)
        xtemp = x_final(j);
        ytemp = y_final(j);
        rmax = 0.5*sqrt((max(xtemp)-min(xtemp))^2+(max(ytemp)-min(ytemp))^2);
        ro = zeros(fix(rmax/dr),1);
        x_c = mean(xtemp);
        y_c = mean(ytemp);
        xtemp = xtemp - x_c;
        ytemp = ytemp - y_c;
        
        %         if (ii == 10)
        %             figure;
        %             scatter(xtemp,ytemp);
        %             str = sprintf('z=%d',i);
        %             title(str);
        %         end
        
        r_temp = sqrt(xtemp.^2+ytemp.^2);
        
        iii = 0;
        for k=0:dr:fix(rmax/dr)*dr
            iii = iii + 1;
            j = find(r_temp >= k & r_temp < k+dr);
            ro(iii) = size(j,1)/(PI*(2*k*dr+dr^2)*dz);
        end
        
        j = find(ro ~= 0);
        
        if (size(j,1) > 1)
            F = @(d,xdata) 0.5*(ro_l+ro_v)-0.5*(ro_l-ro_v)*...
                tanh((xdata-d(1))/d(2));
            d0 = [0,1];
            r = 0:dr:fix(rmax/dr)*dr;
            r = r';
            d = lsqcurvefit(F,d0,r,ro,[0,0],[rmax,5],opts);
            R_temp(ii) = d(1);
            %             if (ii == 10)
            %                 figure;
            %                 scatter(r,ro);
            %                 hold on
            %                 plot(r,F(d,r),'r');
            %
            %             end
        else
            R_temp(ii) = 0;
        end
    else
        R_temp(ii) = 0.0;
    end
    
end

% Fit to circle
F = @(R,xdata) sqrt(R(1)^2-(xdata).^2)-R(2);
R0 = [0,0];
z = fix(zmin/dz)*dz:dz:fix(zmax/dz)*dz;
z = z';
j = find(R_temp > 1);
x = R_temp(j);
y = z(j);
[R,resnorm] = lsqcurvefit(F,R0,x,y,[Rmax/2,-5],[Rmax,max(y)],opts);

if (resnorm < 100)
    x0 = sqrt(R(1)^2-R(2)^2);
    diff = -x0/sqrt(R(1)^2-x0^2);
    theta = abs(atan(diff)*180/PI);
    
%         figure;
%         scatter(x,y);
%         hold on;
%         r1 = 0:0.001:x0;
%         z1 = sqrt(R(1)^2-(r1).^2)-R(2);
%         plot(r1,z1,'r-');
%         % axis([0 max(r1) 0 max(z)]);
else
    theta = 0.0;
end

x = x';
y = y';

end
