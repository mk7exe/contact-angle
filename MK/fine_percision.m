%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies the fine precision droplet identification process
% to remove gas atoms which are very close to the liquid droplet surface
% and was not removed through hit-and-count process.
% 
% Parameters: 
% R_step: number of radial bins
% delta_R: size of the radial bin 
% (see the publication for details)
%
% Authors: 
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
function [x_final,y_final,z_final] = fine_percision(x_final,y_final,...
    z_final,R_step,delta_R)

dx_ave=max(x_final)-min(x_final);
dy_ave=max(y_final)-min(y_final);

D=max(dx_ave,dy_ave);
d=D/2;
Z_max=max(z_final);

%radious calculation
R=(d.^2+Z_max.^2)./(2*Z_max);
x_ave=min(x_final)+dx_ave/2;
y_ave=min(y_final)+dy_ave/2;
z_ave=Z_max-R;
%Find the distance from origin of the remaining points
% after his-and-count
l2=zeros(size(x_final));
for j=1:length(l2)
    l2(j)= sqrt((x_final(j)-x_ave).^2+...
        (y_final(j)-y_ave).^2+(z_final(j)-z_ave).^2);
end
%Define radious search region for evalutation of scattered points
R_interval=R-R_step*delta_R:delta_R:R+R_step*delta_R;
R_interval(end)=R_interval(end)+3;

% define three regions:smaller than search zone radious,
% search zone and bigger than search limit radious.

%First remove bigger than limit radious data
e=find(l2<R_interval(end));
x_final=x_final(e);
y_final=y_final(e);
z_final=z_final(e);
l2=l2(e);
%Second keep data less than minimum radious evaluation
% (do not touch)
e=find(l2<R_interval(1));
x_pass=x_final(e);
y_pass=y_final(e);
z_pass=z_final(e);

%Third evaluate data in the search zone for further action
e=find(l2>=R_interval(1) & l2<=R_interval(end));
x_zone=x_final(e);
y_zone=y_final(e);
z_zone=z_final(e);
l2=l2(e);

x_temp=0;
y_temp=0;
z_temp=0;

for j=1:length(R_interval)-1
    e=find(l2 >= R_interval(j) & l2 <= R_interval(j+1));
    x=x_zone(e);
    y=y_zone(e);
    z=z_zone(e);
    
    for m=1:length(x)
        count=0;
        for n=2:length(x_zone)
            distance=sqrt((x(m)-x_zone(n)).^2+ ...
                (y(m)-y_zone(n)).^2 ...
                +(z(m)-z_zone(n)).^2);
            if (distance<4.5 && distance>0)
                count=count+1;
            end
        end
        if count>=1
            x_temp=[x_temp;x(m)];
            y_temp=[y_temp;y(m)];
            z_temp=[z_temp;z(m)];
        end
    end
end

x_final=[x_pass;x_temp(2:end)];
y_final=[y_pass;y_temp(2:end)];
z_final=[z_pass;z_temp(2:end)];
end
