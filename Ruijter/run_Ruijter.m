%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs the calculation of contact angle of a liquid droplet
% on a solid surface using the method of de Ruijter et al. (1991). The
% boundary points of the liquid droplet at each MD snapshot are found and
% fitted to the circular equation in function R_Ruijter. The contact angle
% is then calculated from the tangent line to the fitted circle.
%
% Author: Mohammad Khalkhali, Sep 2016
% References:
% Khalkhali et al. J. Chem. Phys. (2017)
% M. J. de Ruijter, T. D. Blake, and J. De Coninck,Langmuir 15,7836(1999).
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
fclose ('all');
close('all');

addpath('../MK');
trajName = '../Graphite_Water(3nm)/Graphite_Water.lammpstrj';
% This format_spec is needed when reading lammps trajectories in
% read_LAMMPS_traj function. There are 5 entries at each line of our
% trajectories
format_spec = '%f %f %f %f %f';
StartStep = 4990000;
EndStep = 5000000;
step = 1000;

% Parameters of hit-and-count function (refer to the publication for details)
binsize=2;
numer_atom_in_bin=5;

% Parameters of fine_percision function (refer to the publication for details)
% fine presicion step of identifying droplet limits is usually redundant
% (hit and count step is sufficient).
fine_percision_check = 0; % calls fine_percision function if is 1
delta_R = 5;
R_step = 1;

% Draws graphs corresponding to each step for the first time step. It is
% recommended to check if parameters droplet identification process are
% ajusted properly.
graphcheck = 0;

% Parameters of the radial density profile (unit is angstrom)
dz = 2;
dr = 1;

ii = 0;
theta = zeros((EndStep-StartStep)/step+1,1);
x_ave = [];
y_ave = [];
        
for tt = StartStep:step:EndStep
    
    ii = ii + 1;
    fprintf('Current time step = %d\r', tt);
    Filename = sprintf('%s.%d',trajName,tt); %trajectory file name
    
    % Reading box dimentions to calculates centre of mass of water
    % molecules correctly: some water molecules may pass the boundary.
    % Since it is a NVT simulation we just need to read box dimentions
    % once.
    
    if (tt == StartStep)
        fileID = fopen(Filename,'r');
        if (fileID < 0)
            fprintf('can not open %s file',Filename);
            break;
        end
        %reading box sizes
        header = textscan(fileID, '%s',9,'delimiter', '\n');
        data = textscan(header{1}{6},'%f %f');
        Box_x = data{2}-data{1};
        data = textscan(header{1}{7},'%f %f');
        Box_y = data{2}-data{1};
        data = textscan(header{1}{8},'%f %f');
        Box_z = data{2}-data{1};
        fclose(fileID);
    end
    
    [fileID,x_org,y_org,z_org] = read_LAMMPS_traj(Filename,Box_x,...
        Box_y,Box_z,format_spec);
    
    if (fileID < 0)
        fprintf('can not open %s file',Filename);
        break;
    end
    
    % Drawing the droplet before applying hit-and-run
    if (tt == StartStep && graphcheck == 1)
        figure;
        scatter3(x_org,y_org,z_org,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .25 .25]);
        view(45,45)
        title('original');
        size1 = size(x_org,1);
        axis tight
    end
    
    % Applying hit and count method to remove outliers
    [x_final,y_final,z_final] = hit_and_count(x_org,...
        y_org,z_org,binsize,numer_atom_in_bin);
    
    if (tt == StartStep && graphcheck == 1)
        figure;
        scatter3(x_final,y_final,z_final,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .25 .25]);
        view(45,45)
        size2 = size(x_final,1);
        str = sprintf('After hit-and-count (%d points were removed)', ...
            size1-size2);
        title(str);
        axis tight
    end
    
    % Applying fine precission method to remove near droplet outliers
    if (fine_percision_check == 1)
        [x_final,y_final,z_final] = fine_percision(x_final,y_final,...
            z_final,R_step,delta_R);
    end
    
    if (tt == StartStep && graphcheck == 1)
        figure;
        scatter3(x_final,y_final,z_final,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .25 .25]);
        view(45,45)
        size3 = size(x_final,1);
        str = sprintf('After hit-and-count and fine-percision (%d points were removed)', ...
            size2-size3);
        title(str);
        axis tight
    end
    
    % removing pints close to surface and those at the top of the droplet
    % to improve the fitting 
    j = find(z_final<max(z_final)-10 & z_final>min(z_final)+5);
    
    x_final = x_final(j);
    y_final = y_final(j);
    z_final = z_final(j);
    
    if (tt == StartStep)
        zmax = max(z_final);
        zmin = min(z_final);
    end
    
    [theta(ii),x,y] = ...
        R_Ruijter(x_final,y_final,z_final,zmin,zmax,dz,dr);
    
    x_ave = cat(2,x_ave,x);
    y_ave = cat(2,y_ave,y);
    
end
fprintf('\n');      

theta_ave_1 = mean(nonzeros(theta));

range = 0:180;
[his] = histc(theta,range);

s0=his;
s1=conv(s0,hanning(30),'same');

file(:,1)=0:180;
file(:,2)=s0;
file(:,3)=s1;

dlmwrite('../Results/de_Ruijter.txt',file,'delimiter','\t');

% Fit to circle
F = @(R,xdata) sqrt(R(1)^2-(xdata).^2)-R(2);
R0 = [0,0];
x = x_ave';
y = y_ave';
Rmax = sqrt(max(x_ave)^2+max(y_ave)^2);
opts = optimset('Display','off');
[R,resnorm] = lsqcurvefit(F,R0,x,y,[Rmax/2,-5],[Rmax,max(y)],opts);

x0 = sqrt(R(1)^2-R(2)^2);
diff = -x0/sqrt(R(1)^2-x0^2);
% This average contact angle may be different from the average value
% calculated from the contact angle distribution. This value calculated
% through the circular fitting to all boundry points calculaed for all
% snapshots.
theta_ave_2 = abs(atan(diff)*180/pi);

% drawing angular distribution
figure;
plot(0:180,100*s1/sum(s1));
xlabel('Contact Angle (degree)');
ylabel('Probability (%)');
str(1) = {'de Ruijter et al. method'};
str(2) = {sprintf('\\theta_{ave}(from distribution)=%f',theta_ave_1)};
str(3) = {sprintf('\\theta_{ave}(circular fitting to all points)=%f',theta_ave_2)};
title(str);

figure;
scatter(x,y);
hold on;
r1 = 0:0.001:x0;
z1 = sqrt(R(1)^2-(r1).^2)-R(2);
plot(r1,z1,'r-');

xlabel('r');
ylabel('z');
title('Fitting to ensemble average of droplet boundary points');
