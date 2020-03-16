%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs the calculation of contact angle of a liquid droplet
% on a solid surface using the method of Santiso et al. In this method
% boundary points of the liquid droplet are found through a 3D density
% profile. Then contact angle is calculated as the average value of angles
% between the normal of the planes fitted to boundary points and normal to
% the solid surface (see references for mor details)
%
% Authors: 
% Mohammad Khalkhali
% Nasser Kazemi 
% (Sep 2016)
%
% References:
% Khalkhali et al. J. Chem. Phys. (2017)
% E. Santiso, C. Herdes, and E. Muller, Entropy 15, 3734 (2013).
%
% Functions:
% density_cal: Calculated the 3D density profiles
% boundary_density_thresh: finds the bins located at the droplet boundary
% boundary_points: finds points corresponding to boundary bins
% Santiso_Angle: calculates contact angle through fitting 3D planes to
% boundary points
% fit_plane: fits a 3D plane to a data set
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

addpath('../');
trajName = '../Graphite_Water(3nm)/Graphite_Water.lammpstrj';
% This format_spec is needed when reading lammps trajectories in
% read_LAMMPS_traj function. There are 5 entries at each line of our
% trajectories
format_spec = '%f %f %f %f %f';
StartStep = 4990000;
EndStep = 5000000;
step = 1000;

binsize=5;
%cut off limit for density
ro_c=0.01;

r_cut = 5;%cut-off distance for plane fitting

z_max_coef = 0.25;%portion of droplet height to be considered for fitting

Distribution=[];
ii = 0;

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
    
    %calculate density and normalize it
    density=density_cal(x_org,y_org,z_org,binsize);
    
    density=density./binsize^3;%Calculating number density
    density=density./max(max(max(density)));
    
    % calculate boundary points
    density_thresh=boundary_density_thresh(density,ro_c);
    
    clear density
    
    [x_boundary, y_boundary, z_boundary]=...
        boundary_points(x_org,y_org,z_org,density_thresh,binsize);
    % clear density_thresh
    clear density_thresh
    
    x_cor=x_boundary(:);
    y_cor=y_boundary(:);
    z_cor=z_boundary(:);
    
    % clear density_thresh
    clear x_boundary
    clear y_boundary
    clear z_boundary
    
    %appy lower limit on z coeerdinate
    [ii,pp]=find(z_cor<=z_max_coef*max(z_cor));

    x_cor=x_cor(ii);
    y_cor=y_cor(ii);
    z_cor=z_cor(ii);
    
    % Drawing the selected points to check the patameters
    if (tt==StartStep)
        figure
        scatter3(x_cor,y_cor,z_cor);
        title(['Boundary points calculated at time step =' num2str(tt)]);
    end
    %define angle distribution
    
    angle_fit=Santiso_Angle(x_cor,y_cor,z_cor,r_cut);
    
    clear x_cor
    clear y_cor
    clear z_cor
    
    angle_fit = angle_fit';
    Distribution = cat(2,Distribution,angle_fit(1,:));
    
    clear angle_fit
end
fprintf('\n');      

theta_ave_1 = mean(nonzeros(Distribution));

range = 0:180;
[his] = histc(Distribution,range);

s0=his;
s1=conv(s0,hanning(30),'same');

file(:,1)=0:180;
file(:,2)=s0;
file(:,3)=s1;

dlmwrite('../Results/Suntiso.txt',file,'delimiter','\t');

% drawing angular distribution
figure;
plot(0:180,100*s1/sum(s1));
xlabel('Contact Angle (degree)');
ylabel('Probability (%)');
str(1) = {'Santiso et al. method'};
str(2) = {sprintf('\\theta_{ave}=%f',theta_ave_1)};
title(str);
