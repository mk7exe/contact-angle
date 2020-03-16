%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs the calculation for contact angle of a liquid droplet
% on a solid surface using the method of Khalkhali et al. In this method
% contact angle is calculated through convex hull triangulation (see
% the referenced Publication for more details)
%
% Authors: 
% Mohammad Khalkhali
% Nasser Kazemi 
% (Sep 2016)
%
% References:
% Khalkhali et al. J. Chem. Phys. (2017)
%
% Functions:
% read_LAMMPS_traj: reads a LAMMPS trajectory file and returns the
% positions of centre of masses of water molecules.
% hit_and_count: runs hit and count algorithm to identify points in liquid
% droplet from those in the gas phase.
% fine_precision: applied the fine precision droplet identification process
% to remove near-droplet gas molecules.
% contact_angle: calculates the contact angles along the contact line using
% the convex hull triangulation. This function also returns the area of the
% base of the droplet and the area of corresponding triangles used to
% calculate each contact angle value.
% weighted_distribution: calculates weighted histogram of contact angles
% according to the area of corresponding triangles
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

% Parameters of fine_precision function (refer to the publication for details)
% fine presicion step of identifying droplet limits is usually redundant
% (hit and count step is sufficient). 
fine_precision_check = 0; % calls fine_precision function if is 1
delta_R = 5;
R_step = 1;

% Draws graphs corresponding to each step for the first time step. It is
% recommended to check if parameters droplet identification process are
% ajusted properly.
graphcheck = 0;

% Triangles with atleast one vertex position lower that z_min*max(z) are
% contributing in the contact angle calculation.
z_min = 0.08;

% This parameter is used to recognize triangles contributing to the area of
% the base of the droplet.
precision = 1.0;

Distribution=[];
Weight=[];
BaseA = 0;
            
for tt = StartStep:step:EndStep
    
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
    if (fine_precision_check == 1)
        [x_final,y_final,z_final] = fine_precision(x_final,y_final,...
            z_final,R_step,delta_R);
    end
    
    if (tt == StartStep && graphcheck == 1)
        figure;
        scatter3(x_final,y_final,z_final,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .25 .25]);
        view(45,45)
        size3 = size(x_final,1);
        str = sprintf('After hit-and-count and fine-precision (%d points were removed)', ...
            size2-size3);
        title(str);
        axis tight
    end
            
            
    % Calculating contact angle distribution and interfacial area using convex hull
    [K,V,angle,BaseA_temp] = contact_angle(x_final,y_final,z_final,...
        z_min,precision);
    if (tt == StartStep && graphcheck == 1)
        figure;
        hold on;
        trisurf(K,x_final,y_final,z_final,'facecolor','r','facealpha',0.5)
        scatter3(x_final,y_final,z_final,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .75 .75], ...
            'MarkerFaceAlpha',0.5);
        view(45,45)
        axis tight
    end
    
    BaseA = BaseA + BaseA_temp;
    Distribution = cat(2,Distribution,angle(1,:));
    Weight = cat(2,Weight,angle(2,:));
            
end
        
fprintf('\n');      
% Calculating weighted histogram of angular values
w = weighted_distribution(Distribution,Weight);
BaseA = 0.5*BaseA/((EndStep-StartStep)/step+1);

angle = 0:180;
s0=w;
s1=conv(w,hanning(30),'same');

theta_ave = sum(angle.*w)/sum(w);

% drawing angular distribution
% drawing angular distribution
figure;
plot(0:180,100*s1/sum(s1));
xlabel('Contact Angle (degree)');
ylabel('Probability (%)');
str(1) = {'Khalkhali et al. meathod'};
str(2) = {sprintf('\\theta_{ave}=%f',theta_ave)};
title(str);
file(:,1)=0:180;
file(:,2)=s0;
file(:,3)=s1;

output = sprintf('../Results/MK.txt');
dlmwrite(output,file,'delimiter','\t');
%         toc
