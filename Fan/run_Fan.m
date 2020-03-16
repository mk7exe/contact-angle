%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs the calculation of contact angle of a liquid droplet
% on a solid surface through the method of Fan and Cagin (1991). The volume
% and surface area of the liquid droplet is calculated through convex hull
% triangulation (function Convex_VS)
%
% Author: Mohammad Khalkhali, Sep 2016
% References:
% Khalkhali et al. J. Chem. Phys. (2017)
% C. F. Fan and T. Cain, The Journal of Chemical Physics 103, 9053 (1995).
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

theta = zeros((EndStep-StartStep)/step+1,1);
ii = 0;
V_ave = 0;
S_ave = 0;
Z_max_ave = 0;

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
    
    [V,S] = Convex_VS(x_final,y_final,z_final);
    V_ave = V_ave + V;
    S_ave = S_ave + S;
    Z_max_ave = Z_max_ave + max(z_final);
    
    P = [pi/3 0 -S 4*V];
    h = roots(P);
    h = h(find(imag(h)==0.0) & h < max(z_final)+10 & h > 0);
    R = (S+pi*h^2)/(4*pi*h);
    theta(ii) = acos(1-h/R)*180/pi;
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

dlmwrite('../Results/Fan.txt',file,'delimiter','\t');

S_ave = S_ave/((EndStep-StartStep)/step+1);
V_ave = V_ave/((EndStep-StartStep)/step+1);
Z_max_ave = Z_max_ave/((EndStep-StartStep)/step+1);

P = [pi/3 0 -S_ave 4*V_ave];
h = roots(P);
h = h(find(imag(h)==0.0) & h < Z_max_ave+10 & h > 0);
R = (S_ave+pi*h^2)/(4*pi*h);
% This average contact angle may be different from the average value
% calculated from the contact angle distribution. This value calculated the
% ensemble average of the volume and surface of droplet similar to the
% original publication
theta_ave_2 = acos(1-h/R)*180/pi;

% drawing angular distribution
figure;
plot(0:180,100*s1/sum(s1));
xlabel('Contact Angle (degree)');
ylabel('Probability (%)');
str(1) = {'Fan and Cagin meathod'};
str(2) = {sprintf('\\theta_{ave}(from distribution)=%f',theta_ave_1)};
str(3) = {sprintf('\\theta_{ave}(original publication)=%f',theta_ave_2)};
title(str);
