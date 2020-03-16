%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads a LAMMPS trajectory file and returns the coordinates 
% of the centre of masses of water molecules. The LAMMPS trajectory file
% should contain just one snapshot. Oxygen and Hydrogen atoms of water
% should be the last two atom types in the trajectory file (Hw is the last
% atom type). Box_x, Box_y and Box_z are dimensions of the simulation box.
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

function [fileID,x_org,y_org,z_org] = read_LAMMPS_traj(Filename,Box_x,...
    Box_y,Box_z,format_spec)

    fileID = fopen(Filename,'r'); %opening the trajectory file to read
    
    if (fileID < 0)
        return
    end
    
    % reading the trajectory file
    datacell = textscan(fileID, format_spec, 'HeaderLines', 9);
    
    fclose(fileID);
    
    Data = cell2mat(datacell);
    [dummy, index] = sort(Data(:,1));
    
    ID_org = Data(index,1);
    Type_org = Data(index,2);
    x_org = Data(index,3);
    y_org = Data(index,4);
    z_org = Data(index,5);
    
    % Keeping only the water molecules: Ow and Hw
    max_type=max(Type_org);
    i=find(Type_org>=max_type-1);% Ow and Hw are the last two atom types
    
    ID_org=ID_org(i);
    Type_org=Type_org(i);
    x_org=x_org(i);
    y_org=y_org(i);
    z_org=z_org(i);
      
    j = find(Type_org == max_type-1); %O atoms in water
    
    % update positions of H atoms considering the periodic boundary
    % condition
    x_org(j+1) = x_org(j+1) - Box_x*round((x_org(j+1)-x_org(j))/Box_x);
    x_org(j+2) = x_org(j+2) - Box_x*round((x_org(j+2)-x_org(j))/Box_x);
    y_org(j+1) = y_org(j+1) - Box_y*round((y_org(j+1)-y_org(j))/Box_y);
    y_org(j+2) = y_org(j+2) - Box_y*round((y_org(j+2)-y_org(j))/Box_y);
    z_org(j+1) = z_org(j+1) - Box_z*round((z_org(j+1)-z_org(j))/Box_z);
    z_org(j+2) = z_org(j+2) - Box_z*round((z_org(j+2)-z_org(j))/Box_z);

    % calculating centre of mass of water molecules    
    x_org = (16*x_org(j)+x_org(j+1)+x_org(j+2))/18;
    y_org = (16*y_org(j)+y_org(j+1)+y_org(j+2))/18;
    z_org = (16*z_org(j)+z_org(j+1)+z_org(j+2))/18;
    
end
