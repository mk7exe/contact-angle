%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function applies the hit-and-count algorithm to remove outliers from
% a data set
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

function [x_final,y_final,z_final] = hit_and_count(x_org,y_org,z_org,...
    binsize,numer_atom_in_bin)
% binning of the data using binsize

x_bin=floor(x_org/binsize)+1;
y_bin=floor(y_org/binsize)+1;
z_bin=floor(z_org/binsize)+1;

% applying hit-and-run in X direction
B=unique(x_bin);
C=histc(x_bin(:),B);
e=find(C>=numer_atom_in_bin);
KK=find(x_bin>=B(e(1)) & x_bin<=B(e(end)));

x=x_bin(KK);
y=y_bin(KK);
z=z_bin(KK);
% ID_drop=ID_org(KK);
% Type_drop=Type_org(KK);
x_drop=x_org(KK);
y_drop=y_org(KK);
z_drop=z_org(KK);

% applying hit-and-run in Y direction
B=unique(y);
C=histc(y(:),B);
e=find(C>=numer_atom_in_bin);
KK=find(y>=B(e(1)) & y<=B(e(end)));

x=x(KK);
y=y(KK);
z=z(KK);
% ID_drop=ID_drop(KK);
% Type_drop=Type_drop(KK);
x_drop=x_drop(KK);
y_drop=y_drop(KK);
z_drop=z_drop(KK);

% applying hit-and-run in Z direction
B=unique(z);
C=histc(z(:),B);
e=find(C>=numer_atom_in_bin);
KK=find(z>=B(e(1)) & z<=B(e(end)));

x=x(KK);
y=y(KK);
z=z(KK);
% ID_final=ID_drop(KK);
% Type_final=Type_drop(KK);
x_drop=x_drop(KK);
y_drop=y_drop(KK);
z_drop=z_drop(KK);

x_final=x_drop;
y_final=y_drop;
z_final=z_drop;


end
