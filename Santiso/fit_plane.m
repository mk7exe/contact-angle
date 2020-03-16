function angle=fit_plane(x_cor,y_cor,z_cor, index)
% Second part of Santiso's algorithm
%step (a) was to find molecules within a given cutoff radious from molecule
%i: moluclues within the radius marked in index vector.
R=[x_cor(index) y_cor(index) z_cor(index)];
    r_average=[mean(R(:,1)) mean(R(:,2)) mean(R(:,3))];% This is step (b)
    R=R-repmat(r_average,length(index),1); % This is step (c)
    %%%%%%%%%%%%%%
    % This is step (d)
   W=zeros(3,3);
    for j=1:length(index)
    W=W+R(j,:)'*R(j,:);
    end
    %%%%%%%%%%%%%%%%%%%
    %This is step (e)
    [V,~]=eigs(W);
    n = V(:,3);
    angle= acos(n(3))*180/pi;
