function angle_fit=Santiso_Angle(x_cor,y_cor,z_cor,r_cut)

angle_fit = zeros(size(x_cor,1),1);

for i=1:size(x_cor,1)
    r = sqrt((x_cor-x_cor(i)).^2+(y_cor-y_cor(i)).^2+...
        (z_cor-z_cor(i).^2));
    index = find(r<r_cut);
    if (size(index,1)>1);
        angle=fit_plane(x_cor,y_cor,z_cor, index);
        angle_fit(i)=angle;
    else
        angle_fit(i)=0;
    end
end

angle_fit=nonzeros(angle_fit);
