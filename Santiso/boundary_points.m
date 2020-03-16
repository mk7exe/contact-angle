function [x_boundary, y_boundary, z_boundary]=boundary_points(x,y,z,density_thresh,bin_size)

x_boundary=[];
y_boundary=[];
z_boundary=[];

dx=bin_size;
x_bin=round(x./dx);
y_bin=round(y./dx);
z_bin=round(z./dx);
z_bin=[z_bin;1;1];

minx=min(x_bin);

miny=min(y_bin);

minz=min(z_bin);

for j=1:size(density_thresh,1)
    
     wwx=find(x>minx*dx+(j-1)*dx & x<minx*dx+j*dx);
   for k=1:size(density_thresh,2)
       
       wwy=find(y>miny*dx+(k-1)*dx & y<miny*dx+k*dx);
       for l=1:size(density_thresh,3)
           
          wwz=find(z>minz*dx+(l-1)*dx & z<minz*dx+l*dx);
          if density_thresh(j,k,l)~=0
           pp1=intersect(wwx,wwy);
          pp2=intersect(pp1,wwz);
              x_boundary=[x_boundary;x(pp2)];
              y_boundary=[y_boundary;y(pp2)];
              z_boundary=[z_boundary;z(pp2)];
          end
       end
   end
end
           
           

end


