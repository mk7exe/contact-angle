function density= density_cal(x,y,z,bin_size)
%x,y,z are data points in 3D and dx,dy,dz are bin size for density
%calculation

dx=bin_size;
x_bin=round(x./dx);
y_bin=round(y./dx);
z_bin=round(z./dx);
z_bin=[z_bin;1;1];

minx=min(x_bin);
maxx=max(x_bin);

miny=min(y_bin);
maxy=max(y_bin);

minz=min(z_bin);
maxz=max(z_bin);


density=zeros(maxx-minx+1,maxy-miny+1,maxz-minz+1);
for j=1:size(density,1)
    
     wwx=find(x>minx*dx+(j-1)*dx & x<minx*dx+j*dx);
   for k=1:size(density,2)
       
       wwy=find(y>miny*dx+(k-1)*dx & y<miny*dx+k*dx);
       for l=1:size(density,3)
           
          wwz=find(z>minz*dx+(l-1)*dx & z<minz*dx+l*dx);
          pp1=intersect(wwx,wwy);
          pp2=intersect(pp1,wwz);
          density(j,k,l)=length(pp2);
       end
   end
end

end
