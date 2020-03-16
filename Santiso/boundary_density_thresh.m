function density_thresh=boundary_density_thresh(density,ro_c)

[m,n,q]=size(density);

%boundary recognition

density_thresh=density;
for j=2:m-1
    for p=2:n-1
        for k=2:q-1
            if density(j,p,k)<ro_c
                density_thresh(j,p,k)=0;
            else
                v=[density(j-1,p,k-1);density(j-1,p-1,k-1);density(j,p-1,k-1);density(j+1,p+1,k-1);...
                    density(j-1,p,k+1);density(j-1,p-1,k+1);density(j,p-1,k+1);density(j+1,p+1,k+1);...
                    density(j-1,p,k);density(j-1,p-1,k);density(j,p-1,k);density(j+1,p+1,k)];
                if (min(v)>ro_c || max(v) < ro_c)
                    density_thresh(j,p,k)=0;
                end
                
                if (k==2 && density(j,p,k-1)<ro_c)
                    density_thresh(j,p,k-1)=0;
                else
                    v1=[density(j-1,p,k-1);density(j+1,p,k-1);...
                        density(j,p-1,k-1);density(j,p+1,k-1)];
                    if (min(v1)>ro_c || max(v1) < ro_c)
                        density_thresh(j,p,k-1)=0;
                    end
                end
                
            end

        end
    end
end

%remove the borders
density_thresh(1,:,:)=0;
density_thresh(end,:,:)=0;
density_thresh(:,1,:)=0;
density_thresh(:,end,:)=0;
density_thresh(:,:,1)=0;
density_thresh(:,:,end)=0;
