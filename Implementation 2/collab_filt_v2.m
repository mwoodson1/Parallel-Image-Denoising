function [filt_group] = collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,group_y,mode)
gsize=group.size;
level=floor(log2(gsize))+1;
act_size=2^(level-1);
group.size=act_size;
for j=1:N % filter 1D
    for k=1:N
        aux=reshape(group.data(j,k,1:act_size),act_size,1);
        group.data(j,k,1:act_size)=reshape(T1D_all{level}*aux,1,1,act_size);
        if(mode==1)
             aux=reshape(group_y.data(j,k,1:act_size),act_size,1);
            group_y.data(j,k,1:act_size)=reshape(T1D_all{level}*aux,1,1,act_size);
        end
    end
end
if(mode==0)
    group.data=wthresh(group.data,'h',ht_thresh); % Hard Thresholding   
    Nhar=nnz(group.data);
    if(Nhar>0)
        group.w=1/(estim_var*Nhar);
    else
        group.w=1;
    end
else
    W=abs(group_y.data).^2./(abs(group_y.data).^2+estim_var); % Wiener filter coeffs
    group.data=W.*group.data; % do the Wiener Filtering
    sum=0;
    for j=1:act_size
       sum=sum+norm(W(:,:,j),2); 
    end
    group.w= 1/(estim_var*sum^2);
end
for j=1:N % Invert 1D transform
    for k=1:N
        aux=reshape(group.data(j,k,1:act_size),act_size,1);
        group.data(j,k,1:act_size)=reshape(T1D_all{level}\aux,1,1,act_size);
    end
end
for j=1:act_size % Invert 2D Transform
    group.data(:,:,j)=T2D\group.data(:,:,j)/T2D';
end
filt_group=group;
return