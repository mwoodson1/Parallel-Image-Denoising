function [test_groups] = collab_filt2(T2D,N,ht_thresh,T1D_all,groups,estim_var,groups_y,mode)

test_groups=groups;
test_groups_y=groups_y;
% collaborative filtering
for i=1:numel(test_groups)
    gsize=test_groups(i).size;
    level=floor(log2(gsize))+1;
    act_size=2^(level-1);
    test_groups(i).size=act_size;
    for j=1:N % filter 1D
        for k=1:N
            aux=reshape(test_groups(i).data(j,k,1:act_size),act_size,1);
            test_groups(i).data(j,k,1:act_size)=reshape(T1D_all{level}*aux,1,1,act_size);
            if(mode==1)
                 aux=reshape(test_groups_y(i).data(j,k,1:act_size),act_size,1);
                 test_groups_y(i).data(j,k,1:act_size)=reshape(T1D_all{level}*aux,1,1,act_size);
            end
        end
    end
    if(mode==0)
        test_groups(i).data=wthresh(test_groups(i).data,'h',ht_thresh); % Hard Thresholding   
        if(max(max(isnan(test_groups(i).data)))~=0)
            i
        end
        Nhar=nnz(test_groups(i).data);
        if(Nhar>0)
            test_groups(i).w=1/(estim_var*Nhar);
        else
            test_groups(i).w=1;
        end
    else
        W=abs(test_groups_y(i).data).^2./(abs(test_groups_y(i).data).^2+estim_var); % Wiener filter coeffs
        test_groups(i).data=W.*test_groups(i).data; % do the Wiener Filtering
        sum=0;
        for j=1:act_size
           sum=sum+norm(W(:,:,j),2); 
        end
        test_groups(i).w= 1/(estim_var*sum^2);
    end
    for j=1:N % Invert 1D transform
        for k=1:N
            aux=reshape(test_groups(i).data(j,k,1:act_size),act_size,1);
            test_groups(i).data(j,k,1:act_size)=reshape(T1D_all{level}\aux,1,1,act_size);
        end
    end
    for j=1:act_size % Invert 2D Transform
        test_groups(i).data(:,:,j)=T2D\test_groups(i).data(:,:,j)/T2D';
    end
end
return