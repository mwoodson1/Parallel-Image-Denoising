function [test_groups] = collab_filt(T2D,N,ht_thresh,T1D_all,groups,estim_var)

test_groups=groups;

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
        end
    end
    test_groups(i).data=wthresh(test_groups(i).data,'h',ht_thresh); % Hard Thresholding
    Nhar=nnz(test_groups(i).data);
    if(Nhar>0)
        test_groups(i).w=1/(estim_var*Nhar);
    else
        test_groups(i).w=1;
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