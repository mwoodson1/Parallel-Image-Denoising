close all;
y = imread('cameraman256.png');
y = im2double(y);
sigma=25;
H=size(y,1); % height |
W=size(y,2); % width --
N=8;
Nstep=3;
Nmax=16;
NS=39;
NSh=(NS-1)/2;
imshow(y)
z=y+(sigma/255)*randn(size(y));
figure, imshow(z)
PSNR=10*log10(1/mean((y(:)-z(:)).^2))
% 2D-Bior1.5 transform matrix
T2D=[0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
     0.219417649252501   0.449283757993216   0.449283757993216   0.219417649252501  -0.219417649252501  -0.449283757993216  -0.449283757993216  -0.219417649252501;
     0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846  -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284;
    -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284   0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846;
     0.707106781186547  -0.707106781186547                   0                   0                   0                   0                   0                   0;
                     0                   0   0.707106781186547  -0.707106781186547                   0                   0                   0                   0;
                     0                   0                   0                   0   0.707106781186547  -0.707106781186547                   0                   0;
                     0                   0                   0                   0                   0                   0   0.707106781186547  -0.707106781186547];   
for h=1:log2(Nmax)+1;
    powh=2^(h-1);
    T1D = zeros(powh);
    for i = 1:powh
        T1D(:,i)=wavedec(circshift([1 zeros(1,powh-1)],[0 i-1]), log2(powh), 'haar');  %% construct transform matrix
    end
    T1D = (T1D' * diag(sqrt(1./sum(T1D.^2,2))))';
    T1D_all{h}=T1D;
end
thresh = 3000*(N^2)/(255*255);               
beta=2.0;
Wwin2D=kaiser(N,beta)*kaiser(N,beta)';
%%
tic
tsize=(H-N+1)*(W-N+1);
tBlocks=cell(tsize,1);
for i=1:H-N+1
    for j=1:W-N+1
        tBlocks{(i-1)*(H-N+1)+j}=T2D*z(i:i+N-1,j:j+N-1)*T2D';
    end
end

tvar=zeros(1,tsize*16);
for i=1:tsize
    tvar(1,16*i-15:16*i)=reshape(tBlocks{i}(5:8,5:8),1,16);
end
estim_var=1.4826*median(abs(tvar));
estim_var*255 % just to check if its near the noise value added
ht_thresh = 2.7*estim_var;

im_buf=zeros(H,W);
w_buf=zeros(H,W);
for i=1:Nstep:H-N+1
    for j=1:Nstep:W-N+1
        %group matching
        [group,~]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,0,0);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,0,0);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
end
% Corner Processing 
if(mod(H-N,3)~=0) 
    i=H-N+1;
    for j=1:Nstep:W-N+1
        %group matching
        [group,~]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,0,0);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,0,0);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
    j=H-N+1;
    for i=1:Nstep:W-N+1
        %group matching
        [group,~]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,0,0);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,0,0);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
    i=H-N+1;
    j=H-N+1;
    %group matching
    [group,~]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,0,0);
    %collaborative filtering
    [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,0,0);
    %aggregation
    for k=1:filt_group.size
        sBlock=filt_group.data(:,:,k);
        xPos=filt_group.pos(1,1,k);
        yPos=filt_group.pos(1,2,k);
        im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
        w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
    end
end

y_hat=im_buf./w_buf;
figure, imshow(y_hat)
BASIC_PSNR=10*log10(1/mean((y(:)-y_hat(:)).^2))
toc
% END OF FIRST STEP (BASIC ESTIMATE)
% 2D DCT transform
T2D =[ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
       0.490392640201615   0.415734806151273   0.277785116509801   0.097545161008064  -0.097545161008064  -0.277785116509801  -0.415734806151273  -0.490392640201615;
       0.461939766255643   0.191341716182545  -0.191341716182545  -0.461939766255643  -0.461939766255643  -0.191341716182545   0.191341716182545   0.461939766255643;
       0.415734806151273  -0.097545161008064  -0.490392640201615  -0.277785116509801   0.277785116509801   0.490392640201615   0.097545161008064  -0.415734806151273;
       0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274   0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274;
       0.277785116509801  -0.490392640201615   0.097545161008064   0.415734806151273  -0.415734806151273  -0.097545161008064   0.490392640201615  -0.277785116509801;
       0.191341716182545  -0.461939766255643   0.461939766255643  -0.191341716182545  -0.191341716182545   0.461939766255643  -0.461939766255643   0.191341716182545;
       0.097545161008064  -0.277785116509801   0.415734806151273  -0.490392640201615   0.490392640201615  -0.415734806151273   0.277785116509801  -0.097545161008064];
thresh = 400*(N^2)/(255*255);  
Nmax=32;
for h=1:log2(Nmax)+1;
    powh=2^(h-1);
    T1D = zeros(powh);
    for i = 1:powh
        T1D(:,i)=wavedec(circshift([1 zeros(1,powh-1)],[0 i-1]), log2(powh), 'haar');  %% construct transform matrix
    end
    T1D = (T1D' * diag(sqrt(1./sum(T1D.^2,2))))';
    T1D_all{h}=T1D;
end
%%
tsize=(H-N+1)*(W-N+1);
yBlocks=cell(tsize,1);
for i=1:H-N+1
    for j=1:W-N+1
        yBlocks{(i-1)*(H-N+1)+j}=T2D*y_hat(i:i+N-1,j:j+N-1)*T2D';
    end
end

im_buf=zeros(H,W);
w_buf=zeros(H,W);

for i=1:Nstep:H-N+1
    for j=1:Nstep:W-N+1
        %group matching
        [group,group_y]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,yBlocks,1);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,group_y,1);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
end
% Corner Processing 
if(mod(H-N,3)~=0) 
    i=H-N+1;
    for j=1:Nstep:W-N+1
        %group matching
        [group,group_y]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,yBlocks,1);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,group_y,1);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
    j=H-N+1;
    for i=1:Nstep:W-N+1
        %group matching
        [group,group_y]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,yBlocks,1);
        %collaborative filtering
        [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,group_y,1);
        %aggregation
        for k=1:filt_group.size
            sBlock=filt_group.data(:,:,k);
            xPos=filt_group.pos(1,1,k);
            yPos=filt_group.pos(1,2,k);
            im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
            w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
        end
    end
    i=H-N+1;
    j=H-N+1;
    %group matching
    [group,group_y]=group_matching_v2(N,NSh,Nmax,thresh,H,W,[i j],tBlocks,yBlocks,1);
    %collaborative filtering
    [filt_group]=collab_filt_v2(T2D,N,ht_thresh,T1D_all,group,estim_var,group_y,1);
    %aggregation
    for k=1:filt_group.size
        sBlock=filt_group.data(:,:,k);
        xPos=filt_group.pos(1,1,k);
        yPos=filt_group.pos(1,2,k);
        im_buf(xPos:xPos+N-1,yPos:yPos+N-1)=im_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D.*sBlock;
        w_buf(xPos:xPos+N-1,yPos:yPos+N-1)=w_buf(xPos:xPos+N-1,yPos:yPos+N-1)+filt_group.w*Wwin2D;
    end
end

y_final=im_buf./w_buf;
toc
figure, imshow(y_final)
FINAL_PSNR=10*log10(1/mean((y(:)-y_final(:)).^2))