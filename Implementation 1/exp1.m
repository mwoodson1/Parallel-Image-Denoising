close all;
y = imread('Cameraman256.png');
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
estim_var=var(z,0,1)*var(z,0,2);
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
thresh = 2500*(N^2)/(255*255);               
ht_thresh = 2.7*estim_var/N;
tic
[groups,~] = group_matching5(T2D,N,Nstep,NSh,z,Nmax,thresh,H,W,0,0);
toc
[filt_groups] = collab_filt2(T2D,N,ht_thresh,T1D_all,groups,estim_var,0,0);
toc
[y_hat] = aggregation(filt_groups,N,H,W);
toc
figure, imshow(y_hat)
BASIC_PSNR=10*log10(1/mean((y(:)-y_hat(:)).^2)) 
% END OF FIRST STEP (BASIC ESTIMATE)
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
[groups,groups_y] = group_matching5(T2D,N,Nstep,NSh,z,Nmax,thresh,H,W,y_hat,1);
toc
[filt_groups] = collab_filt2(T2D,N,ht_thresh,T1D_all,groups,estim_var,groups_y,1);
toc
[y_final] = aggregation(filt_groups,N,H,W);
toc
figure, imshow(y_final)
FINAL_PSNR=10*log10(1/mean((y(:)-y_final(:)).^2))