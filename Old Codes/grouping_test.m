close all;
y = imread('lena_crop.png');
y = im2double(y);
imshow(y)
H=size(y,1); % height |
W=size(y,2); % width --
N=8;
Nstep=3;
Nmax=16;
NS=35;
NSh=(NS-1)/2;
sigma=25;
z=y+(sigma/255)*randn(size(y));
figure, imshow(z)
hold on;
%PSNR=10*log10(1/mean((y(:)-z(:)).^2))
T2D=[0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
     0.219417649252501   0.449283757993216   0.449283757993216   0.219417649252501  -0.219417649252501  -0.449283757993216  -0.449283757993216  -0.219417649252501;
     0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846  -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284;
    -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284   0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846;
     0.707106781186547  -0.707106781186547                   0                   0                   0                   0                   0                   0;
                     0                   0   0.707106781186547  -0.707106781186547                   0                   0                   0                   0;
                     0                   0                   0                   0   0.707106781186547  -0.707106781186547                   0                   0;
                     0                   0                   0                   0                   0                   0   0.707106781186547  -0.707106781186547];   
                 

groups=struct('data',zeros(N,N,Nmax),'pos',zeros(1,2,Nmax),'dist',zeros(1,Nmax),'w',1,'size',1);
thresh=2500*(N^2)/(255*255);

i=20; j=60;
rBlock=T2D*z(i:i+N-1,j:j+N-1);
groups.data(:,:,1)=rBlock;
lPos=j+N/2-NSh;
if(lPos<1)
    lPos=1;
end
rPos=j+N/2+NSh;
if(rPos>W)
    rPos=W;
end
uPos=i+N/2-NSh;
if(uPos<1)
    uPos=1;
end
dPos=i+N/2+NSh;
if(dPos>W)
    dPos=H;
end
rectangle('Position',[j,i,N,N],'FaceColor','r');
filled=0;
for k=uPos:dPos-N+1
    if(filled==1)
        break;
    end
    for l=lPos:rPos-N+1
        if(filled==1)
            break;
        end
        if(l==j && k==i)
            continue;
        end
        cBlock=T2D*z(k:k+N-1,l:l+N-1);
        dist=norm((rBlock-cBlock))^2;
        if(dist<thresh) % if cBlock is close enough to rBlock, group it
            groups.size=groups.size+1;
            groups.dist(groups.size)=dist;
            groups.data(:,:,groups.size)=cBlock;
           % rectangle('Position',[l,k,N,N],'FaceColor','b');
            groups.pos(:,:,groups.size)=[k l];
        end
    end
end
[~,ord]=sort([groups.dist]);
test_groups=struct('data',groups.data(:,:,ord),'pos',groups.pos(:,:,ord),'dist',groups.dist(:,ord),'w',1,'size',groups.size);
for i=2:Nmax
    k=test_groups.pos(1,1,i);
    l=test_groups.pos(1,2,i);
    rectangle('Position',[l,k,N,N],'FaceColor','b');
end