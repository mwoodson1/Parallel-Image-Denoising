function [groups] = group_matching(T2D,N,Nstep,NS,NSh,z,Nmax,thresh,H,W)

tBlocks=struct('tBlock',zeros(N),'pos',[0 0]);

% pre compute 2D transform of all blocks
for i=NS-1:NS-1+H-N
    for j=NS-1:NS-1+W-N
        tsize=numel(tBlocks)+1;
        if(i==NS-1 && j==NS-1)
            tsize=tsize-1;
        end
        tBlocks(tsize).tBlock=T2D*z(i:i+N-1,j:j+N-1);
        tBlocks(tsize).pos=[i-NS+2 j-NS+2];
    end
end

groups=struct('data',zeros(N,N,Nmax),'pos',zeros(1,2,Nmax),'w',1,'size',0);

% group matching
for i=NS-1:Nstep:NS-1+H-N
    for j=NS-1:Nstep:NS-1+W-N
        pos=(i-NS+1)*(W-N+1)+j-NS+2;
        rBlock=tBlocks(pos).tBlock;
        gsize=size(groups,2)+1;
        if(j==NS-1 && i==NS-1)
            gsize=gsize-1;
        end
        groups(gsize).size=1;
        groups(gsize).data(:,:,1)=rBlock;
        groups(gsize).pos(:,:,1)=tBlocks(pos).pos;
        lPos=j+N/2-NSh;
        rPos=j+N/2+NSh;
        uPos=i+N/2-NSh;
        dPos=i+N/2+NSh;
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
                if(l<NS-1 || k<NS-1 || l>NS-1+H-N || k>NS-1+H-N)
                    continue;
                end;
                cBlock=tBlocks((k-NS+1)*(W-N+1)+l-NS+2).tBlock;
                cPos=tBlocks((k-NS+1)*(W-N+1)+l-NS+2).pos;
                if(norm((rBlock-cBlock))^2<thresh) % if cBlock is close enough to rBlock, group it
                    groups(gsize).size=groups(gsize).size+1;
                    groups(gsize).data(:,:,groups(gsize).size)=cBlock;
                    groups(gsize).pos(:,:,groups(gsize).size)=cPos;
                    if(groups(gsize).size==Nmax)
                        filled=1;
                    end
                end
            end
        end
    end
end
return