function [groups] = group_matching2(T2D,N,Nstep,NSh,z,Nmax,thresh,H,W)

tBlocks=struct('tBlock',zeros(N),'pos',[0 0]);

% pre compute 2D transform of all blocks
for i=1:H-N+1
    for j=1:W-N+1
        tsize=numel(tBlocks)+1;
        if(i==1 && j==1)
            tsize=tsize-1;
        end
        tBlocks(tsize).tBlock=T2D*z(i:i+N-1,j:j+N-1);
        tBlocks(tsize).pos=[i j];
    end
end

groups=struct('data',zeros(N,N,Nmax),'pos',zeros(1,2,Nmax),'dist',zeros(1,Nmax),'w',1,'size',0);

% group matching
for i=1:Nstep:H-N+1
    for j=1:Nstep:W-N+1
        pos=(i-1)*(W-N+1)+j;
        rBlock=tBlocks(pos).tBlock;
        gsize=size(groups,2)+1;
        if(j==1 && i==1)
            gsize=gsize-1;
        end
        groups(gsize).size=1;
        groups(gsize).data(:,:,1)=rBlock;
        groups(gsize).pos(:,:,1)=tBlocks(pos).pos;
        lPos=j+N/2-NSh;
        rPos=j+N/2+NSh;
        uPos=i+N/2-NSh;
        dPos=i+N/2+NSh;
        if(lPos<1)
            lPos=1;
        end
        if(rPos>W)
            rPos=W;
        end
        if(uPos<1)
            uPos=1;
        end 
        if(dPos>W)
            dPos=H;
        end
        for k=uPos:dPos-N+1
            for l=lPos:rPos-N+1
                if(l==j && k==i)
                    continue;
                end
                cBlock=tBlocks((k-1)*(W-N+1)+l).tBlock;
                cPos=tBlocks((k-1)*(W-N+1)+l).pos;
                dist=norm((rBlock-cBlock))^2;
                if(dist<thresh) % if cBlock is close enough to rBlock, group it
                    groups(gsize).size=groups(gsize).size+1;
                    groups(gsize).dist(groups(gsize).size)=dist;
                    groups(gsize).data(:,:,groups(gsize).size)=cBlock;
                    groups(gsize).pos(:,:,groups(gsize).size)=cPos;
                end
            end
        end
        [~,ord]=sort([groups(gsize).dist]);
        groups(gsize)=struct('data',groups(gsize).data(:,:,ord),'pos',groups(gsize).pos(:,:,ord),'dist',groups(gsize).dist(:,ord),'w',1,'size',groups(gsize).size);
        if(groups(gsize).size>Nmax)
            groups(gsize)=struct('data',groups(gsize).data(:,:,1:Nmax),'pos',groups(gsize).pos(:,:,1:Nmax),'dist',groups(gsize).dist(:,1:Nmax),'w',1,'size',groups(gsize).size);
            groups(gsize).size=Nmax;
        end
    end
end
return