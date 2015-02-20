function [groups,groups_y] = group_matching5(T2D,N,Nstep,NSh,z,Nmax,thresh,H,W,y_hat,mode)

tsize=(H-N+1)*(W-N+1);
tBlocks=cell(tsize,1);
if(mode==1)
    yBlocks=cell(tsize,1);
end
for i=1:H-N+1
    for j=1:W-N+1
        tBlocks{(i-1)*(H-N+1)+j}=T2D*z(i:i+N-1,j:j+N-1)*T2D';
        if(mode==1)
            yBlocks{(i-1)*(H-N+1)+j}=T2D*y_hat(i:i+N-1,j:j+N-1)*T2D';
        end
    end
end

gsize=ceil((H-N+1)/Nstep)*ceil((W-N+1)/Nstep);
group=struct('data',zeros(N,N,Nmax),'pos',zeros(1,2,Nmax),'w',1,'size',0);
groups=repmat(group,gsize,1);
if(mode==1)
    groups_y=repmat(group,gsize,1);
else
    groups_y=0;
end
pos=1;
% group matching
for i=1:Nstep:H-N+1
    for j=1:Nstep:W-N+1
        if(mode==0)
            rBlock=tBlocks{(i-1)*(W-N+1)+j};
        else
            rBlock=yBlocks{(i-1)*(W-N+1)+j};
        end
        positions=zeros(1,2);
        distances=zeros(1,1);
        size=0;
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
                if(mode==0)
                    cBlock=tBlocks{(k-1)*(W-N+1)+l};
                else
                    cBlock=yBlocks{(k-1)*(W-N+1)+l};
                end
                dist=norm((rBlock-cBlock))^2;
                if(dist<thresh) % if cBlock is close enough to rBlock, group it
                    size=size+1;
                    distances(size)=dist;
                    positions(size,:)=[k l];
                end
            end
        end
        [~,ord]=sort(distances);
        positions=positions(ord,:);
        if(size>Nmax)
           val=Nmax;
        else
           val=size;
        end
        groups(pos).size=val;
        if(mode==1)
            groups_y(pos).size=val;
        end
        for m=1:val
            bpos=(positions(m,1)-1)*(W-N+1)+positions(m,2);
            groups(pos).data(:,:,m)=tBlocks{bpos};
            groups(pos).pos(:,:,m)=positions(m,:);
            if(mode==1)
                groups_y(pos).data(:,:,m)=yBlocks{bpos};
                groups_y(pos).pos(:,:,m)=positions(m,:);
            end
        end
        pos=pos+1;
    end
end
return