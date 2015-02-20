function [group,group_y] = group_matching_v2(N,NSh,Nmax,thresh,H,W,curr_pos,tBlocks,yBlocks,mode)

group=struct('data',zeros(N,N,Nmax),'pos',zeros(1,2,Nmax),'w',1,'size',0);
if(mode==1)
    group_y=group;
else
    group_y=0;
end

if(mode==0)
    rBlock=tBlocks{(curr_pos(1)-1)*(W-N+1)+curr_pos(2)};
else
    rBlock=yBlocks{(curr_pos(1)-1)*(W-N+1)+curr_pos(2)};
end
positions=zeros(1,2);
distances=zeros(1,1);
size=0;
lPos=curr_pos(2)+N/2-NSh;
rPos=curr_pos(2)+N/2+NSh;
uPos=curr_pos(1)+N/2-NSh;
dPos=curr_pos(1)+N/2+NSh;
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
group.size=val;
if(mode==1)
    group_y.size=val;
end
for m=1:val
    bpos=(positions(m,1)-1)*(W-N+1)+positions(m,2);
    group.data(:,:,m)=tBlocks{bpos};
    group.pos(:,:,m)=positions(m,:);
    if(mode==1)
        group_y.data(:,:,m)=yBlocks{bpos};
        group_y.pos(:,:,m)=positions(m,:);
    end
end
return