function [y_hat] = aggregation(test_groups,N,H,W)
my_y_hat=zeros(H,W);
agg_weig=zeros(H,W);
beta=2.0;
Wwin2D=kaiser(N,beta)*kaiser(N,beta)';
for i=1:numel(test_groups)
    for j=1:test_groups(i).size
        sBlock=test_groups(i).data(:,:,j);
        xPos=test_groups(i).pos(1,1,j);
        yPos=test_groups(i).pos(1,2,j);
        my_y_hat(xPos:xPos+N-1,yPos:yPos+N-1)=my_y_hat(xPos:xPos+N-1,yPos:yPos+N-1)+test_groups(i).w*Wwin2D.*sBlock;
        agg_weig(xPos:xPos+N-1,yPos:yPos+N-1)=agg_weig(xPos:xPos+N-1,yPos:yPos+N-1)+test_groups(i).w*Wwin2D;
    end
end
y_hat=my_y_hat./agg_weig;
return