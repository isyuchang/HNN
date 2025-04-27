function [LL,HL,LH,HH] = Haart2(x)
% forward haar transform 2d
[M,N,~] = size(x);
if mod(M,2)==1 && mod(N,2)==1
    x = [x,x(:,end,:);x(end,:,:),x(end,end,:)];
elseif mod(M,2)==1 && mod(N,2)==0
    x = [x;x(end,:,:)];
elseif mod(M,2)==0 && mod(N,2)==1
    x = [x,x(:,end,:)];
end
[m,n,~] = size(x);
LL = ( x(1:2:m,1:2:n,:)+x(2:2:m,1:2:n,:)+x(1:2:m,2:2:n,:)+x(2:2:m,2:2:n,:))/2; % a
HL = (-x(1:2:m,1:2:n,:)-x(2:2:m,1:2:n,:)+x(1:2:m,2:2:n,:)+x(2:2:m,2:2:n,:))/2; % v
LH = (-x(1:2:m,1:2:n,:)+x(2:2:m,1:2:n,:)-x(1:2:m,2:2:n,:)+x(2:2:m,2:2:n,:))/2; % h
HH = ( x(1:2:m,1:2:n,:)-x(2:2:m,1:2:n,:)-x(1:2:m,2:2:n,:)+x(2:2:m,2:2:n,:))/2; % d
end