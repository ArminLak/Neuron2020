function map=make2color(c1,c2,varargin)


if nargin > 2
    N = varargin{1};
else
    N = 64;
end

%c1=[1 1 0.2];
%c2=[0.2 0.2 1];
cd=[(c2(1)-c1(1))/(N-1) (c2(2)-c1(2))/(N-1) (c2(3)-c1(3))/(N-1)];
for (i=1:N)
   map(i+1,:)=c1+cd*i;
end
map(1,:)  = c1;
map(N+2,:) = c2;
