%------------------------------------------------------------------------------------------------------------------------------------
% - functionality: calculates correlation between two vectors with same
%   lengths
% - input: 2 vectors with same lengths
% - output: correlation of the two vectors
%------------------------------------------------------------------------------------------------------------------------------------
function [ k ] = correlation( x,y )
mx=mean(x);
my=mean(y);
e=0;
for i=1:size(x,2)
    e=e+ (x(i)-mx).*(y(i)-my);
end
xbar=x-mx;
ybar=y-my;
k=e/(norm(xbar)*norm(ybar));
if isnan(k)
    k=0;
end

end

