function [u,alpha]=eig_vec_householder(x,k)
%function [u,alpha]=householder_vec(x,k)
% y=H*x, H=I-2*u*u'/(u'*u), u=x-y
% i<k   y(i)=x(i)   u(i)=0
% i=k   y(k)=s*eta  u(k)=x(k)-y(k)=x(k)-s*eta=-s*xi
% i>k   y(i)=0      u(i)=x(i)
% eta=norm2(x(k:end))
% s=sgn(x(k))
% xi=-u(k)/s=-(x(k)-s*eta)/s=eta-x(k)/s=eta-|x(k)|
%   =(eta^2-|x(k)|^2)/(eta+|x(k)|)=(|x(k+1)|^2+..+|x(n)|^2)/(eta+|x(k)|)
% u'*u=2*xi*eta
% alpha=2/(u'*u)=1/(xi*eta)

%------------ init
if nargin<2
    k=1;
end
x=x(:);
n=length(x);
u=zeros(n,1);
alpha=0;

%----------- norm
xi=sum(abs(x((k+1):end)).^2);
eta=sqrt(abs(x(k))^2+xi);
if eta==0
    xi=eta-abs(x(k));
else
    xi=xi/(abs(x(k))+eta);
end
u((k+1):end)=x((k+1):end);
if x(k)==0
   u(k)=-xi;
else
   u(k)=-xi/abs(x(k))*x(k);
end
%----------- alpha
if xi==0 || eta==0
    alpha=0;
else
    alpha=1.0/(xi*eta);
end