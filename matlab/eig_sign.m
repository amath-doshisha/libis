function [A]=eig_sign(a,b)
%#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
if b>=0.0
   A=abs(a);
else
   A=-abs(a);
end