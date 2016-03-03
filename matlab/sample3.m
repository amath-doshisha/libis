% machine epsilon of multi precision
clear all;
prec=512;
set_default_prec(prec);
disp(sprintf('The machine epsilon on prec=%d is e=%.4e',prec,2^(-prec+1)));
disp('Finding it by calculations...')
x=2;
x=multi(x);
while x~=1
    e=(x-1)/2;
    x=1+e;
end
e=e*2
disp('Check: (1+e)-1=e')
(1+e)-1
disp('Check: (1+e/2)-1=0')
(1+e/2)-1