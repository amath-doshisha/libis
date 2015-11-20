clear all
close all

set_default_prec(128);

x=multi(2);
y=x-1;
while x>1
    disp(sprintf('x=%.20e y=%.3e',double(x),double(y)));
    y=y/2;
    x=1+y;
end
y=y*2;
disp(sprintf('machine epsilon is %.3e',double(y)));
