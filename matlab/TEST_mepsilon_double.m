clear all
close all

x=2;
y=x-1;
while x>1
    disp(sprintf('x=%.20e y=%.3e',x,y));
    y=y/2;
    x=1+y;
end
y=y*2;
disp(sprintf('machine epsilon is %.3e',y));
