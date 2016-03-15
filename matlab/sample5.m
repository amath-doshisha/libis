%function sample5
clear all
close all
n=10;
A=rand(n)
lambda=sort(eig(A))
plot(real(lambda),imag(lambda),'o')
axis square
grid on