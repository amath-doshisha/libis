clear all
close all

A=rrand(512,5,5)*2-1
b=rround(512,[1 2; 3 -1; 4 5; -1 2; 2 2])
x=A\b
r=A*x-b
rmax(rmax(rabs(r)))

B=czeros(512,5,5)
B=cones(512,5,5)
B=crand(512,5,5)
B=cround(512,[1 2 3; 4 5 6; 7 8 9])
B=cround(512,[1-i 2-2*i 3-3*i; 4-4*i 5-5*i 6-6*i; 7-7*i 8-8*i 9-9*i])
B=ceye(512,3,3)