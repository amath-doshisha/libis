clear all
close all

A=rrand(512,5,5)*2-1
b=rround(512,[1 2; 3 -1; 4 5; -1 2; 2 2])
x=A\b
r=A*x-b
rmax(rmax(rabs(r)))