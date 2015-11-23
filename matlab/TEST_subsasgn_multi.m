clear all
close all

A=zeros(3,4,2,'multi')
A(:)=1:numel(A)
A(1:2:end)=0
A(:,5)=ones(3,1)
A(:,5,1)=ones(3,1)*2
A(:,6,1)=ones(3,1)*3
