%function sample4
clear all

%% initializations
tic
N=2;
M=4;
disp('�y�p�����[�^�̐ݒu�z');
disp(sprintf('�ѕ�: N=%d M=%d',N,M));
lambda=[1 2 3 4 5 6 7 8 9 10];
disp(sprintf('�w��ŗL�l: [%s]',num2str(lambda)));
%prec=53;lambda=double(lambda);
prec=512;set_default_prec(prec);lambda=multi(lambda);
disp(sprintf('�v�Z���x: %d bits',prec));

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
A=sample_matTN(lambda,N,M);
disp('�������ꂽ�s��: A=');
disp(num2str(A,'%.2f     '));
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp('�{���x�ɃL���X�g���ꂽ�s��double(A)�̌ŗL�l��{���x��eig�Ōv�Z: ');
l=sort(eig(double(A)));
disp(num2str(l,'%.40f'))
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp(sprintf('�������ꂽ�s��A�𑽔{��(%d bits)��eig�Ōv�Z:',prec*64));
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
toc

%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------');
tic
disp(sprintf('�{���x�ɃL���X�g���ꂽ�s��double(A)�𑽔{��(%d bits)��eig�Ōv�Z:',prec*64));
A=double(A);
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
toc