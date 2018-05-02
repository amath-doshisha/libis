clear all
set_default_prec(53)
x=multi(665857);
y=multi(470832);
disp('The variables x and y are storaged by 53-bit.');
disp('x='); disp(x);
disp('y='); disp(y);
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 53-bit.');
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is NOT true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 64-bit.');
set_default_prec(64)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is NOT true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 128-bit.');
set_default_prec(128)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is true.');
disp('------');
disp('Compute z=x^4-4*y^2-4*y^4 by 256-bit.');
set_default_prec(256)
z=x^4-4*y^2-4*y^4;
disp('z='); disp(z);
disp('This result is true.');
disp('------');
