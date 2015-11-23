clear all
close all

set_default_prec(64)
auto_prec_enabled;

x=multi(665857)
y=multi(470832)
z=x*x*x*x-4*y*y-4*y*y*y*y
get_prec(z)
% exact: z=1
