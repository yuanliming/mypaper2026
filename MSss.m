function [A,B2] = MSss(dm, dd)
m=1/400;
d=1/200;
a1=-125;a2=-75;a3=-15;

A=[0 1 0 0 0 0;
   0 0 1 0 0 0;
   a1 a2 a3 0 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   a1*dm/(m+dm) a2*dm/(m+dm) (a3*dm+dd)/(m+dm) 0 0 -(d+dd)/(m+dm)];
B2=[0;0;0;0;0;-1/(m+dm)];
end
