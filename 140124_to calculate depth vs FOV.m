clear all

Depth=0.2;      %mm

R=15;          %mm

X0=500;

X=0:0.01:1000;     %mm

Z=R-real((R^2-(X-X0).^2).^0.5);

FOV=max(X(Z<Depth))-min(X(Z<Depth));

plot(X,Z);
