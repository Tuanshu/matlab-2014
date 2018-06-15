clear all

Beam_size=8;        %mm (diameter)
f=80;               %mm

X=0:0.1:20;
X0=10;
Gaussian=gaussmf(X,[Beam_size/8 X0]);
plot(X,Gaussian);

Sig=Beam_size/8;

