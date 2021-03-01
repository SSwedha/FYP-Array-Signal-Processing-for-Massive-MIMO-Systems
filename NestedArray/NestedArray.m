% 2 level Nested array
fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=5;
M2=5;
d=lambda/2;
d1=lambda/2;
d2=(M1+1)*d1;
SNR=20;
angle=[10 40 80; 10 30 40].*pi/180;
A1=getManifoldMatrixA1(angle,M1,M2,d1,lambda);
Phi=getPhi(angle,d,lambda);
A2=A1*Phi;


function A=getManifoldMatrixA1(angle,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sin(angle(1,:)).*sin(angle(2,:))));
end
function Phi=getPhi(angle,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cos(angle(1,:)).*sin(angle(2,:)))));
end