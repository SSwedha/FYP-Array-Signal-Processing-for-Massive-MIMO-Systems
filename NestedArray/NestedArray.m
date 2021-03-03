% 2 level Nested array
fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=10;
M2=10;
M=M2+M1;
d=lambda/2;
d1=lambda/2;
d2=(M1+1)*d1;
SNR=20;
fs=10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*cos(2*pi*fc*t)'];
D=size(S,2);
sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
n2=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
angle=[10 40 80; 10 30 40].*pi/180;
A1=getManifoldMatrixA1(angle,M1,M2,d1,lambda);
Phi=getPhi(angle,d,lambda);
A2=A1*Phi;
x1=S*A1'+n1;
x2=S*A2'+n2;

R11=cov(x1,1);
z1=R11(:);
z1=z1(getIndexOfUniqueElements(M1,M2));
R12=cov(x2,1);
z2=R12(:);
z2=z2(getIndexOfUniqueElements(M1,M2));


function uidx=getIndexOfUniqueElements(M1,M2)
narray=[0:M1 ((2:M2).*(M1+1)-1)];
cr=[];
for i=1:length(narray)
    cr=[cr -narray+narray(i)];
end
[~,uidx,~]=unique(cr);
end
function A=getManifoldMatrixA1(angle,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sin(angle(1,:)).*sin(angle(2,:))));
end
function Phi=getPhi(angle,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cos(angle(1,:)).*sin(angle(2,:)))));
end
