% DOA Matrix for Two parallel Nested Arrays

% Parameters
fc=2.4e9; % frequency of narrowband signal
c=physconst("lightspeed");
lambda=(c/fc); % wavelength
L1=20; % no.of elements in level1
L2=12; % no.of elements in level2
L=L2+L1; % total elemnts in subarray1
dx=lambda/2; % distance between subarrays
dy=lambda/2; % separation in level 1
Dy=(L1+1)*dy; % seperation in level 2
K=1; % no.of sources
sig_angle=[30,50]; % Angle of arrival (theta_1,phi_1;...;theta_k,phi_k)
SNR=200; % SNR of signal
fs=50*fc; % sampling frequency
N=2000; % no.of snapshots
t=0:1/fs:(N-1)*(1/fs);

% Signal Model
S=exp(1i*(2*pi*fc*t+(0:(pi/K):(K-1)*(pi/K))'));

A1=getManifoldMatrixA1(sig_angle(:,1),sig_angle(:,2),L1,L2,dy,lambda);
A2=A1*getPhi(sig_angle(:,1),sig_angle(:,2),dx,lambda);
sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(L,N)+1i*randn(L,N));
n2=sqrt(sigmansq)*(randn(L,N)+1i*randn(L,N));

x1=A1*S + n1;
x2=A2*S + n2;

% DOA Matrix formulation

Ra=(x1*x1')/N;
Rb=(x2*x1')/N;
nv=eigs(Ra,L);
sigmaest=mean(nv(K+1:L));
Ra=Ra-sigmaest*eye(L);
% vectorization
va=Ra(:);
vb=Rb(:);
% Removing redundancies
unqidx=getIndexOfUniqueElements(L1,L2);
vab=va(unqidx);
vbb=vb(unqidx);
% overlapping subvectors
Lb=L2*(L1+1);
Va=hankel(vab(1:Lb),vab(Lb:2*Lb-1));
Vb=hankel(vbb(1:Lb),vbb(Lb:2*Lb-1));

D=Vb*pinv(Va);
% Estimation
[F,Psi]=eigs(D,K);
uk=getUk(F,K,Lb,lambda,dx);
vk=angle(diag(Psi)).*(lambda/(2*pi*dy));

thetak=asind(sqrt(uk.^2 + vk.^2))
phik=(angle(uk+1i*vk))*180/pi


function uk=getUk(F,K,Lb,lambda,dx)
uk=zeros(K,1);
for k=1:K
    uk(k)=mean(angle((F(2:Lb))./(F(1:Lb-1))));
end
uk=uk*(lambda/(2*pi*dx));
end
function uidx=getIndexOfUniqueElements(L1,L2)
narray=[0:L1 ((2:L2).*(L1+1)-1)];
cr=[];
for i=1:length(narray)
    cr=[cr -narray+narray(i)];
end
[~,uidx,~]=unique(cr);
end

function A=getManifoldMatrixA1(theta,phi,L1,L2,dy,lambda)
A=exp(1i*(2*pi*dy/lambda)*[0:L1 ((2:L2).*(L1+1)-1)]'*(sind(theta').*sind(phi')));
end

function Phi=getPhi(theta,phi,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cosd(theta').*sind(phi'))));
end