% 2 level Nested array
fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=3;
M2=3;
d=lambda/2;
d1=lambda/2;
d2=(M1+1)*d1;
SNR=20;
angle=[10 40 80; 10 30 40].*pi/180;
A1=getManifoldMatrixA1(angle,M1,M2,d1,lambda);
Phi=getPhi(angle,d,lambda);
A2=A1*Phi;

B1 = kr(conj(A1),A1) %khati roa product redudancy occurs here
% P to be declared here

function A=getManifoldMatrixA1(angle,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sin(angle(1,:)).*sin(angle(2,:))));
end

function Phi=getPhi(angle,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cos(angle(1,:)).*sin(angle(2,:)))));
end

function AB = kr(A,B);
%KR Khatri-Rao product
%
% The Khatri - Rao product
% For two matrices with similar column dimension the khatri-Rao product
% is kr(A,B) = [kron(A(:,1),B(:,1)) .... kron(A(:,F),B(:,F))]
% 
% I/O AB = kr(A,B);
%
% kr(A,B) equals ppp(B,A) - where ppp is the triple-P product = 
% the parallel proportional profiles product which was originally 
% suggested in Bro, Ph.D. thesis, 1998

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.01 $ May 2001 $ Error in helpfile - A and B reversed $ RB $ Not compiled $

disp('KR.M is obsolete and will be removed in future versions. ')
disp('use KRB.M instead.')


[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in kr.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end
end