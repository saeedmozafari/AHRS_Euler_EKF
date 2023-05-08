function [xhatkp,Pkp]=ExtKalFil(wBLBkm1,ST,xhatkm1p,Pkm1p,yk,Qc,Rk)
%=================================================================================================
% Linearized Continues Model dx' = Fc*dx + Lc*w
% dxhatkm1p=xhatkm1p - xhatkn;
%p=wBLBkm1(1);
q=wBLBkm1(2);
r=wBLBkm1(3);
phi=xhatkm1p(3);
theta=xhatkm1p(2);
Fc=[0   sin(theta)/(cos(theta))^2*(q*sin(phi)+r*cos(phi))   (q*cos(phi)-r*sin(phi))/cos(theta); ...
    0   0                                                  -(q*sin(phi)+r*cos(phi)); ...
    0   (q*sin(phi)+r*cos(phi))/(cos(theta))^2               tan(theta)*(q*cos(phi)-r*sin(phi))];
%==========================================================================
% w=[wx;wy;wz]=[dp;dq;dr]
Lc=[0   sin(phi)*sec(theta)   cos(phi)*sec(theta); ...
    0   cos(phi)             -sin(phi); ...
    1   sin(phi)*tan(theta)   cos(phi)*tan(theta)];
%=================================================================================================
% attention: this procedure can be done by the symbolic procedure to reduce
%            the computational efforts.
% discretizing the Linearized continues model
% X'(t)=A*X(t)+B*U(t)
% X(k+1)=G*X(k)+H*U(k)
% t=kT
% T is sampling period
% d is approximation index
% [Fkm1,Lkm1]=Con2Dis(Fc,Lc,ST,4);
%==========================================================================
% Discretizing the continues model to x(k)=Fkm1(k-1)*x(k-1)+Qkm1
Fkm1 = eye(3)+Fc*ST;
Qkm1_tilda=Lc*Qc*Lc'*ST;
%==========================================================================
% Time Update
% attention: if Qk=Q, i.e. Q=constant, then calculation of Lkm1*Qkm1*Lkm1.'
%            should be performed only one time in the initialization.
Pkn=Fkm1*Pkm1p*Fkm1.'+Qkm1_tilda; % Qkm1_tilda ~= Lkm1*Qkm1*Lkm1.'
xhatkn=Euler_Propag(ST,xhatkm1p,wBLBkm1);
%==========================================================================
% y(k)=H(k)*x(k)+v(k)
Hk=eye(3);
Mk=eye(3);
%==========================================================================
% Measurement Update
Kk=Pkn*Hk.'/(Hk*Pkn*Hk.'+Mk*Rk*Mk.');
xhatkp=xhatkn+Kk*(yk-Hk*xhatkn);
Pkp=(eye(3)-Kk*Hk)*Pkn;
%==========================================================================
end