close all
clear
clc
format compact
format short g
%==========================================================================
%Choose scenario
scenario=3;
% scenario = 1      Stationary        (60 seconds)
% scenario = 2      Ground path       (175 seconds)
% scenario = 3      Arial path        (220 seconds)
if scenario==1
    load('path1.mat')
end
if scenario==2
    load('path2.mat')
end
if scenario==3
    load('path3.mat')
end
%==========================================================================
% Constants
D2R=pi/180;
R2D=180/pi;
g=9.79;
%==========================================================================
ST=.01;
TF=(size(true,1)*ST)-0.01;
%==========================================================================
%pick up the IMU.
IMU_Type = 3;
%                   Type      Grade
% IMU_Type = 1        MEMS      consumer
% IMU_Type = 2        MEMS      low-end tactical
% IMU_Type = 3        FOG       navigation
% IMU_Type = 4        FOG/RLG   strategic
[IMU] = Select_IMU(IMU_Type);
%==========================================================================
GPSF=10;    % GPS Frequency (Hz)
GPSU=(1/ST)/GPSF;   % GPS Update
GPS_vBELkm1=zeros(3,1);
GPS_vBELk=zeros(3,1);
%==========================================================================
Error0=[5;5;5]*D2R;
Est.xhat0p=True.EULER(:,1)+Error0;
% dx=         dpsi  dtheta  dphi
Est.P0p=diag([7e-3   7e-3   7e-3]);     %(5*D2R)^2~=7e-3
Est.Qc=diag ([1e-6   1e-6   1e-6]);     % ???
Est.Rk=diag ([1e-3   3e-4   3e-4]);     %(1*D2R)^2~=3e-4 & (3*D2R)^2~=3e-3
%==================================
Est.xhatkm1p=Est.xhat0p;	% x_hat_k_minus_1_positive
Est.Pkm1p=Est.P0p;	        % P_k_minus_1_positive
%==================================
Est.Xhatp=Est.xhat0p;
Est.Pp=diag(Est.P0p);

%quantization error
qra = zeros(3,1);
qrg = zeros(3,1);
qrm = zeros(3,1);

% buffer the last two records of heading
psi_buffer=[True.EULER(1,1)+Error0(1); True.EULER(1,1)+Error0(1)];
Y=[];
%==========================================================================
for k=1:1:TF/ST       % k=0 means t=0 & k=1 means t=0.01 & ...
    %======================================================================
    %################################AHRS#################################%
    %======================================================================
    [meas_fBIB,qra] = acc_model(ST,True.ABLB_ng(:,k+1),IMU,qra);
    axk=meas_fBIB(1);    % True.ABLB_ng(1,2) means ax1
    ayk=meas_fBIB(2);    % True.ABLB_ng(2,2) means ay1
    azk=meas_fBIB(3);    % True.ABLB_ng(3,2) means az1
    
    [meas_wBIB,qrg] = Gyro_model(ST,True.ABLB_ng(:,k),True.WBLB(:,k),IMU,qrg);
    pkm1=meas_wBIB(1);        % True.WBLB(1,1) means p0
    qkm1=meas_wBIB(2);        % True.WBLB(2,1) means q0
    rkm1=meas_wBIB(3);        % True.WBLB(3,1) means r0
    
    wBLBkm1=[pkm1;qkm1;rkm1];
    
    [meas_wBIB,qrg] = Gyro_model(ST,True.ABLB_ng(:,k+1),True.WBLB(:,k+1),IMU,qrg);
    pk=meas_wBIB(1);        % True.WBLB(1,2) means p1
    qk=meas_wBIB(2);        % True.WBLB(2,2) means q1
    rk=meas_wBIB(3);        % True.WBLB(3,2) means r1
    
    wBLBk=[pk;qk;rk];
    %======================================================================
    if mod(k,GPSU)==0
        if k/GPSU==1
            GPS_vBELk=GPS_model(True.VBLL(:,k));
            GPS_vBELkm1=GPS_vBELk;
        else
            GPS_vBELk = GPS_model(True.VBLL(:,k));
        end
    end
    
    TBL=A2TM(Est.xhatkm1p(1),Est.xhatkm1p(2),Est.xhatkm1p(3));
    VBLB=TBL * GPS_vBELkm1;
    ukm1=VBLB(1);        % means u0
    vkm1=VBLB(2);        % means v0
    wkm1=VBLB(3);        % means w0
    
    prop_one_step=Euler_Propag(ST,Est.xhatkm1p,wBLBk); % Propagate one step  ????????? wBLBk or wBLBkm1???
    TBL=A2TM(prop_one_step(1),prop_one_step(2),prop_one_step(3));
    VBLB=TBL * GPS_vBELk;
    uk=VBLB(1);        %  means u1
    vk=VBLB(2);        %  means v1
    wk=VBLB(3);        %  means w1
    %======================================================================
    %     udk=(uk-ukm1)/ST;
    %     vdk=(vk-vkm1)/ST;
    %     wdk=(wk-wkm1)/ST;
    
    udk=0;
    vdk=0;
    wdk=0;
    %======================================================================
    % theta_m=asin(ax/GEMT.g);
    % phi_m=asin(ay/(-GEMT.g*cos(theta_m)));
    % phi_m=atan(ay/az/GEMT.g);
    % Measurments=[Measurments [0;theta_m;phi_m]];
    %======================================================================
    theta_m=asin((udk-vk*rk+wk*qk-axk)/-g);
    %     phi_m=atan((vdk+uk*rk-wk*pk-ayk)/(wdk-uk*qk+vk*pk-azk));
    phi_m=asin((vdk+uk*rk-wk*pk-ayk)/(g*cos(theta_m)));
    %======================================================================
    %     M_L=wrldmagm(-True.SBLL(3,k),35+True.SBLL(1,k)/Re,55+True.SBLL(2,k)/Re,2011,'2010');
    M_L = igrfmagm(True.SBLL(3,k),True.SBLL(1,k)* R2D,True.SBLL(2,k)* R2D,2011)'; %magnetometer(nT)
    D=atan2(M_L(2),M_L(1));
    TBL=A2TM(True.EULER(1,k),True.EULER(2,k),True.EULER(3,k));
    M_B=TBL*M_L;
    mx=M_B(1);
    my=M_B(2);
    mz=M_B(3);
    dpsi_m=atan2((-cos(phi_m)*my+sin(phi_m)*mz),(cos(theta_m)*mx+sin(theta_m)*sin(phi_m)*my+sin(theta_m)*cos(phi_m)*mz));
    
    %heading signal phase unwraping
    psi_buffer=[psi_buffer(2);dpsi_m+D];
    psi_buffer = Phase_Unwrap(psi_buffer,pi);
    psi_m=psi_buffer(2);
    % D will be read from lookup-table
    %======================================================================
    yk=[psi_m;theta_m;phi_m];
    Y=[Y yk];
    %======================================================================
    [Est.xhatkp,Est.Pkp]=ExtKalFil(wBLBkm1,ST,Est.xhatkm1p,Est.Pkm1p,yk,Est.Qc,Est.Rk);
    %======
    % euler
    Est.xhatkm1p = Est.xhatkp;  %for next iteration
    Est.Pkm1p=Est.Pkp;          %for next iteration
    %======
    Est.Xhatp=[Est.Xhatp Est.xhatkp];
    Est.Pp=[Est.Pp diag(Est.Pkp)];
end
%================================================================================================================================================================
% t=ST:ST:TF;
% figure;plot(t,True.VBLL(1,:));xlabel('time (s)');ylabel('v_N (m/s)');
% figure;plot(t,True.VBLL(2,:));xlabel('time (s)');ylabel('v_E (m/s)');
% figure;plot(t,True.VBLL(3,:));xlabel('time (s)');ylabel('v_D (m/s)');
% figure;plot(t,True.SBLL(1,:));xlabel('time (s)');ylabel('s_N (m)');
% figure;plot(t,True.SBLL(2,:));xlabel('time (s)');ylabel('s_E (m)');
% figure;plot(t,-True.SBLL(3,:));xlabel('time (s)');ylabel('h (m)');
% figure;plot(t,True.EULER(1,:)*R2D);xlabel('time (s)');ylabel('\psi (deg)');
% figure;plot(t,True.EULER(2,:)*R2D);xlabel('time (s)');ylabel('\theta (deg)');
% figure;plot(t,True.EULER(3,:)*R2D);xlabel('time (s)');ylabel('\phi (deg)');
%==========================================================================
% t=0:ST:TF;
% Y=[zeros(3,1) Y];
% figure;plot(t,-True.SBLL(3,:));xlabel('time (s)');ylabel('h (m)')
% figure;plot(t,True.EULER(1,:)*R2D,'r',t,Y(1,:)*R2D,'b-.',t,Est.Xhatp(1,:)*R2D,'g--');
% xlabel('time (s)');ylabel('\psi (deg)');legend('True','Measured','Estimated')
% figure;plot(t,True.EULER(2,:)*R2D,'r',t,Y(2,:)*R2D,'b-.',t,Est.Xhatp(2,:)*R2D,'g--');
% xlabel('time (s)');ylabel('\theta (deg)');legend('True','Measured','Estimated')
% figure;plot(t,True.EULER(3,:)*R2D,'r',t,Y(3,:)*R2D,'b-.',t,Est.Xhatp(3,:)*R2D,'g--');
% xlabel('time (s)');ylabel('\phi (deg)');legend('True','Measured','Estimated')
%==============================================
t=0:ST:TF;
Y=[zeros(3,1) Y];
figure;plot(t,True.SBLL(3,:));xlabel('time (s)');ylabel('h (m)')
figure;plot(t,True.EULER(1,:)*R2D,'r',t,Y(1,:)*R2D,'b-.',t,Est.Xhatp(1,:)*R2D,'g--');
xlabel('time (s)');ylabel('\psi (deg)');legend('True','Measured','Estimated')
figure;plot(t,True.EULER(2,:)*R2D,'r',t,Y(2,:)*R2D,'b-.',t,Est.Xhatp(2,:)*R2D,'g--');
xlabel('time (s)');ylabel('\theta (deg)');legend('True','Measured','Estimated')
figure;plot(t,True.EULER(3,:)*R2D,'r',t,Y(3,:)*R2D,'b-.',t,Est.Xhatp(3,:)*R2D,'g--');
xlabel('time (s)');ylabel('\phi (deg)');legend('True','Measured','Estimated')
%==============================================