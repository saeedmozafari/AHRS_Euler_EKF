function [IMU] = Select_IMU(IMU_Type)
%% Assign IMU Specifications
%   IMU
%      .pout                   pure out of IMU 
%      .b_a                    Accelerometer bias (m/s^2)
%      .b_g                    Gyro biases (rad/s)
%      .b_m                    Magnetometer bias (nT)
%      .SM_a                   Accelerometer scale factor and cross coupling errors
%      .SM_g                   Gyro scale factor and cross coupling errors
%      .SM_m                   Magnetometer scale factor and cross coupling errors
%      .G_g                    Gyro g-dependent biases (rad-sec/m)
%      .accel_VRW              Accelerometer noise PSD (m s^-1.5)
%      .gyro_ARW               Gyro noise PSD (rad s^-0.5)
%      .mag_RW                 Magnetometer noise PSD (nT)
%      .ql_a                   Accelerometer quantization level (m/s^2)
%      .ql_g                   Gyro quantization level (rad/s)
%      .ql_m                   Magnetometer quantization level (nT)
%==========================================================================
%Constant Parameters
D2R = 0.01745329252;%         convert degree to radian
R2D = 1/D2R;%                 convert radian to degree
Mug2mps2 = 9.80665E-6;%       convert micro-g to meter per second.^2
%==========================================================================
%pure output of IMU (acceleration & angular rate)
IMU.No=IMU_Type;
% IMU.pout = IMU_pout;
if IMU.No==1
    %Accelerometers
    IMU.b_a = [9000; -13000; 8000] * Mug2mps2;
    IMU.SM_a = [ 50000,  -15000,  10000;...
                 -7500,  -60000,  12500;...
                -12500,    5000,  20000] * 1E-6;
    IMU.accel_VRW = 1000 * Mug2mps2;
    IMU.ql_a = 1E-1;
    %Gyroscopes 
    IMU.b_g = [-180; 260; -160] * D2R / 3600;
    IMU.SM_g = [40000,  -14000,   12500;...
                    0,  -30000,   -7500;...
                    0,       0,  -17500] * 1E-6;
    IMU.G_g = [90, -110,  -60;...
              -50,  190, -160;...
               30,  110, -130] * D2R / (3600 * 9.80665);
    IMU.gyro_ARW = 1 * D2R / 60;
    IMU.ql_g = 1E-2;
    %Magnetometers
    IMU.b_m = [3300;3300;3300]; 
    sm = 2E-2; %scale factor 
    m_m = 1 * D2R /3600; %input axis misalignment
    IMU.SM_m =  [ sm,  m_m,   m_m;...
                 m_m,   sm,   m_m;...
                 m_m,   m_m,  sm] ;
    IMU.mag_RW = 370; 
    IMU.ql_m = 1E2;
%--------------------------------------------------------------------------
elseif IMU.No==2
    %Accelerometers
    IMU.b_a = [900;-1300;800] * Mug2mps2;
    IMU.SM_a = [ 500,   -300,  200;...
                -150,  -600,  250;...
                -250,   100,  450] * 1E-6;
    IMU.accel_VRW = 100 *Mug2mps2;
    IMU.ql_a = 1E-2;
    %Gyroscopes
    IMU.b_g = [-9;13;-8] * D2R / 3600;
    IMU.SM_g = [400, -300,  250;...
                  0, -300, -150;...
                  0,    0, -350] * 1E-6;
    IMU.G_g = [0.9, -1.1, -0.6;...
              -0.5,  1.9, -1.6;...
               0.3,  1.1, -1.3] * D2R / (3600 * 9.80665);
    IMU.gyro_ARW = 0.01 * D2R / 60;
    IMU.ql_g = 1E-3;
    %Magnetometers
    IMU.b_m = [50; 50; 50]; 
    sm = 0.09E-2; %scale factor 
    m_m = 0.1 * D2R /3600; %input axis misalignment
    IMU.SM_m =  [ sm,  m_m,   m_m;...
                 m_m,   sm,   m_m;...
                 m_m,   m_m,  sm] ;
    IMU.mag_RW = 1; 
    IMU.ql_m = 1E1;
%--------------------------------------------------------------------------    
elseif IMU.No==3
    %Accelerometers
    IMU.b_a = [223; 223; 223] * Mug2mps2;
    sa = 223 * 1E-6; %scale factor 
    m_a = 22 * D2R /3600; %input axis misalignment
    IMU.SM_a =  [ sa,  m_a,   m_a;...
                 m_a,   sa,   m_a;...
                 m_a,   m_a,  sa] ;
    IMU.accel_VRW = 56 * Mug2mps2;
    IMU.ql_a = 1E-3;
    %Gyroscopes
    IMU.b_g = [0.11; 0.11; 0.11] * D2R / 3600;
    sg = 112 * 1E-6; %scale factor
    m_g = 22 * D2R / 3600; %input axis misalignment
    IMU.SM_g =  [ sg,  m_g,   m_g;...
                 m_g,   sg,   m_g;...
                 m_g,   m_g,  sg] ;
    IMU.G_g = zeros(3);
    IMU.gyro_ARW = 0.078 * D2R / 60;
    IMU.ql_g = 1E-5;
    %Magnetometers
    IMU.b_m = [0.2;0.2;0.2]; 
    sm = 0.02E-2; %scale factor 
    m_m = 0.05 * D2R /3600; %input axis misalignment
    IMU.SM_m =  [ sm,  m_m,   m_m;...
                 m_m,   sm,   m_m;...
                 m_m,   m_m,  sm] ;
    IMU.mag_RW = 0.3; 
    IMU.ql_m = 1E-2;
%--------------------------------------------------------------------------    
elseif IMU.No==4
    %Accelerometers
    IMU.b_a = [30;-45;26] * Mug2mps2;
    IMU.SM_a = [ 100,  -120,   80;...
                 -60,  -120,  100;...
                -100,    40,   90] * 1E-6;
    IMU.accel_VRW = 20 * Mug2mps2;
    IMU.ql_a = 1E-5;
    %Gyroscopes
        IMU.b_g = [-0.0009;0.0013;-0.0008] * D2R / 3600;
    IMU.SM_g = [8, -120, 100;...
                0,   -6, -60;...
                0,    0,  -7] * 1E-6;
    IMU.G_g = zeros(3);
    IMU.gyro_ARW = 0.002 * D2R / 60;
    IMU.ql_g = 1E-6;
    %Magnetometers
    IMU.b_m = [0.1;0.1;0.1]; 
    sm = 0.01E-2; %scale factor 
    m_m = 0.005 * D2R /3600; %input axis misalignment
    IMU.SM_m =  [ sm,  m_m,   m_m;...
                 m_m,   sm,   m_m;...
                 m_m,   m_m,  sm] ;
    IMU.mag_RW = 0.01; 
    IMU.ql_m = 1E-3;
end
end
%==========================================================================