function [meas_wBIB,qr1] = Gyro_model(Ts,true_fBIB,true_wBIB,IMU,qr0)
%% Simulates an inertial measurement unit model
%--------------------------------------------------------------------------
% Inputs:
%   Ts             time interval between epochs (s)
%   true_fBIB      true specific force of Body w.r.t. ECI expressed in the Body frame (m/s^2)
%   true_wBIB      true angular rate of Body w.r.t. ECI expressed in the Body frame (rad/s)
%   IMU            error Sources of IMU Sensors 
%   qr0            residuals of previous output quantization process
%--------------------------------------------------------------------------
% Outputs:
%   meas_wBIB      measured angular rate of Body w.r.t. ECI expressed in the Body frame (rad/s)
%   qr1            residuals of output quantization process
%==========================================================================
%% Generate noise
if Ts > 0
    gyro_noise  = randn(3,1) * IMU.gyro_ARW / sqrt(Ts);  
else
    gyro_noise  = [0;0;0];
end 
%==========================================================================
%% Calculate gyroscope outputs 
%gyroscope output
uq_wBIB  = IMU.b_g + (eye(3) + IMU.SM_g) * true_wBIB + IMU.G_g * true_fBIB + gyro_noise;
%==========================================================================    
%% Quantize gyro outputs
%Gyro quantization level (rad/s)
q_g = IMU.ql_g;
if q_g>0
    meas_wBIB = q_g * round((uq_wBIB + qr0(1:3)) / q_g);
    qr1(4:6,1) = uq_wBIB + qr0(1:3) - meas_wBIB;
else
    meas_wBIB = uq_wBIB;
    qr1(4:6,1) = [0;0;0];
end
%==========================================================================
end
