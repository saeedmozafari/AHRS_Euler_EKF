function [meas_fBIB,qr1] = acc_model(Ts,true_fBIB,IMU,qr0)
%% Simulates an Accelerometer unit model
%--------------------------------------------------------------------------
% Inputs:
%   Ts             time interval between epochs (s)
%   true_fBIB      true specific force of Body w.r.t. ECI expressed in the Body frame (m/s^2)
%   IMU            error Sources of IMU Sensors 
%   qr0            residuals of previous output quantization process
%--------------------------------------------------------------------------
% Outputs:
%   meas_fBIB      measured specific force of Body w.r.t. ECI expressed in the Body frame (m/s^2)
%   qr1            residuals of output quantization process
%==========================================================================
%% Generate noise
if Ts > 0
    accel_noise = randn(3,1) * IMU.accel_VRW / sqrt(Ts);  
else
    accel_noise = [0;0;0];
end 
%==========================================================================
%% Calculate accelerometer outputs 
%accelerometer output
uq_fBIB  = IMU.b_a + (eye(3) + IMU.SM_a) * true_fBIB + accel_noise;
%==========================================================================   
%% Quantize accelerometer outputs
%Accelerometer quantization level (m/s^2)
q_a = IMU.ql_a;
if q_a>0
    meas_fBIB = q_a * round((uq_fBIB + qr0(1:3)) / q_a);
    qr1(1:3,1) = uq_fBIB + qr0(1:3) - meas_fBIB;
else
    meas_fBIB = uq_fBIB;
    qr1(1:3,1) = [0;0;0];
end 
end
