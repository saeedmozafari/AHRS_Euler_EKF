function [meas_VBLL] = GPS_model(True_VBLL)
meas_VBLL=True_VBLL+random('Normal',[0;0;0],[0.2;0.2;0.2]);
end