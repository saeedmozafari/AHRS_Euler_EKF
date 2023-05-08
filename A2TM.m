function     T=A2TM(psi,theta,phi)
% A2TM       := convert euler Angles to Transformation Matrix
%==========================================================================
% T     (3*3):= Transformation matrix constructed from euler angles
% psi   (1,1):= psi Angle
% theta (1,1):= theta Angle
% phi   (1,1):= phi Angle
%============================================================================================================================
T=[ cos(psi)*cos(theta),                                                         sin(psi)*cos(theta),           -sin(theta);
    cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi),   sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi),   cos(theta)*sin(phi);
    cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi),   sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi),   cos(theta)*cos(phi)];   
%============================================================================================================================
end