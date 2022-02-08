function eulAng = quat2Euler(q)
%   Converts quaternions to euler angles 
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

eulAng = zeros(3,1);

eulAng(1) = -atan2(....
    2*(q0*q1+q2*q3),...
    1-2*(q1^2+q2^2));

eulAng(2) = -asin(2*(q0*q2-q3*q1));

eulAng(3) = -atan2(...
    2*(q0*q3+q1*q2),...
    1-2*(q2^2+q3^2));
end

