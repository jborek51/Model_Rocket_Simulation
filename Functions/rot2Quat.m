function q = rot2Quat(R)
%   Converts the body-ground rotation matrix into quaternions
q = zeros(4,1);

% Extracting elements
R11 = R(1,1);   R12 = R(1,2);   R13 = R(1,3);
R21 = R(2,1);   R22 = R(2,2);   R23 = R(2,3);
R31 = R(3,1);   R32 = R(3,2);   R33 = R(3,3);

% Calculating squares elements
q0 = (1 + R11 + R22 + R33)/4;
q1 = (1 + R11 - R22 - R33)/4;
q2 = (1 - R11 + R22 - R33)/4;
q3 = (1 - R11 - R22 + R33)/4;

% Choosing the max of the squares
[~,qIdx] = max([q0 q1 q2 q3]);

switch qIdx-1
    case 0
        q0 = sqrt(q0);        q1 = (R32-R23)/(4*q0);        q2 = (R13-R31)/(4*q0);        q3 = (R21-R12)/(4*q0);
    case 1
        q1 = sqrt(q1);        q0 = (R32-R23)/(4*q1);        q2 = (R21+R12)/(4*q1);        q3 = (R13+R31)/(4*q1);
    case 2
        q2 = sqrt(q2);        q0 = (R13-R31)/(4*q2);        q1 = (R21+R12)/(4*q2);        q3 = (R32+R23)/(4*q2);
    case 3
        q3 = sqrt(q3);        q0 = (R21-R12)/(4*q3);        q1 = (R13+R31)/(4*q3);        q2 = (R32+R23)/(4*q3);
end
q(1) = q0;  q(2) = q1;  q(3) = q2;  q(4) = q3;
end

