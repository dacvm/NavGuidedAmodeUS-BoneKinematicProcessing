function T = param2tform(param)
    % param(1:3) -> translation vector
    % param(4:6) -> rotation parameters (Euler angles, for example)

    % 1) Extract translation and Euler angles
    tx = param(1); ty = param(2); tz = param(3);
    rx = param(4); ry = param(5); rz = param(6);

    % 2) Build rotation from Euler angles
    % Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
    % Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
    % Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
    % R  = Rz*Ry*Rx; % or in some chosen convention
    R = eul2rotm([rx, ry, rz], "ZYX");
    t = [tx, ty, tz]';

    % 3) Form the 4x4 homogeneous transform
    T = [ R, t; 0, 0, 0, 1 ];
end
