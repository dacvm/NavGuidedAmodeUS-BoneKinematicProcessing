function [t_err, r_err] = calculateRBerror(T1, T2)
    % calculateRBerror - Computes the relative translation and rotation error 
    % between two transformation matrices.
    %
    % Syntax:
    %   [t_err, r_err] = calculateRBerror(T1, T2)
    %
    % Inputs:
    %   T1 - 4x4 homogeneous transformation matrix (ground truth or reference)
    %   T2 - 4x4 homogeneous transformation matrix (estimated or observed)
    %
    % Outputs:
    %   t_err - 3x1 vector representing translation error [dx; dy; dz]
    %   r_err - 1x3 vector representing rotation error in Euler angles [rx, ry, rz] (in degrees)
    %
    % Description:
    %   The function calculates the rigid body transformation error by finding
    %   the relative transformation from T1 to T2. It then extracts the 
    %   translational difference and computes the rotational difference in 
    %   Euler angles (XYZ order) converted to degrees.
    
    % Compute the relative transformation: T2 relative to T1
    T_boneEst_boneGT = inv(T1) * T2;

    % Extract translation error (x, y, z components)
    t_err = T_boneEst_boneGT(1:3, 4)';

    % Extract rotation error in Euler angles (XYZ order) and convert to degrees
    r_err = rad2deg(rotm2eul(T_boneEst_boneGT(1:3, 1:3), 'XYZ'));
end

