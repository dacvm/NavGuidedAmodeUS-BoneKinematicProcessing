function T_avg = averageRigidTransforms(T1, T2)
%AVERAGERIGIDTRANSFORMS Averages two 4x4 rigid-body transformations.
%
%   T_avg = averageRigidTransforms(T1, T2) computes the average of the
%   transformations T1 and T2 by averaging the translation components and
%   the rotation components via quaternion averaging.
%
%   Inputs:
%     T1, T2 - 4x4 homogeneous transformation matrices.
%
%   Output:
%     T_avg  - The averaged 4x4 transformation matrix.
%
%   Note: This function uses the MATLAB functions rotm2quat and quat2rotm,
%   available in the Robotics System Toolbox.
%
%   Example:
%       T1 = eye(4);
%       T2 = [quat2rotm([0.9239 0 0.3827 0]), [1;2;3]; 0 0 0 1];
%       T_avg = averageRigidTransforms(T1, T2);

    % Extract rotation matrices and translation vectors.
    R1 = T1(1:3, 1:3);
    t1 = T1(1:3, 4);
    R2 = T2(1:3, 1:3);
    t2 = T2(1:3, 4);

    % Convert rotations to quaternions.
    % MATLAB's rotm2quat returns a 1x4 quaternion in the format [w x y z].
    q1 = rotm2quat(R1);
    q2 = rotm2quat(R2);

    % Ensure the quaternions are in the same hemisphere.
    if dot(q1, q2) < 0
        q2 = -q2;
    end

    % Average the quaternions and normalize.
    q_avg = (q1 + q2) / norm(q1 + q2);

    % Convert the averaged quaternion back to a rotation matrix.
    R_avg = quat2rotm(q_avg);

    % Average the translation vectors.
    t_avg = (t1 + t2) / 2;

    % Reconstruct the homogeneous transformation matrix.
    T_avg = eye(4);
    T_avg(1:3, 1:3) = R_avg;
    T_avg(1:3, 4) = t_avg;
end
