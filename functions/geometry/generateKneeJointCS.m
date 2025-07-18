function [T, r_rad, t_mm] = generateKneeJointCS(T_femur_global, T_tibia_global, knee_side)
% GENERATEKNEEJOINTCS generate the 6 degree of freedom of knee joint from
% two 4x4 rigid body transformation of femur and tibia. This function is my
% implementation of Grood and Sunday paper. Some misinterpretation might
% occur, please be more understanding. okey?

    % define body fixed axis e1
    e_1_fix = T_femur_global(1:3, 1); % I
    e_1_ref = T_femur_global(1:3, 2); % J
    e_1_flo = T_femur_global(1:3, 3); % K, or, cross(e_1_fix, e_1_ref);

    % define body fixed axis e3
    e_3_fix = T_tibia_global(1:3, 3); % k
    e_3_ref = T_tibia_global(1:3, 2); % j
    e_3_flo = T_tibia_global(1:3, 1); % i, or, cross(e_3_fix, e_3_ref);

    % calculate floating axis e2
    e_2 = cross(e_1_fix, e_3_fix);

    % calculate the rotation
    alpha = asin(dot(-e_2, e_1_flo));

    if(knee_side=='r')
        beta  = acos(dot(e_1_fix, e_3_fix)) - pi/2;
    else
        beta  = pi/2 - acos(dot(e_1_fix, e_3_fix));
    end

    if(knee_side=='r')
        gamma = asin(dot(-e_2, e_3_flo));
    elseif(knee_side=='l')
        gamma = asin(dot(e_2, e_3_flo));
    end
    
    % flexion, adduction, external
    r_rad = [alpha, beta, gamma];

    % calculate the translation
    H  = T_tibia_global(1:3, 4) - T_femur_global(1:3, 4);
    q1 = dot(H, e_1_fix);
    q2 = dot(H, e_2);
    q3 = dot(-H, e_3_fix);
    % tibial thrust, tibial drawer, distraction
    t_mm = [q1, q2, q3];

    % Still don't know how to implement, i can just convert the r_rad to
    % rotm, but there must be more efficient way to do that.
    T = eye(1);

end

