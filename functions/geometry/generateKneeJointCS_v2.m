function [T, r_rad, t_mm] = generateKneeJointCS_v2(T_femur_global, T_tibia_global, knee_side)
% GENERATEKNEEJOINTCS generate the 6 degree of freedom of knee joint from
% two 4x4 rigid body transformation of femur and tibia. This function is my
% implementation of Grood and Sunday paper. Some misinterpretation might
% occur, you can edit it, please be more understanding. okey?
%
% v2 is basically similar with v1, but in my case, i got T_femur_global and
% T_tibia_global from the script from ERC project. This CS has:
% - x-axis as anterior-posterior, 
% - y-axis as proximal-distal, 
% - z-axis as medial-lateral.
% In the Grood and Sunday's paper, they describe the axis differently:
% - x-axis as medial-lateral, 
% - y-axis as anterior-posterior, 
% - z-axis as medial-lateral.
% 
% Here in v2, i adjust the axis accordingly

    % define body fixed axis e1
    e_1_fix = T_femur_global(1:3, 3); % I (X in GS paper, Z in ERC)
    e_1_ref = T_femur_global(1:3, 1); % J (Y in GS paper, X in ERC)
    e_1_flo = T_femur_global(1:3, 2); % K (Z in GS paper, Y in ERC), or, cross(e_1_fix, e_1_ref);

    % define body fixed axis e3
    e_3_fix = T_tibia_global(1:3, 2); % k (z in GS paper, Y in ERC)
    e_3_ref = T_tibia_global(1:3, 1); % j (y in GS paper, X in ERC)
    e_3_flo = T_tibia_global(1:3, 3); % i (x in GS paper, Z in ERC), or, cross(e_3_fix, e_3_ref);

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

