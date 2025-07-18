function costVal = pinRecoveryCostFunction(x_pro_bone, Tinit_PRO_BONE, Tinit_DIS_BONE, Tfinal_DIS_PRO)
    % 1) Build the final T_PRO_BONE
    Tfinal_PRO_BONE = param2tform(x_pro_bone);

    % 2) Build the final T_DIS_BONE
    Tfinal_DIS_BONE = Tfinal_PRO_BONE * Tfinal_DIS_PRO;

    % 3) Compute distances
    % For example, measure translation error (Euclidean distance)
    % and rotation error (e.g., angle difference) for each transform.

    dist_PRO = transformDistance(Tfinal_PRO_BONE, Tinit_PRO_BONE);
    dist_DIS = transformDistance(Tfinal_DIS_BONE, Tinit_DIS_BONE);

    % 4) Sum up or weigh them as you like
    costVal = dist_PRO + dist_DIS;
end
