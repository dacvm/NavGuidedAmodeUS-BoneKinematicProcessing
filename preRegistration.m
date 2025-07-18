function [T_final, boneCTpoints_afterPR2] = preRegistration(boneCTpoints_original, boneQUSpoints_ref, boneQUSlabels, bonePrereg, is_display)

if (is_display)
    fig1 = figure('Name', 'Pre-registration');
    ax1  = axes(fig1);
    axis(ax1, 'equal');
    hold(ax1, 'on');
    grid(ax1, 'on');
    xlabel(ax1, 'X');
    ylabel(ax1, 'Y');
    zlabel(ax1, 'Z');
    view(ax1, 30,15);
end

fprintf('Coarse Registration... '); tic

% PART 1: Translating with mean -------------------------------------------
% For pre registration with pre-defined area, it seems like it works
% better if i translate my boneCTpoints_ref closer to boneUSpoints_ref.
%
% So i will compute the mean of those two point clouds, and translate 
% boneCTpoints_ref to boneUSpoints_ref. But I don't like to
% just "subtract" the boneCTpoints_ref with the difference of two point
% cloud means. I want to do it elegantly. I defined the mean of the two
% point clouds as a rigid body transformation with rotation defined as
% identity and translation defined as mean.
%
% By doing this, i am generalizing the procedure of transformation,
% similar to the other registration step.

% setup the T for boneCTpoints_ref
t_meanBoneCT_ref = mean(boneCTpoints_original(1:3, :), 2);
T_meanBoneCT_ref = [eye(3), t_meanBoneCT_ref; 0 0 0 1];

% setup the T for boneUSpoints_ref
t_meanBoneUS_ref = mean(boneQUSpoints_ref(1:3, :), 2);
T_meanBoneUS_ref = [eye(3), t_meanBoneUS_ref; 0 0 0 1];

% transform
T_meanBoneUS_meanBoneCT = T_meanBoneUS_ref * inv(T_meanBoneCT_ref);
boneCTpoints_afterPR1   = T_meanBoneUS_meanBoneCT * boneCTpoints_original;

if (is_display)
    % display the amode 3d point
    scatter3( ax1, ...
              boneQUSpoints_ref(1,:), ...
              boneQUSpoints_ref(2,:), ...
              boneQUSpoints_ref(3,:), ...
              50, 'red', 'filled', ...
              'Tag', 'amode_3d');

    % show PR1
    scatter3( ax1, ...
              boneCTpoints_afterPR1(1,:), ...
              boneCTpoints_afterPR1(2,:), ...
              boneCTpoints_afterPR1(3,:), ...
              1, "magenta", "filled", ...
              "Tag", "bonepoint_prereg1");
end


% PART 2: Registration with Preregistration Area --------------------------

% Variable to populates the area centroids
allbonePRAcentroids_afterPR1 = [];
allboneUSGpoints_ref         = [];

% preregistration area contains of multiple areas, so loop over all of
% those areas
for area_idx = 1:length(bonePrereg.areas)

    % get the name of the area
    area_name   = bonePrereg.areas(area_idx).name;

    % get the preregistration area points
    bonePRApoints_ref        = bonePrereg.areas(area_idx).points;
    bonePRApoints_ref        = [bonePRApoints_ref; ones(1, size(bonePRApoints_ref,2))];
    bonePRApoints_afterPR1   = T_meanBoneUS_meanBoneCT * bonePRApoints_ref;

    % get the corresponding a-mode measurement
    group_idx         = find(ismember(boneQUSlabels, area_name));
    boneUSGpoints_ref = boneQUSpoints_ref(:, group_idx);

    % take a mean and knn (to get the centroid of the area)
    area_mean   = mean(bonePRApoints_afterPR1(1:3, :), 2);
    [nn_idx, ~] = knnsearch(bonePRApoints_afterPR1(1:3, :)', area_mean', "K", 1);
    bonePRAcentroid_afterPR1 = bonePRApoints_afterPR1(:, nn_idx);

    % since there are multiple sample point in one preregistration area,
    % and estimation algorithm requires 1-to-1 point pairs, so we
    % duplicate the centroid as the number
    allbonePRAcentroids_afterPR1 = [allbonePRAcentroids_afterPR1, repmat(bonePRAcentroid_afterPR1(1:3), 1, size(boneUSGpoints_ref, 2));];
    allboneUSGpoints_ref         = [allboneUSGpoints_ref, boneUSGpoints_ref(1:3, :)];

    if (is_display)
        % display the pre-registration area
        scatter3(  ax1, ...
                   bonePRApoints_afterPR1(1,:), ...
                   bonePRApoints_afterPR1(2,:),  ...
                   bonePRApoints_afterPR1(3,:), ...
                   1, 'b', ...
                  'filled', ...
                  'Tag', "bonepoint_prereg1");
        % display the pre-registration area centroid
        scatter3(  ax1, ...
                   bonePRAcentroid_afterPR1(1,:), ...
                   bonePRAcentroid_afterPR1(2,:),  ...
                   bonePRAcentroid_afterPR1(3,:), ...
                   20, 'k', ...
                  'filled', ...
                  'Tag', "bonepoint_prereg1");
    end
end


if (is_display)
    % for all the corresponding points...
    for i=1:size(allbonePRAcentroids_afterPR1,2)
        % make the point pair
        tmp_pair = [allbonePRAcentroids_afterPR1(:, i), allboneUSGpoints_ref(:,i)];
        % display the pair with a line
        plot3( ax1, tmp_pair(1,:), tmp_pair(2,:), tmp_pair(3,:), '-m', 'Tag', 'bonepoint_prereg1');
    end
end

% estimate the transformation using linear regressiong
tform = estgeotform3d(allboneUSGpoints_ref', allbonePRAcentroids_afterPR1', "rigid", 'MaxDistance', 100);

% in the estimation function above, i use pc1 as the boneUS and pc2 as
% the prereg centroid. The pair matching is more reliable this way. But
% it means i need to inverse the transformation, because we actually
% want from prereg centroid to bone us points
T_prereg_meanBoneUS     = inv(tform.A);

% propagate the transformation
boneCTpoints_afterPR2  = T_prereg_meanBoneUS * T_meanBoneUS_meanBoneCT * boneCTpoints_original;
T_final = T_prereg_meanBoneUS * T_meanBoneUS_meanBoneCT;

if (is_display)
% plot
scatter3( ax1, ...
          boneCTpoints_afterPR2(1,:), ...
          boneCTpoints_afterPR2(2,:),  ...
          boneCTpoints_afterPR2(3,:), ...
          1, [0.9290 0.6940 0.1250], ...
          'filled', ...
          'Tag', "bonepoint_prereg2");
end

% -------------------------------------------------------------------------

% display duration
time = toc;
fprintf('(%.2fs)\n', time);

end

