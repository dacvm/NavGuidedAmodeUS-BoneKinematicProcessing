%% HELOOOO
% -- This script is intended to perform dpeth projection. I want to get all of
%    the depth ground truth and then project it back to the mmode images.
% -- This script requires geom3d for calculating the depth projection.
%    which can be installed from here:
%    https://nl.mathworks.com/matlabcentral/fileexchange/24484-geom3d

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT] window configuration file
csvfile_windowconfig = 'transducerconfig_v8a_window2024-12-20_14-37-59_edited2025-04-10_15-53-56.csv';

% [EDIT] holder configuration file
csvfile_holderconfig = 'transducerconfig_v8a.csv';

% [EDIT] Specify folder index
folder_idx = 1;

% [EDIT] Select bone and pin
idx_bone = 2;
idx_pin  = 1;

% [EDIT] Select the UST
ust_numbers = 14:26;

% [EDIT] Select beam type
beamType = 'tube';

% [EDIT]
is_display = true;
is_saveMat = false;


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% declare some of important path of the data
path_bonescan     = fullfile(path_data, dir_trial, 'bonescan');
path_intermediate = fullfile(path_data, dir_trial, 'intermediate');
path_measurement  = fullfile(path_data, dir_trial, 'measurement');
path_bonestl  = fullfile(path_root, "data", "ct", "bone");

% 1) Generate path to function directory
addpath(genpath(path_function));

% 2) Get all the folders of the intermediates
items_dir = dir(path_intermediate);
folders_intermediate = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};

% 3) Get all the folders of the measurements
items_dir = dir(path_measurement);
folders_measurement = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};
% the last folder for measurements is always empty (it was like that from
% the navigation system) we don't need that.
folders_measurement(end) = [];

% 4) Load window configuration file
csv_fullpath = fullfile(path_intermediate, csvfile_windowconfig);
ust_windowconfig = table2struct(readtable(csv_fullpath));

% 5) Load holder configuration file
csv_fullpath = fullfile(path_data, 'configs', csvfile_holderconfig);
ust_holderconfig = table2struct(readtable(csv_fullpath));
% Find indices of structs where Group is not equal to 0, means that we are
% not using the ust
indices = arrayfun(@(x) x.Group ~= 0, ust_holderconfig);
% Create a new array of structs excluding those with Group equal to 0
ust_holderconfig = ust_holderconfig(indices);



%% INITIALIZE SOME RIGID BODY METADATA

% 1) Create base rotation to trasnform qualisys base vector to matlab
% -- Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

% 2.a) Load and organize the struct for the CT data (allBone_CT)
run('extra_structCTdata.m');

% 2.b) Get the relevant point cloud of the bone from CT scan
boneCTstl_original     = allBone_CT(idx_bone).stl;
boneCTpc_original      = pointCloud(boneCTstl_original.Points);
boneCTnormals_original = STLVertexNormals(boneCTstl_original.ConnectivityList, boneCTstl_original.Points)';
boneCTpoints_original  = [boneCTpc_original.Location'; ones(1, length(boneCTpc_original.Location))];

% 2.c) Get the relevant transformation of the bone from CT scan
T_bone_ct = allBone_CT(idx_bone).T;
bone_name = allBone_CT(idx_bone).name;
T_pin_ct  = allBone_CT(idx_bone).pin(idx_pin).T;
pin_name  = allBone_CT(idx_bone).pin(idx_pin).name;

% 3) There is this thing called REF. 
% -- It is a calibration tool for B-mode ultrasound. We were using this for 
%    bone surface reconstruction as well.
% -- They were all represented in respect of this REF. So naturally we are 
%    using this REF as our reference for everything. 
% -- Well, it is actually not really necessary to do that, we can just use 
%    Qulisys coordinate frame, but because in the data processing in the 
%    earlier experiment we are using REF coordinate frame, so let's just 
%    use it again.
ref_name   = 'B_N_REF';



%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

fig1 = figure('Name', 'UST in 3d');
ax1  = axes(fig1);
axis(ax1, 'equal');
hold(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
zlabel(ax1, 'Z');
view(ax1, 30,15);



%% MAIN PROGRAM

% % Load the rigid body data that is stored in a .csv File
% % -- Use the dir function to find all files with the .csv extension in the
% %    directory. 
% % -- This .csv file is the recording of the rigid bodies from Qualisys. 
% % -- Here, we assume that there is always only one csv file the directory
% fileList = dir(fullfile(path_measurement, folders_measurement{folder_idx}, '*.csv'));
% 
% % Check if the .csv file was found
% if ~isempty(fileList)
%     % Assuming there's only one CSV file, get its full name
%     fullFileName = fullfile(path_measurement, folders_measurement{folder_idx}, fileList(1).name);
%     % Read the CSV file and get the rigid bodies data (Table object)
%     all_rigidbodies_table = readCSV_qualisysRigidBodies(fullFileName);
% else
%     % Display a message if no CSV file was found
%     disp('No CSV file found in the specified directory.');
%     return;
% end

% This one just for shortcut
load('all_rigidbodies_table.mat');

% Generate discrete points to represent the origins of the hypothetical 
% beam. This script is based on Paper 2, experiment with Maxime.
% xy_circle will be used in the depth gt calculation
radius = 3;
step   = 1;
x      = -radius : step : radius;
[X, Y] = meshgrid(x, x);
xy     = [X(:), Y(:)];
xy_circle = xy( sqrt(sum(xy.^2, 2)) <= radius, :);

% Some initial variables to control the loops
n_ust_selected  = size(ust_numbers, 2);
n_timestamp     = size(all_rigidbodies_table, 1);

% Some variable to store all the depth data
all_depthgtmean_matrix = zeros(n_timestamp, n_ust_selected);
all_depthgtstd_matrix  = zeros(n_timestamp, n_ust_selected);

% initialize waitbar, so the user not get bored
h = waitbar(0, 'Please wait...');

for idx_t = 1:n_timestamp
    %% LOOP FOR TIMESTAMPS

    % show the waitbar progress
    waitbar(idx_t/n_timestamp, h, sprintf('Timestamp processed: %d/%d', idx_t, n_timestamp));
    
    % 0) delete every object in the plot before we go
    delete(findobj('Tag', 'cs_bonepin'));
    delete(findobj('Tag', 'cs_bonegt'));
    delete(findobj('Tag', 'cs_ust'));
    delete(findobj('Tag', 'bonesurface_ct'));
    delete(findobj('Tag', 'amode_3dgt'));    

    % 1) Get the current T_pin_Q
    T_pin_Q = all_rigidbodies_table.(pin_name)(idx_t);
    T_pin_Q = T_pin_Q{1}.T;
    T_ref_Q = all_rigidbodies_table.(ref_name)(idx_t);
    T_ref_Q = T_ref_Q{1}.T;

    % 2) Calculate transformation of pin in ref
    T_pin_ref = inv(T_ref_Q) * T_pin_Q;

    % 3) Calculate transformation of boneGT in ref
    T_bone_pin    = inv(T_pin_ct) * T_bone_ct;
    T_boneGT_ref  = T_pin_ref * T_bone_pin;

    % 4) Transform the point cloud
    T_ct_bone        = inv(T_bone_ct);
    boneCTpoints_GT  = T_boneGT_ref           * inv(T_bone_ct)      * boneCTpoints_original;
    boneCTnormals_GT = T_boneGT_ref(1:3, 1:3) * T_ct_bone(1:3, 1:3) * boneCTnormals_original;
    boneCTtri_GT     = triangulation(boneCTstl_original.ConnectivityList, boneCTpoints_GT(1:3,:)');

    usts_depthgtmean = [];
    usts_depthgtstd  = [];

    % loop for all ust
    for idx_ust = ust_numbers
        %% LOOP FOR USTS

        % 1) Get the current ust and its group name -----------------------
        current_uststruct    = ust_holderconfig(idx_ust);
        current_ustnumber    = ust_holderconfig(idx_ust).Number;
        current_ustgroupname = ust_holderconfig(idx_ust).GroupName;
    

        % 2) Transform the euler to rotation matrix -----------------------
        % -- Here i used the equivalent lines as in the Qt C++ code. 
        % -- Why i do this, so that i know i am doing it exactly like i am 
        %    doing in the Qt C++ code (check VolumeAmodeVisualizer::updateTransformations)
        local_r_euler  = [current_uststruct.Local_Rx, current_uststruct.Local_Ry, current_uststruct.Local_Rz];
        local_R_matrix = rotx(local_r_euler(1)) * ...   % rotate about X by rx
                         roty(local_r_euler(2)) * ...   %   then about Y by ry
                         rotz(local_r_euler(3));        %   then about Z by rz
        % -- Get the translation vector
        local_t = [current_uststruct.Local_tx, current_uststruct.Local_ty, current_uststruct.Local_tz];
        % -- Ccreate the T matrix
        T_probe_holder = [local_R_matrix, local_t'; 0 0 0 1];


        % 3) Get the current T_holder_Q for relevant idx_ust --------------
        current_RBstruct = all_rigidbodies_table.(current_ustgroupname)(idx_t);
        T_holder_Q       = current_RBstruct{1}.T;
        % -- Transform the Probe to REF coordinate frame
        T_probe_Q        = inv(T_ref_Q) * T_holder_Q * T_probe_holder;


        % 4) Calculate distance -------------------------------------------
        % -- This script is based on Paper 2, experiment with Maxime
        % -- There are several type of the beam that we use here:
        % 1) line: A hypothetical line that span from the origin point of
        %          transducers to the bone mesh.
        % 2) tube: A hypotethical tube that span from the entire surface of the
        %          transducer to the bone mesh. For simplicity, we define the
        %          origin point as a circle that lie on the transducer surface.
        %          So basically, it is a line but they are multiple.

        % if the user specify line
        if(strcmp(beamType, 'line'))

            % 4.a) Create the 3d line from the markerstick.
            % ---- The line is represented in parametric form :
            %      [origin, direction] = [x0 y0 z0 dx dy dz]
            % ---- The translation vector of T_probe_Q can serve as origin
            % ---- The z-component (column) of the R of the T_probe_Q can serve
            %      as the direction of the line
            origins   = [T_probe_Q(1,4), T_probe_Q(2,4), T_probe_Q(3,4)];
            direction = [T_probe_Q(1,3), T_probe_Q(2,3), T_probe_Q(3,3)];
            line      = [ origins, -direction ];
    
            % Prepare some variables to store the intersection point and
            % intersection distances. Actually, since in the "line" mode, we
            % don't need to initialize the variable, you can directly assign
            % the value to a variable when you find the intersection. But, i
            % just wanted to be consitent with "tube" mode. There, since i use
            % parallel for, it is required to initialize the variables first.
            lineplane_intersectpoints = [];
            lineplane_intersectdists  = [];
    
            % 4.b) Find the intersection point
            % ---- Param required: the 3d line, the vertices (= stl points) and faces (= stl connectivity list).
            tmp_inters = intersectLineMesh3d(line, boneCTtri_GT.Points, boneCTtri_GT.ConnectivityList);
    
            % 4.c) In the case of multiple intersections... (i.e. one when the 
            % ---- line enters the mesh and one when the line leaves the mesh)
            if ~isempty(tmp_inters) 
    
                % Calculate the distances for both intersections and take the 
                % smallest distance as the distance we are interested in
                distancetointers = zeros(1,length(tmp_inters(:,1)));
                for i = 1:length(tmp_inters(:,1))
    
                    % Clculate the distance between the intersection point and 
                    % transducer tip (origin), for both intersections
                    distancetointers(i) = distancePoints3d(tmp_inters(i,:), origins); % distancePoints3d(P1, P2) returns distance between points P1 and P2
                end
    
                % Get distance ground truth data, using the smallest distance, 
                % because intersections further away do not represent the bone 
                % surface that faces the ultrasound transducer 
                [lineplane_intersectdists, ind] = min(distancetointers); 
                % collect intersection closest to transducer origin in INTERS
                lineplane_intersectpoints = (tmp_inters(ind,:));
                
            % 4.d) In the case of no intersection (probe is not facing the bone)
            else
               lineplane_intersectdists  = [];
               % define intersection point as an empty array
               lineplane_intersectpoints = [];
            end

        % if the user specify tube
        elseif(strcmp(beamType, 'tube'))
            
            % 4.a) Instead of one point, origin is represented as multiple points on
            %      the probe beam surface. We do this by making a circle point grid
            %      in global frame then transform it to transducer rigid body.
            origins    = homogeneous2cartesian( T_probe_Q * cartesian2homogeneous( [xy_circle, zeros(size(xy_circle, 1), 1)] ) );
            direction = [T_probe_Q(1,3), T_probe_Q(2,3), T_probe_Q(3,3)];

            % Prepare some variables, this is required if are using parfor
            lineplane_intersectpoints = zeros(size(origins));
            lineplane_intersectdists  = zeros(size(origins, 1), 1);

            parfor i=1:length(origins)

                % For each point from the circle point grid, the the
                % intersection line to the bone mesh, represented in parametric form : [x0 y0 z0 dx dy dz]
                line = [origins(i,:), -direction];

                % 4.b) Find the intersection point
                % ---- Param required: the 3d line, the vertices (= stl points) and faces (= stl connectivity list).
                tmp_inters = intersectLineMesh3d(line, boneCTtri_GT.Points, boneCTtri_GT.ConnectivityList);
        
                % 4.c) In the case of multiple intersections... (i.e. one when the 
                % ---- line enters the mesh and one when the line leaves the mesh)
                if ~isempty(tmp_inters) 
        
                    % Calculate the distances for both intersections and take the 
                    % smallest distance as the distance we are interested in
                    distancetointers = zeros(1,length(tmp_inters(:,1)));
                    for j = 1:length(tmp_inters(:,1))
        
                        % Clculate the distance between the intersection point and 
                        % transducer tip (origin), for both intersections
                        distancetointers(j) = distancePoints3d(tmp_inters(j,:), origins(i,:)); % distancePoints3d(P1, P2) returns distance between points P1 and P2
                    end
        
                    % Get distance ground truth data, using the smallest distance, 
                    % because intersections further away do not represent the bone 
                    % surface that faces the ultrasound transducer 
                    [lineplane_intersectdist, ind] = min(distancetointers); 
                    % collect intersection closest to transducer origin in INTERS
                    lineplane_intersectpoint = (tmp_inters(ind,:));
                    
                    % 4.d) Store the data to our variable that is declared 
                    % outside the parfor loop
                    lineplane_intersectpoints(i,:) = lineplane_intersectpoint;
                    lineplane_intersectdists(i)    = lineplane_intersectdist;
                end

            end

            % 4.e) If a point doesn't have any intersection, it will return nothing,
            % and since we avoid storing nothing to a variable (which will give 
            % us an error), the initial value (zero) is preserved. We can 
            % actually ignore them but it is annoying when we want to display 
            % the result. There will be points at zero coordinate, so let's
            % delete the delete rows with zero
            lineplane_intersectpoints(~any(lineplane_intersectpoints,2), :) = [];
            lineplane_intersectdists(lineplane_intersectdists==0) = [];
        
        end


        % 5) Display ------------------------------------------------------
        if(is_display)

            % display pin axes
            display_axis(ax1, T_pin_ref(1:3, 4), T_pin_ref(1:3, 1:3), 30, 'T_pin_ref', 'Tag', 'cs_bonepin');
            % display bone axes
            display_axis(ax1, T_boneGT_ref(1:3, 4), T_boneGT_ref(1:3, 1:3), 30, 'T_boneGT_ref', 'Tag', 'cs_bonegt');
            % display bone
            trisurf(boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'Tag', 'bonesurface_ct');

            % display ust axes
            axisname = sprintf('%d_%d', current_uststruct.Group, current_uststruct.Number);
            display_axis(ax1, T_probe_Q(1:3, 4), T_probe_Q(1:3, 1:3), 10, axisname, 'Tag', 'cs_ust');
    
            % Display intersection points
            if(~isempty(lineplane_intersectpoints))
                scatter3( ax1, ...
                          lineplane_intersectpoints(:,1), ...
                          lineplane_intersectpoints(:,2), ...
                          lineplane_intersectpoints(:,3), ...
                          10, 'red', 'filled', ...
                          'Tag', 'amode_3dgt');
            end
        end


        % 6) Store --------------------------------------------------------
        % calculate the mean and the std of the depth (relevant for tube)
        ust_depthgtmean = mean(lineplane_intersectdists);
        ust_depthgtstd  = std(lineplane_intersectdists);
        % storing temporarily the values to an array
        usts_depthgtmean = [usts_depthgtmean, ust_depthgtmean];
        usts_depthgtstd  = [usts_depthgtstd, ust_depthgtstd];
    
    % end ust loop
    end

    % Store the depth mean and std to the big matrix
    all_depthgtmean_matrix(idx_t, :) = usts_depthgtmean;
    all_depthgtstd_matrix(idx_t, :)  = usts_depthgtstd;

% end timestamp loop
end

% Close the waitbar
close(h);

% Initialize a table (with ust numbers as the column name) to organize 
% the depth data better
ust_cell     = arrayfun(@num2str, ust_numbers, 'UniformOutput', false);
table_header = [{'Timestamp_idx'}, ust_cell];

% Convert the big matrix to table
ust_timestampidx_vector = 1:n_timestamp;
all_depthgtmean_table   = array2table([ust_timestampidx_vector', all_depthgtmean_matrix], 'VariableNames', table_header);
all_depthgtstd_table    = array2table([ust_timestampidx_vector', all_depthgtstd_matrix], 'VariableNames', table_header);

% save the data
if(is_saveMat)
    % get the session name
    tmp_str = split(dir_trial, '_');
    sess_str = tmp_str{3};
    meas_str = tmp_str{4};
    % save
    mat_filename = sprintf('all_depthgt_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputs, 'output_allgtdepths', mat_filename);
    save(mat_fullpath, 'all_depthgtmean_table', 'all_depthgtstd_table');
end