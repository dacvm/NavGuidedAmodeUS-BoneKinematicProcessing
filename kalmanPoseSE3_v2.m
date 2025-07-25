function [T_optWorld, ukf] = kalmanPoseSE3(Ts_est, t, T_meas, t_new, ukf)
% Pose UKF on SE(3) with constant linear & angular acceleration
%
% Ts_est : 4x4x3   previous poses  (worldbody)
% t      : 1x3     their time stamps
% T_meas : 4x4     new measurement (worldbody)
% t_new  : scalar  time of the measurement
% ukf    : (opt.)  filter handle for looping
%
% T_optWorld : 4x4 optimal pose after update
% ukf        :    updated filter object
% -------------------------------------------------------------------------

%% 1. Initialisation (first call only) 
if nargin < 5 || isempty(ukf)
    % --- extract last three poses ---
    p = zeros(3,3);  R = zeros(3,3,3);
    for k = 1:3
        p(k,:)   = tform2trvec(Ts_est(:,:,k));          % world coords
        R(:,:,k) = tform2rotm(Ts_est(:,:,k));
    end
    q = rotm2quat(R);                                   % 34 [w x y z]

    % --- finitediff kinematics (world frame) ---
    dt12   = t(2)-t(1);  dt23 = t(3)-t(2);
    v12    = (p(2,:) - p(1,:)) / dt12;
    v23    = (p(3,:) - p(2,:)) / dt23;
    a_lin  = (v23     - v12  ) / dt23;

    w12    = rotDiff2omega(R(:,:,1),R(:,:,2),dt12);     % world frame
    w23    = rotDiff2omega(R(:,:,2),R(:,:,3),dt23);
    a_ang  = (w23     - w12  ) / dt23;

    x0 = [p(3,:)'; v23'; a_lin'; q(3,:)'; w23'; a_ang'];   % 191 column

    % --- covariances & filter object ---
    P0 = blkdiag(1.0*eye(3), ...
                 1.0*eye(3), ...
                 1.0*eye(3), ...
                 1.0*eye(4), ...
                 1.0*eye(3), ...
                 1.0*eye(3));

    % Process Noise
    Q  = blkdiag(1.0*eye(3), ...
                 2.0*eye(3), ...
                 3.0*eye(3), ...
                 1.0*eye(4), ...
                 2.0*eye(3), ...
                 3.0*eye(3));

    % Measurement Noise
    sigma_pos = 1.5;            % mm
    sigma_ang = deg2rad(2.0);   % degree
    Rm = blkdiag(sigma_pos.^2 * eye(3), ...
                 sigma_ang.^2 * eye(4));

    ukf = unscentedKalmanFilter( ...
        @fState, @fMeas, x0, ...
        'StateCovariance',P0, ...
        'ProcessNoise',Q, ...
        'MeasurementNoise',Rm, ...
        'Alpha',0.6,'Beta',2,'Kappa',0);
end

%% 2. PREDICT ------------------------------------------------------------------------
dt = t_new - t(end);
if dt <= 0
    warning('kalmanPoseSE3:NonPositiveDt','Skipped predict: nonpositive dt.');
else
    predict(ukf, dt);                                  % uses dt in fState
end

%% --- 3.  CORRECT (robust version) ---------------------------------------
% 3a. Re-orthonormalise the measured rotation and build z
R_meas = T_meas(1:3,1:3);
t_meas = T_meas(1:3,4);

% Orthonormalise with SVD    R_meas  SO(3)
[U,~,V] = svd(R_meas);
R_meas  = U*V.';
if det(R_meas) < 0                   % make sure det = +1 (no reflection)
    U(:,3) = -U(:,3);                % flip third column, then rebuild
    R_meas = U*V.';
end

qMeas = rotm2quat(R_meas).';         % 41 column  [w x y z]

% Align quaternion hemisphere with prediction
qPred = ukf.State(10:13);            % 41 column
if dot(qPred,qMeas) < 0, qMeas = -qMeas; end

z = [t_meas; qMeas];                 % 71 measurement vector

% 3b. Validate & update
if all(isfinite(z))
    correct(ukf, z);                 % <- UKF update
    ukf.State(10:13) = quatnormalize(ukf.State(10:13).').'; % stay on S
else
    warning('kalmanPoseSE3:BadMeasurement', ...
            'Non-finite measurement detected  update skipped, prediction kept.');
end

%% 4. Return updated pose -------------------------------------------------------------
pOpt = ukf.State(1:3).';                               % 13 row
qOpt = ukf.State(10:13).';                             % 14 row
T_optWorld = trvec2tform(pOpt) * quat2tform(qOpt);     % [R p]  pose

end



% -------------------------------------------------------------------------
function xNext = fState(x,dt)
% Constantacceleration motion model (world frame)

% unpack
p=x(1:3); v=x(4:6); a=x(7:9);
q=x(10:13); w=x(14:16); al=x(17:19);

% linear velocity
vNext = v + a*dt;
% position
pNext = p + v*dt + 0.5*a*dt^2;

% angular velocity
wNext = w + al*dt;

% quaternion
wAvg  = 0.5*(w + wNext);
if norm(wAvg)*dt < 1e-12
    dqRow = [1 0 0 0];
else
    axisRow = (wAvg/norm(wAvg)).';              % 13 row
    dqRow   = axang2quat([axisRow norm(wAvg)*dt]); % 14 row
end
qRow  = q.';                                     % 14
qNext = quatnormalize( quatmultiply(dqRow, qRow) ).'; % FIX premultiply

xNext = [pNext; vNext; a; qNext; wNext; al];
end

% -------------------------------------------------------------------------
function z = fMeas(x), z = [x(1:3); x(10:13)]; end     % poseonly

% -------------------------------------------------------------------------
function omega = rotDiff2omega(R1,R2,dt)
% Average angular velocity in world frame over [t,t+dt]
R_rel = R2 * R1.';                                     % world rotation
axang = rotm2axang(R_rel);                             % angle  [0 ] rad
omega = (axang(4)/dt) * axang(1:3).';                  % rad/s, world frame
omega = omega.';
end