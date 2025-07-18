function distVal = transformDistance(T1, T2, wTrans, wRot)
% transformDistance:
%   Compute a weighted sum of translation and rotation differences 
%   between two transforms T1 and T2 in 4x4 homogeneous form.
%
%   T1, T2  : 4x4 homogeneous transforms
%   wTrans  : weight for translation error (defaults to 1 if not supplied)
%   wRot    : weight for rotation error (defaults to 1 if not supplied)

    if nargin < 3, wTrans = 1; end
    if nargin < 4, wRot   = 1; end

    % Extract rotation & translation
    R1 = T1(1:3,1:3);
    R2 = T2(1:3,1:3);
    t1 = T1(1:3,4);
    t2 = T2(1:3,4);

    % Translation error
    transError = norm(t1 - t2);

    % Rotation error (angle between R1 and R2)
    Rdelta = R1' * R2;    

    % The angle can be computed via trace
    val = (trace(Rdelta) - 1)/2;
    % (clamp inside [-1,1] to handle floating-point issues)
    val = min(max(val, -1), 1);
    % compute the angle
    theta = acos(val);

    % Combine
    distVal = wTrans * transError + wRot * theta;
end
