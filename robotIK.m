function qOpts = robotIK(eeTform, enforceJointLimits, sortByDistance, referenceConfig)
%robotIK Function for generating closed-form inverse kinematics solutions to the DH robot given by the parameters specified below
%   $Revision: $ $Date: $
%
%   Generated on 07-Nov-2022 18:30:26


dhParams = [0 -1.5707963267949 0 0;0.3 0 0 0;0 -1.5707963267949 1.47722600099026e-33 0;0 1.5707963267949 0.2 0;0 -1.5707963267949 -6.35179545195301e-35 0;0 0 0.06 0];
thetaOffsets = [0 0 0 0 0 0];
lastThreeAxesSign = [1 -1 1];
jointLimits = [-3.14159265358979 3.14159265358979;-3.14159265358979 3.14159265358979;-3.14159265358979 3.14159265358979;-3.14159265358979 3.14159265358979;-3.14159265358979 3.14159265358979;-3.14159265358979 3.14159265358979];
isJointRevolute = [1 1 1 1 1 1];

% Compute the shifted joint limits, which are the limits during solution, where theta offsets are not yet in play
shiftedJointLimits = jointLimits + repmat(thetaOffsets(:),1,2);

% Convert the end effector pose in the global frame to the end effector described by the DH parameters relative to the DH-described origin
baseToWorldTform = [1 0 0 0;0 1 0 0;0 0 1 -0.1;0 0 0 1];
actEEToDhEETform = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
eePose = baseToWorldTform*eeTform*actEEToDhEETform;

% Parse optional inputs
narginchk(1,4);
if nargin < 4
    referenceConfig = zeros(1,6);
end
if nargin < 3
    sortByDistance = false;
end
if nargin < 2
    enforceJointLimits = true;
end

% If joint limits are not enforced, set the shifted joint limits to have infinite range
if ~enforceJointLimits
    shiftedJointLimits = repmat([-inf inf], size(jointLimits, 1), 1);
end
% Map the desired end effector pose to the pose of the central intersecting joint.
eeToJt5 = [1 0 0 0;0 6.12323399573677e-17 -1 0;0 1 6.12323399573677e-17 -0.06;0 0 0 1];
jt5Pose = eePose*eeToJt5;

% Solve for the position of the first three joints from the pose of joint 5
q123Opts = solveFirstThreeDHJoints(jt5Pose(1:3,4), dhParams);

% Solve for the positions of the intersecting axes
% For each position solution, this configuration of the last three axes produces at least two possible orientation solutions
numRotationSolns = 2;
q456Opts = zeros(numRotationSolns*size(q123Opts,1), size(q123Opts,2));

% The next step seeks to find the orientation, which is entirely governed by the last three joints. This means that rotation from the fourth joint to the end effector can be mapped to three rotations in-place about the fifth joint. Since the rotations are in-place, they can be defined relative to the fourth joint axes, assuming a fixed pose rotation at the end to align with the end effector. The fixed rotation is found using the DH parameters, and corresponds to the rotation of the end effector relative to the fourth joint when the last three joints are all zero.
eeFixedAlpha = dhParams(4,2) + dhParams(5,2) + dhParams(6,2);
eeFixedRotation = [1 0 0; 0 cos(eeFixedAlpha) -sin(eeFixedAlpha); 0 sin(eeFixedAlpha) cos(eeFixedAlpha)];
for jtIdx = 1:size(q123Opts,1)
    % Get the position of the fourth joint at its zero position when the first three joints are positioned for IK
    jt4ZeroPose = getJoint4PoseFromDH(q123Opts(jtIdx,:));
    
    % Compute the rotation matrix needed to get to the end
    % The orientation of the end effector in the world frame can be written:
    %    eeRot = jt4ZeroRot*(Rotation about axes 4-6)*eeFixedRotation
    % Then the goal is to solve for the rotation about the axes and relate them to he known form from the DH parameters, if a valid solution exists:
    %    (Rotation about axes 4-6) = jt4ZeroRot'*eeRot*eeFixedRotation'
    jt4ToEERot = jt4ZeroPose(1:3,1:3)'*eePose(1:3,1:3)*eeFixedRotation';
    
    % This orientation produces at least two configurations for every solution, when joint limits allow
    orientationSolns = convertRotationToZYZAxesAngles(jt4ToEERot, lastThreeAxesSign, shiftedJointLimits(4:6,:));
    q456Opts(jtIdx,:) = orientationSolns(1,:);
    q456Opts(jtIdx + size(q123Opts,1),:) = orientationSolns(2,:);
    
    % Offset theta to reflect the source robot configuration
    q123Opts(jtIdx,:) = q123Opts(jtIdx,:) - thetaOffsets(1:3);
    q456Opts(jtIdx,:) = q456Opts(jtIdx,:) - thetaOffsets(4:6);
    q456Opts(jtIdx + size(q123Opts,1),:) = q456Opts(jtIdx + size(q123Opts,1),:) - thetaOffsets(4:6);
    
    % Remove solutions that violate joint limits
    if enforceJointLimits
        q123Opts(jtIdx,:) = applyJointLimits(q123Opts(jtIdx,:), jointLimits(1:3,:), isJointRevolute(1:3));
        q456Opts(jtIdx,:) = applyJointLimits(q456Opts(jtIdx,:), jointLimits(4:6,:), isJointRevolute(1:3));
        q456Opts(jtIdx + size(q123Opts,1),:) = applyJointLimits(q456Opts(jtIdx + size(q123Opts,1),:), jointLimits(4:6,:), isJointRevolute(1:3));
    end
    
end

% Filter out any remaining rows with NaNs in them by getting the index of the valid rows and only assembling those in the final output
allSolnOpts = [repmat(q123Opts, numRotationSolns, 1) q456Opts];
isValidRowIdx = all(~isnan(allSolnOpts),2);
qOptsAllSolns = allSolnOpts(isValidRowIdx,:);

% Create a copy of the solutions that wraps all revolute joints to pi, then round within solution tolerance.
qOptsWrappedAndRounded = round(robotics.internal.wrapToPi(qOptsAllSolns)*1.000000e+06)/1.000000e+06;

% Find the indices of all unique values after wrapping to pi
[~, isUniqueValidRowIdx] = unique(qOptsWrappedAndRounded, 'rows');

% Select only unique solutions from the original set of solutions
qOpts = qOptsAllSolns(sort(isUniqueValidRowIdx),:);
% Sort results using a distance metric
if sortByDistance
    qOpts = sortByEuclideanDistance(qOpts, referenceConfig(:)', isJointRevolute);
end

% Helper functions
    function sortedSolutions = sortByEuclideanDistance(solutions, referenceConfig, isJointRevolute)
        %sortByEuclideanDistance Sort a matrix of solution configurations relative to a reference configuration by Euclidean norm
        %   This function sorts a matrix of configurations using a pre-defined
        %   distance metric. The computed distance between any state and the
        %   reference state, referenceConfig, is a Euclidean norm of difference
        %   between a revolute joint's values which is then wrapped to [-pi, pi],
        %   and a displacement between a prismatic joint's values.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % Compute the distances between each configuration and the reference
        dist = robotics.manip.internal.RigidBodyTreeUtils.distance(...,
            referenceConfig, solutions, isJointRevolute);
        
        % Sort the outputs
        [~, sortedIdx] = sort(dist);
        sortedSolutions = solutions(sortedIdx,:);
        
    end

    function validConfig = applyJointLimits(inputConfig, jointLimits, isJointRevolute)
        %applyJointLimits Convert solutions with invalid joint limits to NaNs
        %   Given an N-element configuration, an Nx2 set of lower and upper joint
        %   limits, and an N-element vector indicating the joint type (revolute or
        %   prismatic), this function checks whether the configuration is within
        %   the joint limits. If not, the configuration is converted to NaNs.
        
        %   Copyright 2020-2021 The MathWorks, Inc.
        
        % Initialize output
        validConfig = inputConfig;
        
        for i = 1:numel(inputConfig)
            if jointLimits(i,1) > inputConfig(i) || jointLimits(i,2) < inputConfig(i)
                
                % Compute the offset from the lower joint limit and compare that to
                % the total range
                wrappedJointValueOffset = robotics.internal.wrapTo2Pi(inputConfig(i) - jointLimits(i,1));
                
                % If the wrapped value is 2*pi, make sure it is instead registered
                % as zero to ensure this doesn't fall outside the range
                if isEqualWithinTolerance(wrappedJointValueOffset, 2*pi)
                    wrappedJointValueOffset = 0;
                end
                
                jointRange = jointLimits(i,2) - jointLimits(i,1);
                
                if isJointRevolute(i) && ((wrappedJointValueOffset < jointRange) || isEqualWithinTolerance(wrappedJointValueOffset, jointRange))
                    
                    % Make sure the final value is definitively inside the joint
                    % limits if it was on the bound
                    wrappedJointValueOffset = min(wrappedJointValueOffset, jointRange);
                    
                    % Update the configuration
                    validConfig(i) = jointLimits(i,1) + wrappedJointValueOffset;
                else
                    % If any element is NaN, the whole array will be thrown out so
                    % there is no need to continue
                    validConfig = nan(size(validConfig));
                    return;
                end
            end
        end
        
    end

    function [actAngles, jointsInGimbalLock] = convertRotationToZYZAxesAngles(tgtRotation, axesSign, jointLim)
        %convertRotationToZYZAxesAngles Convert desired orientation to rotation about Z-Y-Z
        %   This function is used to three angles of rotation corresponding to
        %   consecutive joint angles whose joint axes intersect at a single, common
        %   point. This function addresses the case where the first joint rotation
        %   is about Z, and the subsequent angles, in order and defined relative to
        %   the first joint frame, are about Y and then Z. The function accepts the
        %   rotation of the last joint relative to the origin of the first one, as
        %   well as the sign of each axes. The second output indicates joints that
        %   are in gimbal lock, where 1 indicates that they are, and zero indicates
        %   that they are not. When joints are in gimbal lock, the affected joint
        %   axes align and an infinite combination of joint angles are possible (as
        %   long as the total rotation still matches the target). The default
        %   assumption is that rotation is divided over the joint along the
        %   affected axis.
        %   Copyright 2020 The MathWorks, Inc.
        eulAngles = rotm2eul(tgtRotation, 'ZYZ');
        
        % The jointsInGimalLock variable indicates redundant joints, i.e. joints
        % that complement each other and can have an infinite pair of values in the
        % directJointAngleMaps output. Initialize this value to zeros (no joints in
        % gimbal lock). This is a flag that is consistent across rotation functions
        % and may be used by the caller function.
        jointsInGimbalLock = [0 0 0];
        % When the middle angle is zero, the first and last joints are co-axial,
        % meaning there are an infinite set of solutions. Use a helper function to
        % distribute the values consistently given knowledge of the joint limits.
        if isEqualWithinTolerance(eulAngles(2), 0)
            
            newTgtRotation = tgtRotation;
            newTgtRotation(1:2,3) = 0;
            newTgtRotation(3,1:2) = 0;
            newTgtRotation(3,3) = 1;
            eulAngles = rotm2eul(newTgtRotation, 'ZYZ');
            variableJtIdx = [1 3];
            jointsInGimbalLock(variableJtIdx) = [1 1];
            totalRotation = sum(eulAngles(variableJtIdx));
            eulAngles(variableJtIdx) = distributeRotationOverJoints(totalRotation, axesSign(variableJtIdx), jointLim(variableJtIdx,:));
            
            % In this case the alternate Euler angles aren't required, as they will
            % also result in a set of co-axial joints. However, to ensure codegen
            % compatibility, the size must stay the same Therefore, return a set of
            % NaN angles (so the second solution may be thrown out). Note that the
            % axes sign is handled inside the distributeRotationOverJoints helper
            % function.
            actAngles = [eulAngles; nan(1,3)];
        else
            % For finite solutions when the middle angle is non-zero, there are two possible solutions to this problem
            % that can be derived from the first solution set
            eulAltUnwrapped = eulAngles;
            eulAltUnwrapped(:,2) = -eulAltUnwrapped(:,2);
            eulAltUnwrapped = eulAltUnwrapped + pi;
            eulAltUnwrapped(:,2) = eulAltUnwrapped(:,2) - pi;
            eulAnglesAlt = robotics.internal.wrapToPi(eulAltUnwrapped);
            
            % Output the angles given the axes signs
            actAngles = [eulAngles; eulAnglesAlt]*diag(axesSign);
        end
    end

    function jointAngles = distributeRotationOverJoints(totalRotation, axesSigns, jointLim)
        %distributeRotationOverJoints Distribute a rotation over several in-line revolute joints
        %   When revolute joints are co-axial, the total rotation can be distributed
        %   over the joints in a number of ways. This function assigns a default
        %   behavior that respects the joint limits. For the case where no joint
        %   limits are required, they should be provided as infinite, i.e [-inf
        %   inf] for each joint. The default behaviors are as follows:
        %
        %      - If any joint limits have a range of at minimum 2*pi, all total
        %        rotation amounts are achievable and the rotation is distributed
        %        evenly over the joints with infinite range, assuming the other
        %        joints are centered within their range.
        %
        %      - If all joints have finite ranges with total range less than 2*pi,
        %        some total rotation amounts may not be feasible and the rotation
        %        is distributed as much as possible on the distal joints (unused
        %        more proximal joints are centered). If the solution is infeasible,
        %        a NaN-vector is returned.
        %
        %   The goal of these distributions is to favor solutions that are
        %   efficient. This function accepts the total rotation (in radians) to be
        %   divided over N joints, the signs of those N joints (whether rotation is
        %   positive or negative about the axes), and the joint limits, given as an
        %   Nx2 matrix.
        %
        %   If joint limits are ignored, they can be provided as infinite; the
        %   behavior is equivalent. This function returns an N-element row vector.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % Get the total number of joints from the joint limit input
        N = size(jointLim, 1);
        
        % Initialize the output
        jointAngles = zeros(1,N);
        
        % Remap the joint limits to fit the assumption that all axes are positive.
        % Since the joint limits can contain infinite elements, it is necessary to
        % use element-wise multiplication, as matrix multiplication can result in
        % NaNs when it causes sums of infinities.
        jointLim = repmat(axesSigns(:),1,2).*jointLim;
        
        % Re-order the joint limits to ensure the lower limit always comes first
        % (in case the of a sign flip in the previous line)
        jointLim = sort(jointLim,2);
        
        % Determine the total ranges of each joint. Since all joints are revolute,
        % a range of 2*pi or greater is equivalent to an infinite range as the IK
        % problem does not distinguish between periodic equivalents. Note that a
        % downstream helper in the IK solution, applyJointLimits, includes a check
        % that wraps periodic equivalents back to their values given the joint
        % limits.
        jointRange = jointLim(:,2) - jointLim(:,1);
        isRevJointFullRange = (jointRange > 2*pi);
        for limIdx = 1:size(jointRange,1)
            % Use a tolerance check on the equality. Since isEqualWithinTolerance
            % returns a scalar, it is necessary to do this inside a for-loop
            isRevJointFullRange(limIdx) = isRevJointFullRange(limIdx) || isEqualWithinTolerance(jointRange(limIdx), 2*pi);
        end
        
        % There are two primary cases: when some joints have full range, any
        % solution is feasible and the variable values are distributed over these
        % joints. When all of the joints have range of less than 2*pi, the problem
        % is more complex, as some solutions may not be feasible.
        if any(isRevJointFullRange)
            % If any of the joint have infinite ranges, use that to distribute the
            % total rotation. First, place the joints with finite ranges in the
            % middle of their respective ranges, then distribute the remaining
            % joint rotation evenly over the joints with at least 2*pi range.
            jointIdxVector = 1:N;
            jointsWithIncompleteRange = jointIdxVector(~isRevJointFullRange);
            for i = 1:numel(jointsWithIncompleteRange)
                jointIdx = jointsWithIncompleteRange(i);
                jointAngles(jointIdx) = sum(jointLim(jointIdx,:))/2;
            end
            
            % Compute the remaining rotation and wrap it to the interval [-pi, pi],
            % then distribute over the joints with complete range
            wrappedRemainder = robotics.internal.wrapToPi(totalRotation - sum(jointAngles));
            jointsWithCompleteRange = jointIdxVector(isRevJointFullRange);
            for j = 1:numel(jointsWithCompleteRange)
                jointIdx = jointsWithCompleteRange(j);
                jointAngles(jointIdx) = wrappedRemainder/numel(jointsWithCompleteRange);
            end
            
        else
            % Use an algorithm that favors loading distal joints, which are
            % typically easier to change: first set all the joints to their
            % mid-range values. Then iterate over the joints from first to last,
            % moving each joint up or down based on the difference in the current
            % total from the desired total, until the desired total is reached.
            % This is essentially a cascaded bang-bang controller.
            
            % Initialize the joint angles to their mid-range values
            jointAngles(:) = (sum(jointLim,2)/2)';
            
            % Iterate over the joints, using a feedback law to move them closer
            % to the desired total
            jointIdxVector = N:-1:1;
            wrappedTotalRotation = robotics.internal.wrapToPi(totalRotation);
            for jointIdx = 1:numel(jointIdxVector)
                diffRotation = robotics.internal.wrapToPi(wrappedTotalRotation - sum(jointAngles));
                jointAngles(jointIdx) = jointAngles(jointIdx) + sign(diffRotation)*min(abs(diffRotation), jointRange(jointIdx)/2);
            end
            
            % Check if the sum of the joint angles reaches the desired total. If
            % not, the solution is infeasible and a vector of NaNs is returned.
            if ~isEqualWithinTolerance(robotics.internal.wrapToPi(sum(jointAngles)), wrappedTotalRotation)
                jointAngles = nan(size(jointAngles));
                return;
            end
        end
        
        % Factor back in the axes signs. Since all valid joint angles are finite,
        % matrix multiplication is the most efficient approach.
        jointAngles = jointAngles*diag(axesSigns);
        
    end

end

function outputThetas = solveFirstThreeDHJoints(jt5Pos, dhParams)
%solveFirstThreeDHJoints Solve for the first three joint angles of a DH-parameterized robot
%   This function computes the first three joint angles of a robot
%   parameterized using Denavit-Hartenberg parameters. The function accepts
%   a matrix of the fixed DH parameters, as well as the position of the
%   fifth joint. The matrix of DH parameters is of size 6x4 for the 6
%   non-fixed joints, where each row has the order [a alpha d 0], where a
%   is a translation along x, alpha is a rotation about x, and d is a
%   translation along z. The last value, which typically refers to theta
%   (the rotation about z) for that joint, is not yet known; this function
%   will solve for theta for the first three joints. When a robot has the
%   last three axes intersecting, the position and orientation of the end
%   effector can be split up: the position is entirely determined by the
%   first three joints, while the orientation is governed by the last three
%   joints (provided the first three are known). Furthermore, the position
%   of any end effector can be related to the position of the fifth joint
%   frame, which corresponds to the joint frame at the midpoint of the
%   three intersecting joints. This frame will have the same position in
%   any valid solution to the particular IK problem (though its orientation
%   may differ), and its translation relative to the base frame is entirely
%   defined by the first three joints. This function solves for those first
%   three joints given the position of the joint 5 frame relative to the
%   base. This solution method and notation follows Chp.3 of Pieper's 1968
%   thesis, but adds two corrections, as well as minor notation changes and
%   the addition of constraints to ensure only feasible solutions are
%   output:
%
%   Pieper, D. The Kinematics Of Manipulators Under Computer Control.
%   Stanford University (1968).

% Extract DH parameters from matrix
[a1, a2, a3] = deal(dhParams(1,1), dhParams(2,1), dhParams(3,1));
[alpha1, alpha2, alpha3] = deal(dhParams(1,2), dhParams(2,2), dhParams(3,2));

% Note that Pieper uses "s" instead of "d" in his solutions
[d1, d2, d3, d4] = deal(dhParams(1,3), dhParams(2,3), dhParams(3,3), dhParams(4,3));

% Three variables derived from jt5Pos
z3 = jt5Pos(3);
R3 = jt5Pos(1)^2 + jt5Pos(2)^2 + (jt5Pos(3) - d1)^2;
z = z3 - d1;

% The first step is to solve for theta3. This is achieved by eliminating
% theta1 and theta2 through a number of substitutions and sum-of-squares
% operations. The resultant equation for theta3, is a function of
% sin(theta3) and cos(theta3), but this can be further mapped to a
% polynomial in h, where h = tan(theta3/2). This substitutions is made
% possible by the Weierstrass transformation, which maps sin(theta3) to
% 2*h/(1 + h^2) and cos(theta3) to (1-h^2)/(1+h^2). The goal is then to
% solve the polynomial for h and map the solutions back to theta3.

% Since a1 = 0, the solution arises from R3 = F3, which produces a quadratic in h
[hSolns, ~, hasPiSoln] = solveForHCaseZeroA1(R3, a2, a3, alpha2, alpha3, d2, d3, d4);
% Initialize the matrix of possible solutions
possThetas = zeros(16,3);
% After all solutions are processed, rows with NaNs will be removed
% Initialize theta3 to NaN and replace based on actual solutions
possThetas(:,3) = NaN;
for hIdx = 1:numel(hSolns)
    % Ensure only real solutions to h are converted
    h3 = replaceImagWithNaN(hSolns(hIdx));
    if isnan(h3)
        % When h3 is imaginary, theta3 = NaN
        possThetas(hIdx,3) = NaN;
        possThetas(4+hIdx,3) = NaN;
    else
        % When h3 is real, there are two possible equivalent values of theta3
        possThetas(hIdx,3) = 2*atan2(h3,1);
        possThetas(4+hIdx,3) = 2*atan2(-h3,-1);
    end
end
if hasPiSoln
    possThetas(numel(hSolns)+1,3) = pi;
end
for theta3Idx = 1:8
    % If theta3 is NaN or imaginary, replace whole row with NaNs and skip to next row
    theta3 = replaceImagWithNaN(possThetas(theta3Idx,3));
    if isnan(possThetas(theta3Idx,3))
        possThetas(theta3Idx,:) = [NaN NaN NaN];
        continue
    end
    
    % Compute key subexpressions f1 to f3 and F1 to F4, which are functions of theta3
    f = computef13SupportingEquations(a3, alpha3, d3, d4, theta3);
    F = computeF14SupportingEquations(a1, a2, alpha1, alpha2, f(1), f(2), f(3), d2);
    % Compute theta2. The exact approach depends on the DH
    % parameters, but the gist is the same: since the equations
    % output multiple solutions, but some are actually just results
    % of the sum of squares, i.e., they solve the local problem,
    % but do not actually solve the overlying problem. Rather than
    % compute all solutions and filter at the end, we filter here
    % by always solving using two different equations. Then we
    % choose only the solution that satisfies both equations.
    
    % Since a1 is zero and sin(alpha1) is nonzero, solve for theta2 using equation 3.26 and the third element of 3.20
    theta2Opts = solveTrigEquations(-F(2)*sin(alpha1), F(1)*sin(alpha1), z - F(4));
    
    rhsConst = jt5Pos(3) - d1 - cos(alpha1)*(sin(alpha2)*f(2) + cos(alpha2)*f(3) + d2);
    lhsCosCoeff = -sin(alpha1)*(-cos(alpha2)*f(2) + sin(alpha2)*f(3));
    lhsSinCoeff = sin(alpha1)*(a2 + f(1));
    theta2Constraint = solveTrigEquations(lhsCosCoeff, lhsSinCoeff, rhsConst);
    
    % Choose the solution(s) that solve both equations
    theta2 = chooseCorrectSolution(theta2Opts, theta2Constraint, 1.000000e-06);
    
    % Theta2 is a 2-element vector with up to two valid solutions (invalid
    % solutions are represented by NaNs). Iterate over the possible values
    % and add the second solution set in the latter half of the matrix (so
    % they aren't overwritten by subsequent loops).
    for theta2Idx = 1:2
        % Update the local index so it's reflective of the indexed value of theta2
        solIdx = theta3Idx + 8*(theta2Idx-1);
        
        % Update the value of theta3 in case it was previously set to NaN,
        % and replace any invalid values of theta2 with NaN
        possThetas(solIdx,3) = theta3;
        possThetas(solIdx,2) = replaceImagWithNaN(theta2(theta2Idx));
        
        % If any of the joint variables in NaN, replace it and all the
        % remaining joints to solve with NaNs and move on to the next loop
        if isnan(possThetas(solIdx, 2))
            possThetas(solIdx, 1:2) = [NaN NaN];
            continue;
        end
        
        % Compute theta1 from the first two elements of eq 3.20
        g = computeg12SupportingEquations(a1, a2, alpha1, alpha2, f(1), f(2), f(3), d2, theta2(theta2Idx));
        theta1Opts = solveTrigEquations(g(1), g(2), jt5Pos(1));
        theta1Constraint = solveTrigEquations(-g(2), g(1), jt5Pos(2));
        theta1Opts = chooseCorrectSolution(theta1Opts, theta1Constraint, 1.000000e-06);
        
        % Since theta1 is the last value that is solved for, only one
        % of the solutions will be valid, and chooseCorrectSolution
        % sorts the results so that if there is only one solution, it
        % is always the first element (and the other element is nan)
        theta1 = theta1Opts(1);
        
        % Update the array of possible theta values
        possThetas(solIdx,1) = replaceImagWithNaN(theta1);
        
    end
    
end

% Now we are left with an 8x3 matrix where some values are NaN. The
% function will only output the rows where all elements are non-NaN.
outputThetas = possThetas(all(~isnan(possThetas),2),:);
% Helper functions
    function f = computef13SupportingEquations(a3,alpha3,s3,s4,theta3)
        %computef13SupportingEquations Compute f1 to f3 supporting equations
        %   This function computes f1 to f3, which are functions of theta3. For a
        %   given robot with three consecutive revolute axes, the position of a
        %   joint a distance s4 along the joint 3 axis can be described as:
        %      P = T1*T2*T3*[0; 0; s4; 1],
        %   where Ti represent transformation matrices associated with links. Then
        %   this equation may be rewritten as P = T1*T2*f. This function computes
        %   the values of f that satisfy the rewritten equation.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % Initialize output
        f = zeros(3,1);
        
        % Compute component terms
        t2 = sin(alpha3);
        t3 = cos(theta3);
        t4 = sin(theta3);
        
        % Assemble outputs. Note that there is technically a fourth output, f(4) =
        % 1, but its value is unused, so it is not computed or returned.
        f(1) = a3.*t3+s4.*t2.*t4;
        f(2) = a3.*t4-s4.*t2.*t3;
        f(3) = s3 + s4.*cos(alpha3);
        
    end


    function F = computeF14SupportingEquations(a1,a2,alpha1,alpha2,f1,f2,f3,s2)
        %computeF14SupportingEquations Compute intermediate variables F1 to F4
        %   This function computes F1 to F4, which are intermediate variables in
        %   Pieper's derivation that are functions of the theta3 joint position
        %   (and constant parameters). The function accepts several DH parameters,
        %   as well as the intermediate variables f1 to f3, which are functions of
        %   theta3, and outputs the four F1 to F4 intermediate variables.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % Initialize output
        F = zeros(1,4);
        
        F(1) = a2+f1;
        
        t2 = cos(alpha2);
        t3 = sin(alpha2);
        F(2) = -f2.*t2+f3.*t3;
        
        t4 = f3.*t2;
        t5 = f2.*t3;
        F(3) = s2.*(t4+t5).*2.0+a2.*f1.*2.0+a1.^2+a2.^2+f1.^2+f2.^2+f3.^2+s2.^2;
        
        F(4) = cos(alpha1).*(s2+t4+t5);
        
    end

    function g = computeg12SupportingEquations(a1,a2,alpha1,alpha2,f1,f2,f3,s2,theta2)
        %computeg12SupportingEquations Compute g1 and g2 supporting equations
        %   This function computes g1 and g2, which are functions of theta2 and
        %   theta3.
        
        %   Copyright 2020 The MathWorks, Inc.
        % Initialize output
        g = zeros(1,2);
        
        % Compute component terms
        t2 = cos(alpha1);
        t3 = cos(alpha2);
        t4 = sin(alpha2);
        t5 = cos(theta2);
        t6 = sin(theta2);
        t7 = a2+f1;
        t8 = f2.*t3;
        t9 = f3.*t4;
        t10 = -t9;
        t11 = t8+t10;
        
        % Assemble outputs
        g(1) = a1+t5.*t7-t6.*t11;
        g(2) = sin(alpha1).*(s2+f2.*t4+f3.*t3)-t2.*t6.*t7-t2.*t5.*t11;
        
    end


    function theta = solveTrigEquations(a,b,c)
        %solveTrigEquations Solve equations of the form a*cos(theta) + b*sin(theta) = c for theta
        %   This function solves the common trigonometric equality by equating the
        %   solution to cos(phi)sin(theta) + sin(phi)cos(theta) = sin(phi + theta).
        %   The function returns two possible solutions for theta.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        theta = nan(1,2);
        % Handle the trivial case
        if isEqualWithinTolerance(a,0) && isEqualWithinTolerance(b,0) && isEqualWithinTolerance(c,0)
            theta(1) = 0;
            return;
        elseif isEqualWithinTolerance(a,0) && isEqualWithinTolerance(b,0) && ~isEqualWithinTolerance(c,0)
            return;
        end
        % As long as a or b are nonzero, a set of general solutions may be found
        d = sqrt(a^2 + b^2);
        cPrime = c/d;
        if cPrime < 1 || isEqualWithinTolerance(cPrime,1)
            % Throw out the imaginary solutions, which occur when cPrime > 1
            phi1 = atan2(a,b);
            phi2 = atan2(-a,-b);
            theta(1) = real(asin(complex(cPrime))) - phi1;
            theta(2) = -real(asin(complex(cPrime))) - phi2;
        end
        
    end

    function correctSolution = chooseCorrectSolution(solutionPair1, solutionPair2, solTolerance)
        %chooseCorrectSolution Choose the solution that appears in both solutionPair1 and solutionPair2
        %   This helper function is used to choose a correct solution when two
        %   solutions are provided, e.g. as the result of a sum of squares. The
        %   function accepts two 2-element vectors, solutionPair1 and
        %   solutionPair2, which represent the solution options from the source and
        %   constraint equations, respectively. The correct solution will be the
        %   solution both of the source equation, as well as a constraint equation
        %   for the same problem. This helper simply chooses the value that occurs
        %   in both the original and constraint solutions, within a tolerance.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % Filter any imaginary values out of the solution pairs by replacing them
        % with NaNs
        realSolutionPair1 = zeros(1,2);
        realSolutionPair2 = zeros(1,2);
        for i = 1:2
            % Have to wrap to pi so that the two values are comparable
            realSolutionPair1(i) = robotics.internal.wrapToPi(replaceImagWithNaN(solutionPair1(i)));
            realSolutionPair2(i) = robotics.internal.wrapToPi(replaceImagWithNaN(solutionPair2(i)));
        end
        
        % To check equivalence, it's insufficient to just check whether the values
        % are equal, because they are periodic. For example, -pi and pi are both
        % valid outcomes of wrapToPi that fail a basic equality test but are
        % equivalent in this context. Therefore, it's necessary to check that the
        % difference of the two values, when wrapped to pi, is inside the expected tolerance.
        correctSolution = nan(1,2);
        for i = 1:2
            for j = 1:2
                isCorrect = abs(robotics.internal.wrapToPi(realSolutionPair1(i) - realSolutionPair2(j))) < solTolerance;
                if isCorrect
                    correctSolution(i) = realSolutionPair1(i);
                end
            end
        end
        
        % Sort the output so that if there is one correct solution it is always in
        % the first element slot
        correctSolution = sort(correctSolution);
        
    end

    function q = replaceImagWithNaN(qToCheck)
        %replaceImagWithNaN Replace imaginary and empty elements with NaNs
        %   This function replaces imaginary values with NaNs. This is useful when
        %   the element is part of a matrix, and rendering one element of the
        %   matrix imaginary will make the entire matrix imaginary. Furthermore, it
        %   may be used to filter invalid solutions.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        if isempty(qToCheck)
            q = NaN;
        elseif ~isEqualWithinTolerance(imag(qToCheck), 0)
            q = NaN;
        else
            q = real(qToCheck);
        end
        
    end

    function [hSolns, hasFiniteNumSol, hasPiSoln] = solveForHCaseZeroA1(R3, a2, a3, alpha2, alpha3, d2, d3, d4)
        %solveForHCaseZeroA1 Solve for h when a1 is zero
        %   To solve for theta3, it is necessary to reparameterize a trigonometric
        %   equation in terms of a new parameter h = tan(theta3/2) using the
        %   Weierstrass equation. This function solves equation 3.25 in the Pieper
        %   inverse kinematics solution for h in the case where DH parameter a1 is
        %   zero. In that case, equation 3.25 becomes:
        %      R3 = F3
        %   Here F1 to F4 are functions of theta3 (and the constant DH parameters),
        %   and R3 is a function of P, a known position input from the IK problem,
        %   and the DH parameter d1:
        %      R3 = P(1)^2 + P(2)^2 + (P(3) – d1)^2
        %   Equation 3.25 may then may be reparameterized in h, producing a
        %   quadratic polynomial in h. This function solves that polynomial for the
        %   values of h given R3 and the DH parameters of the associated serial
        %   manipulator.
        
        %   Copyright 2020 The MathWorks, Inc.
        
        % These solutions solve the equation above, R3 = F3, which becomes a
        % quadratic in h = tan(theta3/2). As the polynomial is quadratic, the
        % solution has a format that matches that of the solutions to the quadratic
        % equation A*h^2 + B*h + C = 0.
        % Compute the terms of the quadratic equation
        A = a2^2 - 2*a2*a3 - R3 + a3^2 + d2^2 + d3^2 + d4^2 + 2*d2*d3*cos(alpha2) + 2*d3*d4*cos(alpha3) + 2*d2*d4*cos(alpha2)*cos(alpha3) + 2*d2*d4*sin(alpha2)*sin(alpha3);
        B = 4*a3*d2*sin(alpha2) + 4*a2*d4*sin(alpha3);
        C = 2*a2*a3 - R3 + a2^2 + a3^2 + d2^2 + d3^2 + d4^2 + 2*d2*d3*cos(alpha2) + 2*d3*d4*cos(alpha3) + 2*d2*d4*cos(alpha2)*cos(alpha3) - 2*d2*d4*sin(alpha2)*sin(alpha3);
        % There are three possible solution cases
        if isEqualWithinTolerance([A B C], [0 0 0])
            % The trivial case happens whenever any value of theta3 solves the
            % quadratic derived from equation 3.25. In that case, any theta3 may
            % solve the problem, though it may be further constrained by equation
            % 3.26. Physically, this happens when theta3 has no impact on the
            % position solution, or when its solution is intertwined with theta2.
            hSolns = 0;
            hasFiniteNumSol = false;
        elseif isEqualWithinTolerance(A, 0)
            % When the first term is zero, the equation is linear
            hSolns = -C/B;
            hasFiniteNumSol = true;
        else
            % The equation is quadratic
            if B^2 - 4*A*C < 0
                % This solution will be complex
                h1 = (-B - sqrt(complex(B^2 - 4*A*C)))/(2*A);
                h2 = (-B + sqrt(complex(B^2 - 4*A*C)))/(2*A);
            else
                % This solution will be real
                h1 = real((-B - sqrt(B^2 - 4*A*C))/(2*A));
                h2 = real((-B + sqrt(B^2 - 4*A*C))/(2*A));
            end
            hSolns = [h1 h2];
            hasFiniteNumSol = true;
        end
        % Check if there is a solution at theta3 = pi, for which h is undefined, by
        % checking if R3 = F3 (eq 3.25) is satisfied for that solution.
        if hasFiniteNumSol
            localA1 = 0;
            localTheta = pi;
            fTerms = computef13SupportingEquations(a3, alpha3, d3, d4, localTheta);
            
            % alpha1 is required to use the standard function that computes F, but
            % it only affects the fourth term, F4, which is unused here. Therefore,
            % a dummy value is used.
            dummyAlpha1 = 0;
            FTerms = computeF14SupportingEquations(localA1, a2, dummyAlpha1, alpha2, fTerms(1), fTerms(2), fTerms(3), d2);
            
            hasPiSoln = isEqualWithinTolerance(FTerms(3), R3);
        else
            hasPiSoln = true;
        end
        
    end

end

function jt4Pose = getJoint4PoseFromDH(q123)
%getJoint4PoseFromDH Get the pose of the fourth joint when the first three joints set to q123 and joint4 angle is zero

dhParams = [0 -1.5707963267949 0 0;0.3 0 0 0;0 -1.5707963267949 1.47722600099026e-33 0;0 1.5707963267949 0.2 0;0 -1.5707963267949 -6.35179545195301e-35 0;0 0 0.06 0];

% Initialize output
jt4Pose = eye(4);
for i = 1:3
    a = dhParams(i,1);
    alpha = dhParams(i,2);
    d = dhParams(i,3);
    theta = q123(i);
    
    Ttheta = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
    TFixed = [1 0 0 a; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) d; 0 0 0 1];
    
    jt4Pose = jt4Pose*Ttheta*TFixed;
end

end

function isEquiv = isEqualWithinTolerance(mat1, mat2)
%isEqualWithinTolerance Check if two matrices are equal within a set tolerance
%   This is a convenience function designed for inputs with up to two
%   dimensions. If the input has 3+ dimensions, a non-scalar output will be
%   returned.

tol = 1.000000e-06;
diffMat = abs(mat1 - mat2);
isEquiv = all(all(diffMat < tol));

end
