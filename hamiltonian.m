function result = estimate_thetas(time, q1, q2, q1dot, q2dot, v, varargin)
% ESTIMATE_THETAS  Estimate parameter vector theta (6x1) using energy method.
%
% result = estimate_thetas(time,q1,q2,q1dot,q2dot,v)
%   Inputs (all vectors, same length N):
%     time   - time vector (seconds)
%     q1     - measured joint 1 angle (rad)
%     q2     - measured joint 2 angle (rad)
%     q1dot  - measured joint 1 velocity (rad/s). If empty, function will compute from q1
%     q2dot  - measured joint 2 velocity (rad/s). If empty, function will compute from q2
%     v      - motor input (voltage) applied (same length)
%
%   Optional Name-Value:
%     'g'         - gravity (default 9.81)
%     'useRows'   - subset of rows to use (default [] meaning use 2:end)
%     'regularize' - Tikhonov regularization lambda (default 0)
%
%   Output (struct):
%     result.theta       - 6x1 estimated parameter vector
%     result.A           - regression matrix used (M x 6)
%     result.d           - data vector (M x 1)
%     result.residuals   - d - A*theta
%     result.rmse        - root mean squared error
%     result.condA       - condition number of A
%     result.info        - textual notes about sign-convention and warnings
%
% Reference: Hamiltonian (energy) method as in your System Identification notes.

% -----------------------
% Input checks & defaults
p = inputParser;
addRequired(p,'time',@isvector);
addRequired(p,'q1',@isvector);
addRequired(p,'q2',@isvector);
addRequired(p,'q1dot',@(x) isvector(x) || isempty(x));
addRequired(p,'q2dot',@(x) isvector(x) || isempty(x));
addRequired(p,'v',@isvector);
addParameter(p,'g',9.81,@isscalar);
addParameter(p,'useRows',[],@(x) isempty(x) || (isvector(x) && all(x>=1)));
addParameter(p,'regularize',0,@(x)isnumeric(x) && x>=0);
parse(p,time,q1,q2,q1dot,q2dot,v,varargin{:});
g = p.Results.g;
useRows = p.Results.useRows;
%lambda = p.Results.regularize;

time = time(:);
q1 = q1(:); q2 = q2(:); v = v(:);
N = numel(time);
if numel(q1)~=N || numel(q2)~=N || numel(v)~=N
    error('time, q1, q2, and v must be same length N');
end

% compute velocities if not provided
if isempty(q1dot)
    q1dot = gradient(q1, time);
else
    q1dot = q1dot(:);
end
if isempty(q2dot)
    q2dot = gradient(q2, time);
else
    q2dot = q2dot(:);
end
if numel(q1dot)~=N || numel(q2dot)~=N
    error('All signals must be length N');
end

% -----------------------
% Compute basis functions h_i(t) (values at each time)
% According to your notes:
% h1 = 1/2 * q1dot^2
% h2 = 1/2 * sin(q2)^2 * q1dot^2 + 1/2 * q2dot^2
% h3 = cos(q2) * q1dot .* q2dot
% h4 = g * cos(q2)
h1 = 0.5 .* (q1dot.^2);
h2 = 0.5 .* ( (sin(q2)).^2 .* q1dot.^2 ) + 0.5 .* (q2dot.^2);
h3 = cos(q2) .* (q1dot .* q2dot);
h4 = g .* cos(q2);

% Integrals (from t0 to tk) using trapezoidal cumulative integration
int_vq1dot = cumtrapz(time, v .* q1dot);        % d entries
int_q1dot2   = cumtrapz(time, q1dot.^2);        % for F(:,1)
int_q2dot2   = cumtrapz(time, q2dot.^2);        % for F(:,2)

% Build H(tk,t0) = h(q(tk),qdot(tk)) - h(q(t0),qdot(t0))
h1_0 = h1(1); h2_0 = h2(1); h3_0 = h3(1); h4_0 = h4(1);
H1 = h1 - h1_0;
H2 = h2 - h2_0;
H3 = h3 - h3_0;
H4 = h4 - h4_0;

% By construction at k=1 (t=t0) all entries are zero; it's fine but provides trivial row.
% It's typical to drop the first row (all zeros) to avoid singularities.
rows = 2:N;
if ~isempty(useRows)
    rows = intersect(rows, useRows);
end

% Create d and A (use selected rows)
d = int_vq1dot(rows);
A = [ H1(rows) H2(rows) H3(rows) H4(rows) int_q1dot2(rows) int_q2dot2(rows) ];

% Solve least squares with optional regularization
%if lambda > 0
    % Tikhonov regularization (lambda * I)
    %theta = (A' * A + lambda * eye(6)) \ (A' * d);
%else
    theta = pinv(A) * d;   % pseudo-inverse (robust)
%end

% Diagnostics
res = d - A * theta;
rmse = sqrt(mean(res.^2));
condA = cond(A);

% Prepare output
result.theta = theta;
result.A = A;
result.d = d;
result.residuals = res;
result.rmse = rmse;
result.condA = condA;
result.info = sprintf(['Notes:\n' ...
    '- Theta is returned as [theta1; theta2; theta3; theta4; theta5; theta6].\n' ...
    '- theta5 and theta6 correspond to the coefficients multiplying integral(q1dot^2) and integral(q2dot^2).\n' ...
    '  In your derivation these terms came from -beta1*q1dot^2 - beta2*q2dot^2, so interpret signs accordingly (if you expect positive viscous betas, theta5/6 may be -beta1/-beta2 depending on sign convention).\n' ...
    '- First sample (t0) was dropped to avoid trivial zero row.\n' ...
    '- If A is ill-conditioned (condA large), consider exciting the system more or using regularization.']);
end


% time = out.tout;
% q1 = squeeze(out.q1_sim.signals.values);
% q2 = squeeze(out.q2_sim.signals.values);
% voltage = out.voltage.signals.values;
time = q1_hil.time;
q1 = squeeze(q1_hil.signals.values);
q2 = squeeze(q2_hil.signals.values);
q1dot = squeeze(q1dot_hil1.signals.values);
q2dot = squeeze(q2dot_hil1.signals.values);
voltage = voltage_hil.signals.values;
%voltage = voltage*(-1);

% Example: if velocities were not saved, pass [] and function computes gradients

% time = ScopeData(:,1);
% voltage = ScopeData(:, 2);
% q1 = ScopeData(:,3);
% q1_dot = ScopeData(:,4);
% q2 = ScopeData(:,5);
% q2_dot = ScopeData(:,6);

%res = estimate_thetas(time, q1, q2, q1_dot, q2_dot, voltage);
res = estimate_thetas(time, q1, q2, q1dot, q2dot, voltage);

disp('Estimated theta:');
disp(res.theta);

