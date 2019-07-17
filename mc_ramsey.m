function [phi_est, y, N_spin_det, N_spin] = mc_ramsey(N,qe,phi,n_mc)
% Monte-carlo simulation of Ramsey interferometry with finite detector efficiency
%
%   [phi_est, y, N_spin_det, N_spin] = mc_ramsey(N_qubits,qe,phi,N_mc)
%   
%   phi_est = estimate of phase
%   y = normalised polarisation: (n_up - n_down) / (n_up + n_down)
%   N_spin_det = detected number of spin: [UP, DOWN] 
%   N_spin = number of spin (ideal detection)
% 
%   N_qubits = number of qubits
%   qe = detector efficiency
%   phi = ramsey phase
%   n_mc = number of sims
%   
%
% DKS
% 2019-07-17

% Ramsey
p_up = @(x) cos(x/2).^2;        % probability to measure spin |+>
y2phi = @(y) acos(y);


% simulate spin projection 
N_0 = binornd(N,p_up(phi),[n_mc,1]);
N_1 = N - N_0;
N_spin = [N_0, N_1];        % number in each spin state

% simulate detection
N_spin_det = binornd(N_spin,qe);

% evaluate measurement result
y = -diff(N_spin_det,1,2)./sum(N_spin_det,2);
phi_est = y2phi(y);

end