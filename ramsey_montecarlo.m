%% Ramsey interferometry with non-ideal detection efficiency
% when measurement is normalised polarisation...
%
% DKS
% 2019-07-17

%% Config
% Ramsey
p_up = @(phi) cos(phi/2).^2;        % probability to measure spin-UP
y2phi = @(y) acos(y);

phi_ramsey = linspace(0,2*pi,1e2);       % phase
N_qubits = 30;         % number of qubits
p_detect = 0.1;         % detection efficiency


% Monte-carlo
N_mc = 1e3;          % repetition


%% Main 
y_mc = NaN(length(phi_ramsey),2);
phi_est = NaN(length(phi_ramsey),2);

for ii=1:length(phi_ramsey)
    [tphi_est, ty, N_spin_det, N_spin] = mc_ramsey(N_qubits,p_detect,phi_ramsey(ii),N_mc);
    
    % get statistics
    ty_mean = mean(ty,'omitnan');
    ty_std = std(ty,0,'omitnan');
    
    tphi_est_mean = mean(tphi_est,'omitnan');
    tphi_est_std = std(tphi_est,0,'omitnan');
    
    y_mc(ii,:) = [ty_mean, ty_std];
    phi_est(ii,:) = [tphi_est_mean, tphi_est_std];
end

%% plot
H_mc = figure;

% N_qubits = 64;         % number of qubits
% p_detect = 0.1;         % detection efficiency

%%%
subplot(3,1,1);


p_MC_y = ploterr(phi_ramsey,y_mc(:,1),[],y_mc(:,2),'o','hhxy',0);
set(p_MC_y(1),'MarkerEdgeColor','k','MarkerFaceColor','w',...
    'LineWidth',1);
set(p_MC_y(2),'Color','k','LineWidth',1);

xlabel('true phase, $\phi$ (rad)');
ylabel('polarisation');

titlestr = sprintf('N qubits = %d; QE = %0.2g',N_qubits,p_detect);
title(titlestr);

%%%
subplot(3,1,2);

p_MC_phi = ploterr(phi_ramsey,phi_est(:,1),[],phi_est(:,2),'o','hhxy',0);
set(p_MC_phi(1),'MarkerEdgeColor','r','MarkerFaceColor','w',...
    'LineWidth',1);
set(p_MC_phi(2),'Color','r','LineWidth',1);

hold on;
plot(phi_ramsey,phi_ramsey,'k--')

xlabel('true phase, $\phi$ (rad)');
ylabel('est. phase, $\phi_{\textrm{est}}$ (rad)');

%%%
subplot(3,1,3);

plot(phi_ramsey,phi_est(:,2),'r-');

xlabel('true phase, $\phi$ (rad)');
ylabel('uncertainty, $\Delta \phi_{\textrm{est}}$ (rad)');