%% Ramsey interferometry with non-ideal detection efficiency
% when measurement is normalised polarisation...
%
% DKS
% 2019-07-17

tic;

%% Config
%%% Ramsey 
% measurement scheme
p_up = @(phi) cos(phi/2).^2;        % probability to measure spin-UP
% y2phi_est = @(y) acos(y);
phi2ymean = @(phi) cos(phi);        % mean y predicted from phi (TODO: this is a GUESS)
dEydphi = @(phi) -sin(phi);         % average y vs. phi

% parameters
phi_ramsey = linspace(0,2*pi,1e2);       % phase
% N_qubits = 1e3;         % number of qubits
% p_detect = 0.1;         % detection efficiency

n_avg = 6.4;        % fixed
p_detect = 0.1;
N_qubits = round(n_avg/p_detect);


%%% Monte-carlo
N_mc = 1e3;          % repetition


%% Main 
y_mc = NaN(length(phi_ramsey),2);
phi_est = NaN(length(phi_ramsey),2);

parfor ii=1:length(phi_ramsey)
    [tphi_est, ty, N_spin_det, N_spin] = mc_ramsey(N_qubits,p_detect,phi_ramsey(ii),N_mc);
    
    % get statistics
    ty_mean = mean(ty,'omitnan');
    ty_std = std(ty,0,'omitnan');
    
    tphi_est_mean = mean(tphi_est,'omitnan');
    tphi_est_std = std(tphi_est,0,'omitnan');
    
    y_mc(ii,:) = [ty_mean, ty_std];
    phi_est(ii,:) = [tphi_est_mean, tphi_est_std];
end

dphi_pred = y_mc(:,2)./abs(dEydphi(phi_ramsey'));    % phase estimate from Dx = Dy/|dEy/dx|

%%% Residuals
% y vs theory
y_res = y_mc(:,1) - phi2ymean(phi_ramsey)';     % residual of y (sim) vs theory

% phi_est vs true phi
phi_ramsey_folded = abs(wrapToPi(phi_ramsey));      % phi wrapped in [0,pi]
phi_res = phi_est(:,1) - phi_ramsey_folded';        % residual



%% VIS
figname = sprintf('mcramsey_N%0.2e_qe%0.2e_S%0.2e',N_qubits,p_detect, N_mc);

H_mc = figure('Name',figname,'Units','centimeters','Position',[0,0,17.6,20]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Polarisation, y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phi vs. y
subplot(3,2,1);
hold on;

% MC
% p_MC_y = ploterr(phi_ramsey,y_mc(:,1),[],y_mc(:,2),'o','hhxy',0);
% set(p_MC_y(1),'MarkerEdgeColor','k','MarkerFaceColor','w',...
%     'LineWidth',1);
% set(p_MC_y(2),'Color','k','LineWidth',1);

mseb_opts.width = 1;
mseb_opts.edgestyle = ':';
mseb_opts.col = {'k'};
mseb_transp = 1;

p_MC_y = mseb(phi_ramsey,y_mc(:,1)',y_mc(:,2)',mseb_opts,mseb_transp);
p_MC_y(1).mainLine.DisplayName = 'data';

% theory
p_theory_phivsy = plot(phi_ramsey,phi2ymean(phi_ramsey),'b--',...
    'LineWidth',2,'DisplayName','$\cos \phi$');

xlabel('true phase, $\phi$ (rad)');
ylabel('polarisation, $y$');

box on;
xlim(minmax(phi_ramsey));
ylim([-1.5,1.5]);

lgd = legend(p_theory_phivsy,'Location','north');

%%% Residual
subplot(3,2,2);
hold on;

% MC
% p_MC_yres = ploterr(phi_ramsey,y_res,[],y_mc(:,2),'o','hhxy',0);
% set(p_MC_yres(1),'MarkerEdgeColor','k','MarkerFaceColor','w',...
%     'LineWidth',1);
% set(p_MC_yres(2),'Color','k','LineWidth',1);

mseb_opts.width = 1;
mseb_opts.edgestyle = ':';
mseb_opts.col = {'k'};
mseb_transp = 1;

p_MC_yres = mseb(phi_ramsey,y_res',y_mc(:,2)',mseb_opts,mseb_transp);
p_MC_yres(1).mainLine.DisplayName = 'data';

xlabel('true phase, $\phi$ (rad)');
ylabel('$y - \cos\phi$');

box on;
xlim(minmax(phi_ramsey));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phi vs. phi_est
subplot(3,2,3);
hold on;

% MC
% p_MC_phi = ploterr(phi_ramsey,phi_est(:,1),[],phi_est(:,2),'o','hhxy',0);
% set(p_MC_phi(1),'MarkerEdgeColor','r','MarkerFaceColor','w',...
%     'LineWidth',1);
% set(p_MC_phi(2),'Color','r','LineWidth',1);

mseb_opts.width = 1;
mseb_opts.edgestyle = ':';
mseb_opts.col = {'r'};
mseb_transp = 0;

p_MC_phi = mseb(phi_ramsey,phi_est(:,1)',phi_est(:,2)',mseb_opts,mseb_transp);
p_MC_phi(1).mainLine.DisplayName = 'data';


hold on;
p_unbiasedest1 = plot(phi_ramsey,phi_ramsey,'k-','DisplayName','$\phi_\textrm{est} = \phi$');
p_unbiasedest2 = plot(phi_ramsey,phi_ramsey_folded,'k--','DisplayName','folded');

xlabel('true phase, $\phi$ (rad)');
ylabel('est. phase, $\phi_{\textrm{est}}$ (rad)');

box on;
xlim(minmax(phi_ramsey));
ylim([-pi/4,pi+pi/4])

lgd = legend([p_unbiasedest1,p_unbiasedest2],'Location','south');


%%% Residual
subplot(3,2,4);
hold on;

% MC
% p_MC_phires = ploterr(phi_ramsey,phi_res(:,1),[],phi_est(:,2),'o','hhxy',0);
% set(p_MC_phires(1),'MarkerEdgeColor','r','MarkerFaceColor','w',...
%     'LineWidth',1);
% set(p_MC_phires(2),'Color','r','LineWidth',1);

mseb_opts.width = 1;
mseb_opts.edgestyle = ':';
mseb_opts.col = {'r'};
mseb_transp = 0;

p_MC_phires = mseb(phi_ramsey,phi_res(:,1)',phi_est(:,2)',mseb_opts,mseb_transp);
p_MC_phires(1).mainLine.DisplayName = 'data';


xlabel('true phase, $\phi$ (rad)');
ylabel('$\phi_{\textrm{est}} - \phi$');

box on;
xlim(minmax(phi_ramsey));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phi vs. sdev of phi_est
subplot(3,2,[5,6]);
hold on;

% MC
p_phiest = plot(phi_ramsey,phi_est(:,2),'r-','DisplayName','$\Delta \phi_\textrm{est}$');
p_phipred = plot(phi_ramsey,dphi_pred,'b-.','DisplayName','$\Delta y/|d\bar{y}/d\phi|$');

% quantum limits
p_detlim = plot(phi_ramsey,1/sqrt(p_detect*N_qubits)*ones(size(phi_ramsey)),'k-','DisplayName','DL');
p_sql = plot(phi_ramsey,1/sqrt(N_qubits)*ones(size(phi_ramsey)),'k:','DisplayName','SQL');
p_hl = plot(phi_ramsey,1/N_qubits*ones(size(phi_ramsey)),'k--','DisplayName','HL');

xlabel('true phase, $\phi$ (rad)');
ylabel('phase uncertainty (rad)');

box on;
xlim(minmax(phi_ramsey));

legend([p_phiest,p_phipred,p_detlim,p_sql,p_hl],'Location','eastoutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global annotation

titlestr = sprintf('N(qubits) = %d; QE = %0.2g; samples = %0.1E',N_qubits,p_detect, N_mc);
mtit(titlestr);

%% end
toc;