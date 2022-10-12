clear

close all

% Define distributions epsilon1 & epsilon2
rng('default')  % For reproducibility

% Define your MCM parameters
n_trial = 100000;
con_step = 50;
un_pct = 5;

% Nominal Values for all variables
kuni=1.6483e+07;%/M/s
DGsp=7.3;%kcal/mol
DG_bp=1.7;%kcal/mol
DG_assoc=1.9;%kcal/mol
h=31;%toehold length
b=23;%branch migration length of base
R=8.314;%J/mol*K
T=310.15;%Kelvin


% True Values for Measured Outputs
keff_true = 364424.57;
n_true=0.99926;

% True Values for Measured Outputs (Calculation for validation purposes)
n=1/exp(DG_assoc/(R*T));
keff=(kuni*exp(-(DGsp-h*abs(DG_bp))/(R*T))*n)/(2*b);

% Normal distribution parameters @ 95% confidence
Xtrue = 0;
sigma = un_pct/200;     % Uncertainty divided by 2 (2-tailed test)


% MCM trials
kuni_mcm   = kuni  + sampleNormDist(Xtrue,kuni*sigma,n_trial);
DGsp_mcm  = DGsp  + sampleNormDist(Xtrue,DGsp*sigma,n_trial);
DG_bp_mcm = DG_bp + sampleNormDist(Xtrue,DG_bp*sigma,n_trial);
DG_assoc_mcm   = DG_assoc   + sampleNormDist(Xtrue,DG_assoc*sigma,n_trial);

n_mcm=1./exp(DG_assoc_mcm./(R.*T));
keff_mcm=(kuni_mcm.*(exp(-((DGsp_mcm-h*abs(DG_bp_mcm))./(R*T))).*n_mcm))./(2*b);

% Calculate Uncertainty
ukeff_mcm  = std(keff_mcm);        
un_mcm = std(n_mcm);  

Ukeff_xp = 2*ukeff_mcm;
Un_xp = 2*un_mcm;

Ukeff_keff = 2*ukeff_mcm/keff_true*100;
Un_n = 2*un_mcm/n_true*100;

% Plot the results
figure(1)
PlotFontSz = 16;
hist(keff_mcm,100), title('keff_m_c_m'),set(gca,'FontSize',PlotFontSz)

figure(2)
histfit(keff_mcm,100), title('keff_m_c_m'),set(gca,'FontSize',PlotFontSz)
hold on
plot([keff_true keff_true],[0 max(ylim)],'g-','LineWidth',3)
hold off

figure(3)
for i = 1:round(n_trial/con_step)
    uf_mcm_con(i)  = std(keff_mcm(1:i*con_step));
    num_con_checks = i;
end

plot([1: num_con_checks]*con_step/1000,uf_mcm_con), title('Convergence of s_m_c_m'),set(gca,'FontSize',PlotFontSz)
xlabel('Number of MC Evaluations [in 1000s]');
ylabel('s_m_c_m');

