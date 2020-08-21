%% Create SDE dipole data
clearvars,close all

%% Settings
gamma = 37;
x0 = 5.6;
epsilon = 1.9;
D = 100;
theta = 15;
tcorr = 0.6;

%% Define drift and noise functions
xDomain = -12:0.1:12;
driftFun = @(x) -gamma*(abs(x) - x0).*sign(x).*(1-exp(-abs(x)/epsilon));
driftLimit = @(x) -gamma*(abs(x) - x0).*sign(x);
noiseFun = @(x) 2*D + 0*x;

plotModel(xDomain,driftFun,driftLimit,noiseFun)

[X_PADM2M,t_PADM2M] = getPADM2M();

%% Simulating SDE
dt = 0.1;
dt_int_max = 0.2;
N = 2E5;

rng(1)
eta_ini = rand();
X_ini = x0*rand();
[t,X,eta] = OUcorr(dt,dt_int_max,N,driftFun,noiseFun,theta,eta_ini,X_ini);
Xconv = gaussConv(X,t,tcorr);

%% Time until reversal

tUntilReversal = timeUntilRev(t,Xconv);
tUntilReversal_PADM2M = timeUntilRev(t_PADM2M,X_PADM2M);

plotTimeUntilReversal(t_PADM2M,X_PADM2M,tUntilReversal_PADM2M,...
    t,Xconv,tUntilReversal)

%% Saving data

% Trim NaNs
save_t = t;
save_tUntilReversal = tUntilReversal;
save_X = Xconv;
save_t(isnan(tUntilReversal)) = [];
save_tUntilReversal(isnan(tUntilReversal)) = [];
save_X(isnan(tUntilReversal)) = [];

%save SDE_reversals.mat save_t save_X save_tUntilReversal

%% Functions 
function plotModel(xDomain,drift,driftLimit,noise)
% Load PADM2M
load PADM2M_recovered.mat xcentre f_est f_est_pct g_est g_est_pct

% Plotting
blue = [3,67,233]/255;
green = [21,176,26]/255;

figure('Position',[440   378   560*2   420])
subplot(1,2,1)
hold on,box on
fill([xcentre,fliplr(xcentre)],[f_est+f_est_pct(1,:),fliplr(f_est+f_est_pct(2,:))],...
    [0.9,0.9,0.9],'linestyle','none','facealpha',0.6)
leg.h1 = plot(xcentre,f_est+f_est_pct(1,:),'Color',blue,'LineStyle','--','LineWidth',2.0);
plot(xcentre,f_est+f_est_pct(2,:),'Color',blue,'LineStyle','--','LineWidth',2.0)
leg.h2 = plot(xcentre,f_est,'Color',blue,'LineStyle','-','LineWidth',3.0);
leg.h3 = plot(xDomain,drift(xDomain),'k-','LineWidth',2);
leg.h4 = plot(xDomain,driftLimit(xDomain),'k--','LineWidth',0.1);
title('Drift function')
xlabel('ADM state, \times 10^{22} Am^2')
ylabel('Drift, \times 10^{22} Am^2/kyrs')
legend([leg.h2,leg.h1,leg.h3,leg.h4],{'PADM2M','PADM2M 95% CI.','Model','Limit'})

subplot(1,2,2)
hold on,box on
fill([xcentre,fliplr(xcentre)],[g_est+g_est_pct(1,:),fliplr(g_est+g_est_pct(2,:))].^2,...
    [0.9,0.9,0.9],'linestyle','none','facealpha',0.6)
leg.h5 = plot(xcentre,(g_est+g_est_pct(1,:)).^2,'Color',green,'LineStyle','--','LineWidth',2.0);
plot(xcentre,(g_est+g_est_pct(2,:)).^2,'Color',green,'LineStyle','--','LineWidth',2.0)
leg.h6 = plot(xcentre,g_est.^2,'Color',green,'LineStyle','-','LineWidth',3.0);
leg.h7 = plot(xDomain,noise(xDomain),'k-','LineWidth',2);
title('Noise function')
xlabel('ADM state, \times 10^{22} Am^2')
ylabel('Noise, \times 10^{44} A^2m^4/kyrs')
legend([leg.h6,leg.h5,leg.h7],{'PADM2M','PADM2M 95% CI.','Model'})
xlim([0,inf])
ylim([0,inf])
end
function plotTimeUntilReversal(t_PADM2M,X_PADM2M,tUntilReversal_PADM2M,...
    t,Xconv,tUntilReversal)
figure
subplot(2,2,1)
plot(t_PADM2M,X_PADM2M,'k-')
ylabel('ADM, \times 10^{22} Am^2')
title('PADM2M')

subplot(2,2,3)
hold on,box on
plot(t_PADM2M,tUntilReversal_PADM2M,'k-')
plotNaNdata(t_PADM2M,tUntilReversal_PADM2M)
xlim([min(t_PADM2M),max(t_PADM2M)])
ylabel('Time until reversal, kyrs')
xlabel('Time, kyrs')

subplot(2,2,2)
hold on,box on
plot(t,Xconv,'k-')
ylabel('ADM, \times 10^{22} Am^2')
title('Simulation')

subplot(2,2,4)
hold on,box on
plot(t,tUntilReversal,'k-')
plotNaNdata(t,tUntilReversal)
xlim([min(t),max(t)])
ylabel('Time until reversal, kyrs')
xlabel('Time, kyrs')
end
function [X_PADM2M,t_PADM2M] = getPADM2M()
load PADM2M_reversing.mat X t
X_PADM2M = X{1};
t_PADM2M = 1000*t{1}';
end
function Xout = gaussConv(X,t,tcorr)
dt = t(2) - t(1);
tspan = 5*tcorr;
tkern = -tspan:dt:tspan;
kern = exp(-tkern.^2/(2*tcorr^2));
kernNorm = kern/sum(kern);
Xout = conv(X,kernNorm,'same');
end
function tUntilReversal = timeUntilRev(t,X)
reversalIndex = diff(sign(X)) ~= 0;
reversalTimes = t(reversalIndex);

tUntilReversal = t - reversalTimes';
tUntilReversal(tUntilReversal > 0) = nan;
tUntilReversal = -max(tUntilReversal);
end
function plotNaNdata(t,tUntilReversal_PADM2M)
nanLocations = isnan(tUntilReversal_PADM2M);
firstNan = find(nanLocations,1);
lastNan = find(nanLocations,1,'last');
t_first = t(firstNan);
t_last = t(lastNan);
yylim = ylim;
fill([t_first,t_last,t_last,t_first],yylim(2)*[0,0,1,1],'r',...
    'linestyle','none','facealpha',0.6)
end