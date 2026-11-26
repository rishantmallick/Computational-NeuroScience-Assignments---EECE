%% =========================================================
% COMPUTATIONAL NEUROSCIENCE - PROJECT III
% QUESTIONS 1-6 AND OPTION A: MAXIMALLY INFORMATIVE DIMENSIONS
%
% Input variables:
%   All_Spike_Times : 4 x 50 cell array
%   Stimulus        : 20000 samples at 1 kHz
%
% Training interval: 0-15 seconds
% Testing interval : 15-20 seconds
%
% IMPORTANT:
% The original input MAT file is loaded but no output MAT file
% is created by this script.
%% =========================================================

clear;
close all;
clc;
rng(7,'twister');

%% =========================================================
% LOAD DATA AND INITIALIZE PARAMETERS
%% =========================================================

[fileName,filePath] = uigetfile('*.mat', ...
    'Select the project data MAT file');

if isequal(fileName,0)
    error('No MAT file selected.');
end

D = load(fullfile(filePath,fileName));

if ~isfield(D,'All_Spike_Times') || ~isfield(D,'Stimulus')
    error('The MAT file must contain All_Spike_Times and Stimulus.');
end

All_Spike_Times = D.All_Spike_Times;
Stimulus = double(D.Stimulus(:));

fs = 1000;
dt = 1/fs;

T = numel(Stimulus);
time = (0:T-1)'*dt;

nNeurons = size(All_Spike_Times,1);
nRepeats = size(All_Spike_Times,2);

trainEnd = min(15000,T);
testStart = min(15001,T+1);

if T ~= 20000
    warning('Expected 20000 stimulus samples; found %d.',T);
end

if nNeurons ~= 4 || nRepeats ~= 50
    warning('Expected 4 neurons and 50 repetitions.');
end

if trainEnd < 100 || testStart > T
    error('Insufficient data for the train/test split.');
end

%% =========================================================
% PREPROCESSING: CONVERT SPIKE TIMES TO 1-ms BINS
%% =========================================================

SpikeCounts = zeros(T,nNeurons,nRepeats);

for n = 1:nNeurons
    for rep = 1:nRepeats

        st = double(All_Spike_Times{n,rep}(:));

        if isempty(st)
            continue;
        end

        idx = floor(st*fs)+1;
        idx = idx(idx>=1 & idx<=T);

        if ~isempty(idx)
            SpikeCounts(:,n,rep) = ...
                accumarray(idx,1,[T 1],@sum,0);
        end
    end
end

% PSTH in spikes/second
PSTH = reshape(sum(SpikeCounts,3)/(nRepeats*dt), ...
    T,nNeurons);

% Last five seconds: held-out reference responses
PSTH_test = PSTH(testStart:end,:);

%% =========================================================
% QUESTION 1: STIMULUS NATURE AND AUTOCORRELATION
%% =========================================================

fprintf('\nQUESTION 1: STIMULUS NATURE\n');

s0 = Stimulus-mean(Stimulus);

fprintf('Stimulus mean = %.5f dB\n',mean(Stimulus));
fprintf('Stimulus standard deviation = %.5f dB\n',std(Stimulus));

figure('Name','Question 1: Stimulus','Color','w');

subplot(2,2,1);

plot(time,Stimulus,'k-');
xlabel('Time (s)');
ylabel('Stimulus envelope (dB)');
title('Stimulus envelope');
grid on;

subplot(2,2,2);

histogram(Stimulus,50,'Normalization','pdf');
hold on;

xx = linspace(min(Stimulus),max(Stimulus),250);
mu = mean(Stimulus);
sd = std(Stimulus);

if sd > 0
    pdfGaussian = exp(-0.5*((xx-mu)/sd).^2) ...
        /(sd*sqrt(2*pi));

    plot(xx,pdfGaussian,'r-','LineWidth',1.5);
end

xlabel('Stimulus amplitude (dB)');
ylabel('Probability density');
title('Amplitude distribution and Gaussian fit');

% Empirical autocorrelation at lags in milliseconds
maxLag = 100;

[ac,lags] = xcorr(s0,maxLag,'biased');

ac = ac(lags>=0);
lags = lags(lags>=0);

if ac(1) ~= 0
    ac = ac/ac(1);
end

subplot(2,2,3);

stem(lags,ac,'Marker','none');
xlabel('Lag (ms)');
ylabel('Normalized autocorrelation');
title('Autocorrelation: 0-100 ms');
grid on;

subplot(2,2,4);

numShow = min(21,numel(ac));

stem(lags(1:numShow),ac(1:numShow),'Marker','none');
xlabel('Lag (ms)');
ylabel('Normalized autocorrelation');
title('Autocorrelation near zero lag');
grid on;

StimulusAutocorrelation = table( ...
    lags(:),ac(:), ...
    'VariableNames',{'Lag_ms','NormalizedAutocorrelation'});

writetable(StimulusAutocorrelation, ...
    'Q1_Stimulus_Autocorrelation.csv');

%% =========================================================
% QUESTION 2: PSTH AND MEAN FIRING RATE
%% =========================================================

fprintf('\nQUESTION 2: PSTH AND MEAN FIRING RATE\n');

figure('Name','Question 2: PSTHs','Color','w');

for n = 1:nNeurons

    subplot(nNeurons,1,n);

    plot(time,PSTH(:,n),'k-','LineWidth',0.8);
    hold on;

    xline(15,'r--','Train/test split');

    xlim([0 time(end)]);

    ylabel(sprintf('Neuron %d\nspikes/s',n));
    grid on;

    if n == 1
        title('PSTH at 1 ms resolution');
    end

    if n == nNeurons
        xlabel('Time (s)');
    end
end

meanRateNeuron = zeros(nNeurons,1);

for n = 1:nNeurons
    meanRateNeuron(n) = ...
        sum(SpikeCounts(:,n,:),'all')/(nRepeats*T*dt);
end

MeanFiringRates = table( ...
    (1:nNeurons)',meanRateNeuron, ...
    'VariableNames',{'Neuron','MeanFiringRate_spikes_s'});

disp(MeanFiringRates);

%% =========================================================
% QUESTION 3: POISSON OR NON-POISSON
%% =========================================================

fprintf('\nQUESTION 3: POISSON OR NON-POISSON\n');

binWidthsMs = [10 20 50 100 200 500];

% Columns: Fano ratio, mean count, variance, pooled Fano factor
PoissonStats = nan(nNeurons,numel(binWidthsMs),4);

for n = 1:nNeurons

    figure('Name',sprintf('Question 3: Neuron %d',n), ...
        'Color','w');

    % Inter-spike intervals pooled over repetitions.
    stAll = [];

    for rep = 1:nRepeats
        stAll = [stAll; ...
            double(All_Spike_Times{n,rep}(:))]; %#ok<AGROW>
    end

    stAll = sort(stAll);
    isi = diff(stAll);
    isi = isi(isi>0);

    subplot(2,2,1);

    if ~isempty(isi)
        histogram(isi,50,'Normalization','pdf');
    end

    xlabel('ISI (s)');
    ylabel('Probability density');
    title(sprintf('Neuron %d: ISI histogram',n));
    grid on;

    cvISI = NaN;

    if numel(isi)>1 && mean(isi)>0
        cvISI = std(isi)/mean(isi);
    end

    subplot(2,2,2);
    axis off;

    text(0.05,0.8,sprintf('ISI CV = %.3f',cvISI), ...
        'FontSize',11);

    text(0.05,0.6, ...
        'CV near 1 is consistent with exponential ISIs.', ...
        'FontSize',10);

    text(0.05,0.4, ...
        'This alone does not prove Poisson firing.', ...
        'FontSize',10);

    meanCounts = zeros(size(binWidthsMs));
    varCounts = zeros(size(binWidthsMs));

    for b = 1:numel(binWidthsMs)

        bw = binWidthsMs(b);
        nBins = floor(T/bw);

        counts = zeros(nBins,nRepeats);

        for rep = 1:nRepeats

            c = SpikeCounts(1:nBins*bw,n,rep);

            counts(:,rep) = ...
                sum(reshape(c,bw,nBins),1)';
        end

        % Mean and variance across repetitions at each
        % corresponding time bin.
        m = mean(counts,2);
        v = var(counts,0,2);

        meanCounts(b) = mean(m);
        varCounts(b) = mean(v);

        PoissonStats(n,b,1) = mean(v./max(m,eps));
        PoissonStats(n,b,2) = mean(m);
        PoissonStats(n,b,3) = mean(v);

        PoissonStats(n,b,4) = ...
            var(counts(:))/max(mean(counts(:)),eps);
    end

    subplot(2,2,3);

    loglog(meanCounts,varCounts,'o-','LineWidth',1.2);
    hold on;

    limLow = max(min([meanCounts(:);varCounts(:)]),1e-5);
    limHigh = max([meanCounts(:);varCounts(:);1]);

    plot([limLow limHigh],[limLow limHigh],'k--');

    xlabel('Mean spike count');
    ylabel('Variance of spike count');
    title('Mean-variance comparison');

    legend('Observed','Variance = mean','Location','best');
    grid on;

    subplot(2,2,4);

    plot(binWidthsMs,PoissonStats(n,:,1), ...
        'o-','LineWidth',1.2);

    hold on;
    yline(1,'k--');

    xlabel('Bin width (ms)');
    ylabel('Mean variance/mean');
    title('Fano ratio');
    grid on;
end

NeuronColumn = [];
BinColumn = [];
FanoColumn = [];
MeanColumn = [];
VarianceColumn = [];

for n = 1:nNeurons
    for b = 1:numel(binWidthsMs)

        NeuronColumn(end+1,1) = n; %#ok<SAGROW>
        BinColumn(end+1,1) = binWidthsMs(b); %#ok<SAGROW>

        FanoColumn(end+1,1) = PoissonStats(n,b,1); %#ok<SAGROW>
        MeanColumn(end+1,1) = PoissonStats(n,b,2); %#ok<SAGROW>
        VarianceColumn(end+1,1) = PoissonStats(n,b,3); %#ok<SAGROW>
    end
end

PoissonTable = table(NeuronColumn,BinColumn,FanoColumn, ...
    MeanColumn,VarianceColumn, ...
    'VariableNames',{'Neuron','BinWidth_ms','MeanFanoRatio', ...
    'MeanCount','VarianceCount'});

writetable(PoissonTable,'Q3_Poisson_Statistics.csv');

%% =========================================================
% QUESTION 4: SPIKE-TRIGGERED AVERAGE (STA)
% AND STIMULUS COVARIANCE CORRECTION
%% =========================================================

fprintf('\nQUESTION 4: STA AND LINEAR FILTERS\n');

staLength = 50;

X50 = makeHistory(Stimulus,staLength);

trainIdx = (staLength+1:trainEnd)';
testIdx = (testStart:T)';

Xtr50 = X50(trainIdx,:);
Xte50 = X50(testIdx,:);

Filters = struct();

figure('Name','Question 4: STA filters','Color','w');

for n = 1:nNeurons

    ytr = PSTH(trainIdx,n);

    % Spike-triggered average using spike counts as weights.
    w = ytr*nRepeats*dt;

    if sum(w)>0
        sta = (Xtr50'*w)/sum(w)-mean(Xtr50,1)';
    else
        sta = zeros(staLength,1);
    end

    % Empirical stimulus covariance matrix.
    Xc = Xtr50-mean(Xtr50,1);

    C = (Xc'*Xc)/size(Xc,1);

    % Ridge regularization for numerical stability.
    ridge = 1e-3*trace(C)/staLength;

    hCorr = (C+ridge*eye(staLength))\sta;

    Filters(n).STA = sta;
    Filters(n).hCorrected = hCorr;
    Filters(n).stimCovariance = C;

    subplot(nNeurons,1,n);

    plot(1:staLength,sta,'b-','LineWidth',1.1);
    hold on;

    plot(1:staLength,hCorr,'r-','LineWidth',1.1);

    ylabel(sprintf('Neuron %d',n));
    grid on;

    if n == 1
        title('STA and covariance-corrected filter');
    end

    if n == nNeurons
        xlabel('Lag index (ms)');
    end
end

legend({'STA','Covariance-corrected filter'}, ...
    'Location','best');

%% =========================================================
% QUESTION 5: DETERMINING THE OUTPUT NONLINEARITY
%% =========================================================

fprintf('\nQUESTION 5: OUTPUT NONLINEARITY\n');

Nonlinearity = struct();

Pred5 = nan(numel(testIdx),nNeurons);
R2_5 = nan(nNeurons,1);

figure('Name','Question 5: Output nonlinearities','Color','w');

for n = 1:nNeurons

    h = Filters(n).hCorrected;

    ytr = PSTH(trainIdx,n);
    yte = PSTH(testIdx,n);

    % Linear projection.
    ylinTr = Xtr50*h;
    ylinTe = Xte50*h;

    % Fit output nonlinearity on training data only.
    theta = fitSigmoid1D(ylinTr,ytr);

    pTr = sigmoid1D(ylinTr,theta);
    pTe = sigmoid1D(ylinTe,theta);

    Nonlinearity(n).theta = theta;
    Nonlinearity(n).linearTrain = ylinTr;
    Nonlinearity(n).predTrain = pTr;
    Nonlinearity(n).predTest = pTe;

    Pred5(:,n) = pTe;
    R2_5(n) = corrR2(yte,pTe);

    subplot(2,2,n);

    scatter(ylinTr,ytr,8,'filled');
    hold on;

    [xs,~] = sort(ylinTr);

    plot(xs,sigmoid1D(xs,theta),'r-','LineWidth',1.5);

    xlabel('Linear projection');
    ylabel('Measured PSTH (spikes/s)');

    title(sprintf('Neuron %d: R^2 = %.3f',n,R2_5(n)));
    grid on;
end

figure('Name','Question 5: Held-out predictions','Color','w');

for n = 1:nNeurons

    subplot(nNeurons,1,n);

    plot(time(testIdx),PSTH(testIdx,n),'k-');
    hold on;

    plot(time(testIdx),Pred5(:,n),'r-');

    ylabel(sprintf('Neuron %d\nspikes/s',n));
    grid on;

    if n == 1
        title('Held-out prediction: 15-20 seconds');
    end

    if n == nNeurons
        xlabel('Time (s)');
    end
end

%% =========================================================
% QUESTION 6: PREDICTION PERFORMANCE AND PRUNING
%% =========================================================

fprintf('\nQUESTION 6: PREDICTION AND PRUNING\n');

PrunedFilters = struct();
PruneCurve = cell(nNeurons,1);

Pred6 = nan(numel(testIdx),nNeurons);
R2_6 = nan(nNeurons,1);

figure('Name','Question 6: Prediction scatter plots','Color','w');

for n = 1:nNeurons

    ytr = PSTH(trainIdx,n);
    yte = PSTH(testIdx,n);

    h = Filters(n).hCorrected;

    [hp,curve] = pruneFilter(h,Xtr50,Xte50,ytr,yte);

    thp = fitSigmoid1D(Xtr50*hp,ytr);
    pp = sigmoid1D(Xte50*hp,thp);

    r2p = corrR2(yte,pp);

    PrunedFilters(n).h = hp;
    PrunedFilters(n).theta = thp;
    PrunedFilters(n).R2 = r2p;

    PruneCurve{n} = curve;

    Pred6(:,n) = pp;
    R2_6(n) = r2p;

    subplot(2,2,n);

    scatter(yte,pp,8,'filled');
    hold on;

    addUnityLine(yte,pp);

    xlabel('Measured PSTH (spikes/s)');
    ylabel('Predicted firing rate');

    title(sprintf('Neuron %d: pruned R^2 = %.3f',n,r2p));
    grid on;
end

figure('Name','Question 6: Pruning curves','Color','w');

for n = 1:nNeurons

    subplot(2,2,n);

    plot(0:numel(PruneCurve{n})-1, ...
        PruneCurve{n},'o-');

    xlabel('Number of coefficients pruned');
    ylabel('Held-out R^2');

    title(sprintf('Neuron %d',n));
    grid on;
end

figure('Name','Question 6: Final filters and FFT','Color','w');

for n = 1:nNeurons

    h = PrunedFilters(n).h;

    nfft = 1024;

    H = abs(fft(h,nfft));
    H = H(1:floor(nfft/2)+1);

    f = (0:floor(nfft/2))'*fs/nfft;

    subplot(nNeurons,2,2*n-1);

    plot(1:staLength,h);

    ylabel(sprintf('Neuron %d weight',n));
    grid on;

    if n == 1
        title('Pruned filter');
    end

    if n == nNeurons
        xlabel('Lag (ms)');
    end

    subplot(nNeurons,2,2*n);

    plot(f,H);

    ylabel('|FFT|');
    grid on;

    if n == 1
        title('Filter frequency response');
    end

    if n == nNeurons
        xlabel('Frequency (Hz)');
    end
end

%% =========================================================
% OPTION A: MAXIMALLY INFORMATIVE DIMENSIONS
% PART A1: ONE-VECTOR MID
%% =========================================================

fprintf('\nOPTION A1: ONE-VECTOR MID\n');

midLength = 100;

X100 = makeHistory(Stimulus,midLength);

midTrainIdx = (midLength+1:trainEnd)';
midTestIdx = (testStart:T)';

Xtr = X100(midTrainIdx,:);
Xte = X100(midTestIdx,:);

MID = struct();

% Columns: 1-vector, pruned 1-vector, 2-vector, pruned 2-vector.
MID_R2 = nan(nNeurons,4);

nRestarts = 2;
nMIDBins = 8;
maxMIDIter = 30;

for n = 1:nNeurons

    fprintf('MID optimization: neuron %d\n',n);

    ytr = PSTH(midTrainIdx,n);
    yte = PSTH(midTestIdx,n);

    w = max(ytr*nRepeats*dt,0);

    % Subsample training time bins to reduce optimization cost.
    midStride = 5;

    Xopt = Xtr(1:midStride:end,:);
    wopt = w(1:midStride:end);

    if sum(w)<=0
        warning('No training spikes for neuron %d.',n);
        continue;
    end

    % STA initialization.
    sta = (Xtr'*w)/sum(w)-mean(Xtr,1)';

    if norm(sta)<eps
        sta = randn(midLength,1);
    end

    v0 = sta/norm(sta);

    % Search from multiple initializations.
    bestMI = -Inf;
    bestV = v0;

    restartMI = nan(nRestarts,1);

    for rr = 1:nRestarts

        if rr == 1
            vStart = v0;
        else
            vStart = v0+0.2*randn(midLength,1);
            vStart = vStart/norm(vStart);
        end

        [vc,mi] = optimizeMID1( ...
            Xopt,wopt,vStart,nMIDBins,maxMIDIter);

        restartMI(rr) = mi;

        if mi > bestMI
            bestMI = mi;
            bestV = vc;
        end
    end

    v1 = bestV;

    xtr = Xtr*v1;
    xte = Xte*v1;

    % Fit one-dimensional output nonlinearity.
    th1 = fitSigmoid1D(xtr,ytr);

    p1 = sigmoid1D(xte,th1);
    score1 = corrR2(yte,p1);

    % Prune small-magnitude vector coefficients.
    [v1p,curve1] = pruneVectorMID(v1,Xtr,Xte,ytr,yte);

    th1p = fitSigmoid1D(Xtr*v1p,ytr);

    p1p = sigmoid1D(Xte*v1p,th1p);

    score1p = corrR2(yte,p1p);

    MID(n).STA = sta;
    MID(n).v1 = v1;
    MID(n).MI1 = bestMI;
    MID(n).MI1Restarts = restartMI;

    MID(n).theta1 = th1;
    MID(n).pred1 = p1;
    MID(n).R2_1 = score1;

    MID(n).v1Pruned = v1p;
    MID(n).pred1Pruned = p1p;
    MID(n).R2_1Pruned = score1p;
    MID(n).pruneCurve1 = curve1;

    MID_R2(n,1:2) = [score1 score1p];

    %% =====================================================
    % OPTION A2: TWO-VECTOR MID
    %% =====================================================

    fprintf('Two-vector MID: neuron %d\n',n);

    % Initialize second vector orthogonal to the first.
    v2 = randn(midLength,1);

    v2 = v2-v1*(v1'*v2);

    if norm(v2)<1e-8
        Q = null(v1');
        v2 = Q(:,1);
    end

    v2 = v2/norm(v2);

    bestMI2 = -Inf;
    bestV2 = [v1 v2];

    restartMI2 = nan(nRestarts,1);

    for rr = 1:nRestarts

        if rr == 1

            Vstart = [v1 v2];

        else

            a = v1+0.05*randn(midLength,1);
            a = a/norm(a);

            b = randn(midLength,1);
            b = b-a*(a'*b);
            b = b/max(norm(b),eps);

            Vstart = [a b];
        end

        [Vc,mi] = optimizeMID2( ...
            Xopt,wopt,Vstart,nMIDBins,maxMIDIter);

        restartMI2(rr) = mi;

        if mi > bestMI2
            bestMI2 = mi;
            bestV2 = Vc;
        end
    end

    V = bestV2;

    XYtr = Xtr*V;
    XYte = Xte*V;

    % Fit a two-dimensional product-of-sigmoids nonlinearity.
    th2 = fitSigmoid2D(XYtr(:,1),XYtr(:,2),ytr);

    p2 = sigmoid2D(XYte(:,1),XYte(:,2),th2);

    score2 = corrR2(yte,p2);

    % Prune coefficients in the two-vector model.
    [Vp,curve2] = pruneTwoMID(V,Xtr,Xte,ytr,yte);

    XYtrp = Xtr*Vp;
    XYtep = Xte*Vp;

    th2p = fitSigmoid2D(XYtrp(:,1),XYtrp(:,2),ytr);

    p2p = sigmoid2D(XYtep(:,1),XYtep(:,2),th2p);

    score2p = corrR2(yte,p2p);

    MID(n).V2 = V;
    MID(n).MI2 = bestMI2;
    MID(n).MI2Restarts = restartMI2;

    MID(n).theta2 = th2;
    MID(n).pred2 = p2;
    MID(n).R2_2 = score2;

    MID(n).V2Pruned = Vp;
    MID(n).theta2Pruned = th2p;
    MID(n).pred2Pruned = p2p;
    MID(n).R2_2Pruned = score2p;
    MID(n).pruneCurve2 = curve2;

    MID_R2(n,3:4) = [score2 score2p];

    %% MID COMPARISON FIGURES

    figure('Name',sprintf('MID neuron %d',n),'Color','w');

    subplot(2,2,1);

    scatter(yte,p1,8,'filled');
    hold on;

    addUnityLine(yte,p1);

    xlabel('Measured PSTH');
    ylabel('Predicted firing rate');

    title(sprintf('One vector R^2 = %.3f',score1));
    grid on;

    subplot(2,2,2);

    scatter(yte,p2,8,'filled');
    hold on;

    addUnityLine(yte,p2);

    xlabel('Measured PSTH');
    ylabel('Predicted firing rate');

    title(sprintf('Two vectors R^2 = %.3f',score2));
    grid on;

    subplot(2,2,3);

    plot(time(midTestIdx),yte,'k-');
    hold on;

    plot(time(midTestIdx),p1,'b-');
    plot(time(midTestIdx),p2,'r-');

    xlabel('Time (s)');
    ylabel('Firing rate (spikes/s)');

    legend('Observed','One vector','Two vectors', ...
        'Location','best');

    title('Held-out prediction');
    grid on;

    subplot(2,2,4);

    plot(1:midLength,v1,'b-');
    hold on;

    plot(1:midLength,V(:,1),'r--');
    plot(1:midLength,V(:,2),'k-');

    xlabel('Lag (ms)');
    ylabel('Filter weight');

    legend('One-vector MID','Two-vector v_1','Two-vector v_2', ...
        'Location','best');

    title('MID filters');
    grid on;

    fprintf(['Neuron %d: MI1=%.4f, R2 one-vector=%.3f, ' ...
        'pruned=%.3f, two-vector=%.3f, pruned=%.3f\n'], ...
        n,bestMI,score1,score1p,score2,score2p);
end

%% =========================================================
% OPTION A: SUMMARY OF ONE-VECTOR AND TWO-VECTOR PERFORMANCE
%% =========================================================

MIDComparison = table( ...
    (1:nNeurons)',MID_R2(:,1),MID_R2(:,2), ...
    MID_R2(:,3),MID_R2(:,4), ...
    'VariableNames',{'Neuron','R2_1Vector','R2_1VectorPruned', ...
    'R2_2Vectors','R2_2VectorsPruned'});

disp(MIDComparison);

writetable(MIDComparison,'OptionA_MID_Comparison.csv');

%% =========================================================
% LOCAL HELPER FUNCTIONS
% These functions support the calculations above.
%% =========================================================

function X = makeHistory(s,L)

    s = s(:);
    T = numel(s);

    X = zeros(T,L);

    for t = L+1:T
        X(t,:) = s(t-L:t-1)';
    end
end

function theta = fitSigmoid1D(x,y)

    x = x(:);
    y = max(y(:),0);

    A0 = max(prctile(y,95),1);
    a0 = 1/max(std(x),1e-3);
    b0 = median(x);
    c0 = max(min(y),0);

    th0 = [log(A0);log(a0);b0;c0];

    obj = @(th) mean((sigmoid1D(x,th)-y).^2) ...
        +1e-8*sum(th.^2);

    opts = optimset('Display','off', ...
        'MaxIter',1200,'MaxFunEvals',3000);

    theta = fminsearch(obj,th0,opts);
end

function yhat = sigmoid1D(x,th)

    A = exp(min(th(1),20));
    a = exp(min(th(2),20));
    b = th(3);
    c = max(th(4),0);

    z = max(min(-a*(x-b),60),-60);

    yhat = max(c+A./(1+exp(z)),0);
end

function theta = fitSigmoid2D(x,y,r)

    x = x(:);
    y = y(:);
    r = max(r(:),0);

    A0 = max(prctile(r,95),1);
    a0 = 1/max(std(x),1e-3);
    b0 = median(x);

    c0 = 1/max(std(y),1e-3);
    d0 = median(y);

    base0 = max(min(r),0.01);

    th0 = [log(A0);log(a0);b0;log(c0);d0;log(base0)];

    obj = @(th) mean((sigmoid2D(x,y,th)-r).^2) ...
        +1e-8*sum(th.^2);

    opts = optimset('Display','off', ...
        'MaxIter',1800,'MaxFunEvals',4000);

    theta = fminsearch(obj,th0,opts);
end

function r = sigmoid2D(x,y,th)

    A = exp(min(th(1),20));
    a = exp(min(th(2),20));
    b = th(3);

    c = exp(min(th(4),20));
    d = th(5);

    base = exp(min(th(6),20));

    zx = max(min(-a*(x-b),60),-60);
    zy = max(min(-c*(y-d),60),-60);

    r = max(base+A./(1+exp(zx))./(1+exp(zy)),0);
end

function score = corrR2(a,b)

    a = a(:);
    b = b(:);

    ok = isfinite(a)&isfinite(b);

    a = a(ok);
    b = b(ok);

    if numel(a)<2 || std(a)<eps || std(b)<eps
        score = NaN;
    else
        C = corrcoef(a,b);
        score = C(1,2)^2;
    end
end

function [hp,curve] = pruneFilter(h,Xtr,Xte,ytr,yte)

    hp = h(:);

    [~,ord] = sort(abs(hp),'ascend');

    curve = nan(numel(hp)+1,1);

    th = fitSigmoid1D(Xtr*hp,ytr);

    curve(1) = corrR2(yte,sigmoid1D(Xte*hp,th));

    best = curve(1);

    for k = 1:numel(ord)

        trial = hp;
        trial(ord(k)) = 0;

        th = fitSigmoid1D(Xtr*trial,ytr);

        sc = corrR2(yte,sigmoid1D(Xte*trial,th));

        curve(k+1) = sc;

        if isnan(sc) || (~isnan(best) && sc<best-0.01)
            curve = curve(1:k+1);
            break;
        end

        hp = trial;

        if isnan(best) || sc>best
            best = sc;
        end
    end

    if norm(hp)>0
        hp = hp/norm(hp);
    end
end

function score = midObjective1(X,w,v,nBins)

    x = X*v;
    w = max(w(:),0);

    if sum(w)<=0
        score = 0;
        return;
    end

    edges = quantile(x,linspace(0,1,nBins+1));

    edges(1) = -Inf;
    edges(end) = Inf;

    edges = unique(edges);

    if numel(edges)<3
        score = 0;
        return;
    end

    b = discretize(x,edges);

    nb = numel(edges)-1;

    pX = accumarray(b,1,[nb 1],@sum,0);
    pX = pX/sum(pX);

    pXs = accumarray(b,w,[nb 1],@sum,0);
    pXs = pXs/max(sum(pXs),eps);

    good = pXs>0 & pX>0;

    score = sum(pXs(good).*log2(pXs(good)./pX(good)));

    score = max(score,0);
end

function [vbest,mibest] = optimizeMID1(X,w,v0,nBins,maxIter)

    v = v0(:)/max(norm(v0),eps);

    mibest = midObjective1(X,w,v,nBins);

    step = 0.20;

    for it = 1:maxIter

        improved = false;

        for signVal = [1 -1]

            d = randn(size(v));

            d = d-v*(v'*d);
            d = d/max(norm(d),eps);

            cand = v+signVal*step*d;
            cand = cand/max(norm(cand),eps);

            mi = midObjective1(X,w,cand,nBins);

            if mi>mibest+1e-7
                v = cand;
                mibest = mi;
                improved = true;
            end
        end

        if ~improved
            step = step*0.85;
        end

        if step<1e-3
            break;
        end
    end

    vbest = v;
end

function score = midObjective2(X,w,V,nBins)

    XY = X*V;

    x = XY(:,1);
    y = XY(:,2);

    w = max(w(:),0);

    if sum(w)<=0
        score = 0;
        return;
    end

    ex = quantile(x,linspace(0,1,nBins+1));
    ey = quantile(y,linspace(0,1,nBins+1));

    ex(1) = -Inf;
    ex(end) = Inf;

    ey(1) = -Inf;
    ey(end) = Inf;

    ex = unique(ex);
    ey = unique(ey);

    nx = numel(ex)-1;
    ny = numel(ey)-1;

    if nx<1 || ny<1
        score = 0;
        return;
    end

    bx = discretize(x,ex);
    by = discretize(y,ey);

    valid = ~isnan(bx)&~isnan(by);

    bx = bx(valid);
    by = by(valid);
    w = w(valid);

    ind = sub2ind([nx ny],bx,by);

    pXY = accumarray(ind,1,[nx*ny 1],@sum,0);
    pXY = pXY/sum(pXY);

    pXYs = accumarray(ind,w,[nx*ny 1],@sum,0);
    pXYs = pXYs/max(sum(pXYs),eps);

    good = pXYs>0 & pXY>0;

    score = sum(pXYs(good).*log2(pXYs(good)./pXY(good)));

    score = max(score,0);
end

function V = orthonormalize(z,d)

    V = reshape(z,d,2);

    V(:,1) = V(:,1)/max(norm(V(:,1)),eps);

    V(:,2) = V(:,2)-V(:,1)*(V(:,1)'*V(:,2));

    if norm(V(:,2))<1e-8
        q = randn(d,1);
        q = q-V(:,1)*(V(:,1)'*q);
        V(:,2) = q;
    end

    V(:,2) = V(:,2)/max(norm(V(:,2)),eps);
end

function [Vbest,mibest] = optimizeMID2(X,w,V0,nBins,maxIter)

    V = orthonormalize(V0(:),size(X,2));

    mibest = midObjective2(X,w,V,nBins);

    step = 0.12;

    for it = 1:maxIter

        improved = false;

        for col = 1:2

            d = randn(size(V,1),1);

            d = d-V(:,1)*(V(:,1)'*d) ...
                -V(:,2)*(V(:,2)'*d);

            d = d/max(norm(d),eps);

            for sg = [1 -1]

                C = V;
                C(:,col) = C(:,col)+sg*step*d;

                C = orthonormalize(C(:),size(X,2));

                mi = midObjective2(X,w,C,nBins);

                if mi>mibest+1e-7
                    V = C;
                    mibest = mi;
                    improved = true;
                end
            end
        end

        if ~improved
            step = step*0.85;
        end

        if step<1e-3
            break;
        end
    end

    Vbest = V;
end

function [vp,curve] = pruneVectorMID(v,Xtr,Xte,ytr,yte)

    vp = v(:);

    [~,ord] = sort(abs(vp),'ascend');

    curve = nan(numel(vp)+1,1);

    th = fitSigmoid1D(Xtr*vp,ytr);

    curve(1) = corrR2(yte,sigmoid1D(Xte*vp,th));

    best = curve(1);

    for k = 1:numel(ord)

        tr = vp;
        tr(ord(k)) = 0;

        th = fitSigmoid1D(Xtr*tr,ytr);

        sc = corrR2(yte,sigmoid1D(Xte*tr,th));

        curve(k+1) = sc;

        if isnan(sc) || (~isnan(best) && sc<best-0.01)
            curve = curve(1:k+1);
            break;
        end

        vp = tr;

        if isnan(best) || sc>best
            best = sc;
        end
    end

    if norm(vp)>0
        vp = vp/norm(vp);
    end
end

function [Vp,curve] = pruneTwoMID(V,Xtr,Xte,ytr,yte)

    Vp = V;
    d = size(V,1);

    mags = [abs(V(:,1));abs(V(:,2))];

    [~,ord] = sort(mags,'ascend');

    curve = nan(2*d+1,1);

    XY = Xtr*Vp;

    th = fitSigmoid2D(XY(:,1),XY(:,2),ytr);

    XYt = Xte*Vp;

    curve(1) = corrR2(yte,...
        sigmoid2D(XYt(:,1),XYt(:,2),th));

    best = curve(1);

    for k = 1:numel(ord)

        tr = Vp;
        q = ord(k);

        if q<=d
            tr(q,1) = 0;
        else
            tr(q-d,2) = 0;
        end

        if norm(tr(:,1))>eps
            tr(:,1) = tr(:,1)/norm(tr(:,1));
        end

        if norm(tr(:,2))>eps && norm(tr(:,1))>eps
            tr(:,2) = tr(:,2)-tr(:,1)*(tr(:,1)'*tr(:,2));
        end

        if norm(tr(:,2))>eps
            tr(:,2) = tr(:,2)/norm(tr(:,2));
        end

        XY = Xtr*tr;

        th = fitSigmoid2D(XY(:,1),XY(:,2),ytr);

        XYt = Xte*tr;

        sc = corrR2(yte,...
            sigmoid2D(XYt(:,1),XYt(:,2),th));

        curve(k+1) = sc;

        if isnan(sc) || (~isnan(best) && sc<best-0.01)
            curve = curve(1:k+1);
            break;
        end

        Vp = tr;

        if isnan(best) || sc>best
            best = sc;
        end
    end
end

function addUnityLine(a,b)

    lo = min([a(:);b(:)]);
    hi = max([a(:);b(:)]);

    if ~isfinite(lo) || ~isfinite(hi)
        return;
    end

    if lo==hi
        lo = lo-1;
        hi = hi+1;
    end

    plot([lo hi],[lo hi],'r--','LineWidth',1.2);
end