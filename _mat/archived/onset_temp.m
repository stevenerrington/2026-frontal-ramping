%% ============================================================
% Piecewise (hinge) regression to estimate ramp onset
%
% Model:
%
%   y = b0 + b1*t + b2*max(0,t-bp)
%
% Ramp onset = bp
%
%% ============================================================

clearvars -except fixon_array ramping_neuron_idx

rng(1)

%% Parameters

cond = 3;                  % Medium foreperiod

time = -2000:2000;         % ms

fitWindow = time >= -500 & time <= 800;
searchWindow = time >= 0 & time <= 1000;

smoothSigma = 100;

nExamples = 10;

%% Random neurons

exampleNeurons = ramping_neuron_idx(randperm(length(ramping_neuron_idx),nExamples));

figure
set(gcf,'Position',[100 100 1500 700])

for ii = 1:nExamples

    neuron = exampleNeurons(ii);

    %% Get PSTH

    raw = squeeze(fixon_array(neuron,:,cond));

    raw = raw(:);

    smooth = smoothdata(raw,'gaussian',smoothSigma);

    %% Restrict fit

    t = time(fitWindow)';
    y = smooth(fitWindow);

    %% Initialise

    bestSSE = inf;
    bestFit = [];
    bestBP = NaN;
    bestB = [];

    %% Candidate breakpoints

    candidateTimes = time(searchWindow);

    %% Search

    for bp = candidateTimes

        hinge = max(0,t-bp);

        X = [ones(length(t),1), ...
             t, ...
             hinge];

        b = X\y;

        yhat = X*b;

        SSE = sum((y-yhat).^2);

        if SSE < bestSSE

            bestSSE = SSE;
            bestFit = yhat;
            bestBP = bp;
            bestB = b;

        end

    end

    %% Slopes

    baselineSlope = bestB(2);
    rampSlope = bestB(2)+bestB(3);

    %% Plot

    subplot(2,5,ii)
    hold on

    plot(time,raw,...
        'Color',[0.8 0.8 0.8],...
        'LineWidth',1)

    plot(time,smooth,...
        'k',...
        'LineWidth',2)

    plot(t,bestFit,...
        'b',...
        'LineWidth',2)

    xline(bestBP,'r--','LineWidth',2)

    xlabel('Time (ms)')
    ylabel('FR')

    title(sprintf('Neuron %d\nOnset=%d ms', ...
        neuron,bestBP))

    xlim([-500 800])

    box off

end

sgtitle('Piecewise linear fit')