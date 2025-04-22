clc
clear
% close all

addpath(fullfile(pwd, 'functions'));

%%

% filesLoc = 'path\to\results';
filesLoc = 'C:\Users\Dima\Documents\MATLAB\SRT\test\next';

calInfo.type = 'sine';
calInfo.f_min = 5;
calInfo.f_step = 5;
calInfo.f_max = 100;
calInfo.V_min = 0.05;
calInfo.V_step = 0.05;
calInfo.V_max = 0.05;
calInfo.cycle = 10;
calInfo.sampling = 20000;
calInfo.fps = 2000;
calInfo.t_delay = 0.1;
calInfo.t_rest = 0.05;

calInfo = stims.calStimGen(filesLoc, calInfo);
writestruct(calInfo, [filesLoc, '\calInfo.json'])

%%

% filesLoc = 'path\to\results';
filesLoc = 'C:\Users\Dima\Documents\MATLAB\SRT\test';
resultName = 'JetCal';

[x, y, t, ~] = resps.getBundleTrace(filesLoc, resultName);
calInfo = readstruct([filesLoc, '\calInfo.json']);

params.drift_wind = 2000;       % window for detranding (movmean)
params.amplitude.outlier = 8;   % median filter window for amplitude outliers
params.amplitude.lo = 400;      % window for low envelope
params.amplitude.up = 10;       % window for high envelope
params.Hart.dip = 0.0001;       % dip threshold for Hartigan's test
params.Hart.p = 0.002;          % p-value threshold for Hartigan's test
params.phase_smooth = 100;      % window for phase smoothing (movmean)

freq_ampl_resp = resps.extractCal(x, y, params, calInfo);

sweepInfo.type = 'triangle';
sweepInfo.f_min = 50;
sweepInfo.f_step = 5;
sweepInfo.f_max = 50;
sweepInfo.x_min = 50;
sweepInfo.x_step = 10;
sweepInfo.x_max = 50;
sweepInfo.cycle = 10;
sweepInfo.sampling = 20000;
sweepInfo.fps = 2000;
sweepInfo.t_delay = 0.1;
sweepInfo.t_rest = 0.05;

sweepInfo = stims.stimFileGen(filesLoc, sweepInfo, calInfo, freq_ampl_resp);
writestruct(sweepInfo, [filesLoc, '\sweepInfo.json'])
