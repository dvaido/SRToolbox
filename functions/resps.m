classdef resps
    %RESPS Class containing all functions for response analysis

    methods (Static, Access = public)
        function [x, y, t, fps] = getBundleTrace(filesLoc, resultName)
            %GETBUNDLETRACE Exports data from a results file
            %   Outputs x-y position, time, and fps
            
            resultPath = sprintf("%s\\%s.mat", filesLoc, resultName);
            resultsParser = ResultsParser(resultPath);
            
            x = resultsParser.getProcessedTrace();
            y = resultsParser.getProcessedTrace2();
            t = resultsParser.getTime();
            fps = resultsParser.getFps();
        end

        function [i_start, i_end] = indexGen(sweepInfo)
            % Finds the starting and ending index for each stimulus
            % i_start and i_end are matrices containing the indices
            if isfield(sweepInfo, 'x_min')
                [f_grid, ~] = stims.fxGridGen(sweepInfo);
            elseif isfield(sweepInfo, 'V_min')
                [f_grid, ~] = stims.fvGridGen(sweepInfo);
            else
                warning('Check your sweepInfo and calInfo.')
            end
            t_grid = (sweepInfo.cycle)./f_grid;
            
            i_start = zeros(size(t_grid));
            i_end = zeros(size(t_grid));
            
            n_V = size(t_grid, 1);
            n_f = size(t_grid, 2);
            
            counter = uint32(sweepInfo.t_delay*sweepInfo.fps);
            for c_V = 1:n_V
                for c_f = 1:n_f
                    i_start(c_V, c_f) = counter;
                    counter = counter + uint32(t_grid(c_V, c_f)*sweepInfo.fps);
                    i_end(c_V, c_f) = counter;
                    counter = counter + uint32(sweepInfo.t_rest*sweepInfo.fps);
                end
            end
        end

        function smooth_and_flat = preprocess(x, drift_wind)
            %PREPROCESS Denoises and detrends the data stored in x
            no_lin_trend = detrend(x);
            flat = no_lin_trend - movmean(no_lin_trend, drift_wind, 1);
            smooth_and_flat = wdenoise(flat);
        end

        function [Ampl, freq] = getAmpl_PSD(x, y, fps)
            %GETPSD Calculates the amplitude from the PSD estimation
            %   Note that if x and y are both present the PSDs are added.
            nfft = length(x);
            
            if size(x, 2) > 1
                warning('getAmpl_PSD assumes that x is a column vector')
            end
    
            if isempty(y)
              [PSD, freq] = periodogram(x, [], nfft, fps);
            else
              [PSD, freq] = periodogram([x, y], [], nfft, fps);
              PSD = PSD(:, 1) + PSD(:, 2);
            end
            
              Ampl = sqrt(2*sum(PSD, 2)*fps/nfft);
        end

        function inst_ampl = instAmpl(x, amplitude)
            %INSTAMPL Calculates instantaneous amplitude using envelopes
            outliers_removed = medfilt1(x, amplitude.outlier, 1);
            [up, ~] = envelope(outliers_removed, amplitude.up, 'peak');
            [~, lo] = envelope(outliers_removed, amplitude.lo, 'peak');
            inst_ampl = up - lo;
            inst_ampl(inst_ampl < 0) = 0;
        end

        function [data_fit, decision_boundary, top_peak, top_width, q_frac] = burstAnalyzer(x, Hart)
            %BIMODALANALYZER Checks for bimodality (bursting) and calculates
            %   the amplitude of the oscilations (with st dev)
            data_fit = zeros(size(x));
            decision_boundary = zeros(1, size(x, 2));
            top_peak = zeros(1, size(x, 2));
            top_width = zeros(1, size(x, 2));
            q_frac = zeros(1, size(x, 2));
            for k = 1:size(x, 2)
                % Hartigan's dip test
                [p, dip, ~, ~] = dipTest(x(:,k));
                is_bimodal = (dip > Hart.dip) & (p < Hart.p);
        
                % Values of x for PDF and CDF
                x_vals = linspace(min(x(:,k)), max(x(:,k)), length(x(:,k)))';
                % options = statset('MaxIter', 500); % fit options
                if is_bimodal
                    % Fit a Gaussian Mixture Model with 2 components
                    gm = fitgmdist(x(:,k), 2);%, 'Options', options);
                    mu1 = gm.mu(1);                 % Mean of the first mode
                    mu2 = gm.mu(2);                 % Mean of the second mode
                    sigma1 = sqrt(gm.Sigma(1));     % Standard deviation of the first mode
                    sigma2 = sqrt(gm.Sigma(2));     % Standard deviation of the second mode
                    pi1 = gm.ComponentProportion(1);% Mixing proportion for the first mode
                    pi2 = gm.ComponentProportion(2);% Mixing proportion for the second mode
    
                    % Compute the decision boundary
                    f1 = @(x) pi1*(1/(sigma1*sqrt(2*pi)))*exp(-((x-mu1).^2)/(2*sigma1^2));
                    f2 = @(x) pi2*(1/(sigma2*sqrt(2*pi)))*exp(-((x-mu2).^2)/(2*sigma2^2));
                    x_search_range = linspace(mu1, mu2, 1000);
                    pdf1 = f1(x_search_range);
                    pdf2 = f2(x_search_range);
                    [~, idx] = min(abs(pdf1 - pdf2)); % Find closest point of equality
                    decision_boundary(k) = x_search_range(idx);
        
                    % Compute the position of the top peak
                    if mu1 > mu2
                        top_peak(k) = mu1;
                        top_width(k) = sigma1;
                    else
                        top_peak(k) = mu2;
                        top_width(k) = sigma2;
                    end
            
                    % Compute the two-state fraction = p1 / (p1 + p2)
                    q_frac(k) = sum(x(:,k) < decision_boundary(k))/length(x(:,k));
        
                    % Evaluate the PDF from the fit
                    data_fit(:,k) = pdf(gm, x_vals);
                else
                    gm = fitgmdist(x(:,k), 1, 'Options', options); % Gaussian model
                    data_fit(:,k) = pdf(gm, x_vals);
                    top_peak(k) = gm.mu(1);
                    top_width(k) = sqrt(gm.Sigma(1));
                    decision_boundary(k) = NaN;
                    q_frac(k) = NaN;
                end
            end
        end

        function inst_freq = instFreq(t, x, inst_ampl, decision_boundary, phase_smooth)
            % INSTFREQ Calculates instantaneous frequency
            dt = t(2) - t(1);
            analytic_signal = hilbert(x);
            phase = unwrap(angle(analytic_signal));  % Unwrap phase
            inst_freq_0 = movmean(diff(phase), phase_smooth)/(2*pi*dt);
            inst_freq = [NaN(1, size(inst_freq_0, 2)); inst_freq_0];
    
            % Set the instantaneous frequency to NaN if the amplitude is below the threshold
            inst_freq(inst_ampl < decision_boundary) = NaN;
        end

        function [freq_fit, mu_freq, sigma_freq] = freqAnalyzer(inst_freq)
            % Determine the average frequency with st dev
            freq_fit = zeros(size(inst_freq));
            mu_freq = zeros(1, size(inst_freq, 2));
            sigma_freq = zeros(1, size(inst_freq, 2));
            options = statset('MaxIter', 500); % fit options
            for k = 1:size(inst_freq, 2)
                gm = fitgmdist(inst_freq(:,k), 1, 'Options', options); % Gaussian model
                mu_freq(k) = gm.mu(1); % Mean of the first mode
                sigma_freq(k) = sqrt(gm.Sigma(1)); % Standard deviation of the first mode
                f_vals = linspace(min(inst_freq(:,k)), max(inst_freq(:,k)), length(inst_freq(:,k)))';
                freq_fit(:, k) = pdf(gm, f_vals);
            end
        end

        function freq_ampl_resp = extractCal(x, y, params, calInfo)
            %EXTRACTCAL 
            if size(x, 2) > size (x, 1)
                x = x';
                y = y';
            end
    
            if ~isempty(y)
                x = [x, y];
                y_is_present = 1;
            else
                y_is_present = 0;
            end
            
            % Denoise and detrend
            x_smooth = resps.preprocess(x, params.drift_wind);
    
            % Get the indices and frequency grid
            [i_start, i_end] = resps.indexGen(calInfo);
            [f_grid, ~] = stims.fvGridGen(calInfo);
    
            n_V = size(f_grid, 1);
            n_f = size(f_grid, 2);
    
            freq_ampl_resp = zeros(n_V, n_f);
    
            for i = 1:n_V
                for k = 1:n_f
                    x_ik = x_smooth(i_start(i, k):i_end(i, k), 1);
                    switch y_is_present
                        case 1
                            y_ik = x_smooth(i_start(i, k):i_end(i, k), 2);
                        otherwise
                            y_ik = [];
                    end
    
                    [x_f, freq] = resps.getAmpl_PSD(x_ik, y_ik, calInfo.fps);
                    freq_ampl_resp(i, k) = interp1(freq, x_f, f_grid(i, k), 'pchip');
                end
            end
        end

        function [ampl_mu, ampl_sigma, ampl_FFT] = analyzeSweep(x, params, sweepInfo)
            %ANALYZESWEEP
            if size(x, 2) > size (x, 1)
                x = x';
            end
            
            % Denoise and detrend
            x_smooth = resps.preprocess(x, params.drift_wind);
    
            % Instantaneous amplitude
            x_inst_ampl = resps.instAmpl(x_smooth, params.amplitude);
    
            % Get the indices and frequency grid
            [i_start, i_end] = resps.indexGen(sweepInfo);
            [f_grid, ~] = stims.gridGen(sweepInfo);
    
            n_V = size(f_grid, 1);
            n_f = size(f_grid, 2);
    
            ampl_mu = zeros(n_V, n_f);
            ampl_sigma = zeros(n_V, n_f);
            ampl_FFT = zeros(n_V, n_f);
    
            for i = 1:n_V
                for k = 1:n_f
                    ampl_ik = x_inst_ampl(i_start(i, k):i_end(i, k), :);
                    ampl_mu(i, k) = sum(mean(ampl_ik, 1));
                    ampl_sigma(i, k) = norm(std(ampl_ik, 1));
    
                    x_ik = x_smooth(i_start(i, k):i_end(i, k), 1);
    
                    [x_f, freq] = resps.getAmpl_PSD(x_ik, [], sweepInfo.fps);
                    ampl_FFT(i, k) = interp1(freq, x_f, f_grid(i, k), 'pchip');
                end
            end
        end
    end
end

