classdef stims
    %STIMS Class containing all functions for stimulus generation
    
    methods (Static, Access = public)
        function [f_grid, V_grid] = fvGridGen(sweepInfo)
            %FVGRIDGEN Creates frequency and voltage grids
            f_array = sweepInfo.f_min:sweepInfo.f_step:sweepInfo.f_max;
            V_array = sweepInfo.V_min:sweepInfo.V_step:sweepInfo.V_max;
            [f_grid, V_grid] = meshgrid(f_array, V_array);
        end

        function [f_grid, x_grid] = fxGridGen(sweepInfo)
            %FXGRIDGEN Creates frequency and position grids
            f_array = sweepInfo.f_min:sweepInfo.f_step:sweepInfo.f_max;
            x_array = sweepInfo.x_min:sweepInfo.x_step:sweepInfo.x_max;
            [f_grid, x_grid] = meshgrid(f_array, x_array);
        end

        function [f_grid, V_grid] = interGridGen(sweepInfo, calInfo, freq_ampl_resp)
            %INTERGRIDGEN Generates a V_grid by interpolating freq_ampl_resp
            [F_cal, V_cal] = gridGen(calInfo);

            C_cal = mean(freq_ampl_resp./V_cal, 1);
            f_0 = F_cal(1,:);
            
            f_array = sweepInfo.f_min:sweepInfo.f_step:sweepInfo.f_max;
            C_array = interp1(f_0, C_cal, f_array);
            
            [f_grid, x_grid] = stims.fxGridGen(sweepInfo);
            x_array = x_grid(:,1);
            V_grid = x_array * (1./C_array);
        end

        function [dur] = durCalc(sweepInfo)
            %DURCALC This function calculates the duration of the stimulus
            if isfield(sweepInfo, 'x_min')
                [f_grid, ~] = stims.fxGridGen(sweepInfo);
            elseif isfield(sweepInfo, 'V_min')
                [f_grid, ~] = stims.fvGridGen(sweepInfo);
            else
                warning('Check your sweepInfo and calInfo.')
            end
            
            t_grid = (sweepInfo.cycle)./f_grid;
            sample_grid = uint32(t_grid*sweepInfo.sampling);
            
            dur.sample = uint32(1/sweepInfo.fps*sweepInfo.sampling)...
                + uint32(2*sweepInfo.t_delay*sweepInfo.sampling)...
                + sum(sample_grid, "all")...
                + uint32(sweepInfo.t_rest*(numel(t_grid)-1)*sweepInfo.sampling);
            dur.time = double(dur.sample)/sweepInfo.sampling; % [s]
            dur.frame = uint32(dur.time*sweepInfo.fps);
        end

        function [x, t] = stimGen(sweepInfo, calInfo, freq_ampl_resp)
        %STIMGEN Generates x(t) for a given sweepInfo
        if nargin < 2
            [f_grid, V_grid] = stims.fvGridGen(sweepInfo);
        else
            [f_grid, V_grid] = stims.interGridGen(sweepInfo, calInfo, freq_ampl_resp);
        end

        t_grid = (sweepInfo.cycle)./f_grid;
        sample_grid = uint32(t_grid*sweepInfo.sampling);
        
        n_V = size(f_grid, 1);
        n_f = size(f_grid, 2);
        
        smp = 1/sweepInfo.sampling;
        
        % generating the initial delay
        x = zeros(1, uint32((sweepInfo.t_delay+1/sweepInfo.fps)*sweepInfo.sampling));
        
        % generating the resting period
        rest = zeros(1, uint32(sweepInfo.t_rest*sweepInfo.sampling));
        
        % generating all stimuli
        switch sweepInfo.type
            case 'sine'
                for i = 1:n_V
                    for k = 1:n_f
                        t_ik = linspace(0,t_grid(i,k)-smp,uint32(sample_grid(i,k)));
                        x_ik = V_grid(i,k)*sin(2*pi*f_grid(i,k)*t_ik);
                        x = [x, x_ik, rest];
                        clear t_ik x_ik
                    end
                end
            case 'square'
                for i = 1:n_V
                    for k = 1:n_f
                        t_ik = linspace(0,t_grid(i,k)-smp,uint32(sample_grid(i,k)));
                        x_ik = V_grid(i,k)*square(2*pi*f_grid(i,k)*t_ik);
                        x = [x, x_ik, rest];
                        clear t_ik x_ik
                    end
                end
            case 'triangle'
                for i = 1:n_V
                    for k = 1:n_f
                        t_ik = linspace(0,t_grid(i,k)-smp,uint32(sample_grid(i,k)));
                        x_ik = V_grid(i,k)*sawtooth((2*pi*f_grid(i,k)*t_ik + pi/2), 1/2);
                        x = [x, x_ik, rest];
                        clear t_ik x_ik
                    end
                end
            otherwise
                warning('Unexpected sweep type. No sweep created.')
        end
        
        % adding the final delay
        end_delay = sweepInfo.t_delay - sweepInfo.t_rest;
        if end_delay < 0
            disp('The delay time is shorter than the rest time.')
        end
        x = [x, zeros(1, uint32(end_delay*sweepInfo.sampling))];
        
        duration = stims.durCalc(sweepInfo);
        t = 0:smp:(duration.time-smp);
        end

        function camFileGen(filesLoc, sweepInfo)    % old function
        %CAMFILEGEN Generates and plots camera file

        duration = tools.durCalc(sweepInfo);
        smp = 1/sweepInfo.sampling;
        time = 0:smp:(duration.time-smp); % already includes one frame shift
        camera = (square(2*pi*sweepInfo.fps*time,50)+1)*5/2;

        disp('The camera file has been created')
        disp(['Total time: ', num2str(duration.time), ' s'])
        disp(['Frame count: ', num2str(duration.frame)])
        disp(['Camera signal: ', num2str(length(camera)), ' points'])

        figure('Name','Camera Signal','NumberTitle','off');
        plot(time, camera);
        axis([0 duration.time -1 6]);
        title('Camera Signal')
        xlabel('time, s')
        ylabel('voltage, V')

        strgfps = string(sweepInfo.fps);
        strgsmp = string(sweepInfo.sampling/1000);
        strgframes = string(duration.frame);
        writematrix(camera', sprintf('%s\\Camera_%sFps_Smp%skHz_%sFrames.txt', ...
            filesLoc, strgfps, strgsmp, strgframes))
        end

        function sweepInfo = stimFileGen(filesLoc, sweepInfo, calInfo, freq_ampl_resp)
        %STIMFILEGEN Generates and plots stimulus file
        %   Lo-hi freq and lo-hi amplitude sweep
        duration = stims.durCalc(sweepInfo);
        [stimulus, time] = stims.stimGen(sweepInfo, calInfo, freq_ampl_resp);

        sweepInfo.t_tot = duration.time;
        sweepInfo.frames = duration.frame;
        sweepInfo.points = length(stimulus);
        
        disp('The stimulus file has been created')
        disp(append("Total time: ", string(sweepInfo.t_tot), " s"))
        disp(append("Frame count: ", string(sweepInfo.frames)))
        disp(append("Calibration: ", string(sweepInfo.points), " points"))
        
        figure('Name','Stimulus','NumberTitle','off');
        plot(time, stimulus);
        title('Stimulus')
        xlabel('time, s')
        ylabel('voltage, V')
        
        strgfmin = string(sweepInfo.f_min);
        strgfmax = string(sweepInfo.f_max);
        strgfstep = string(sweepInfo.f_step);
        strgXmin = string(sweepInfo.x_min);
        strgXmax = string(sweepInfo.x_max);
        strgXstep = string(sweepInfo.x_step);
        writematrix(stimulus', sprintf('%s\\%s_%s-%s(%s)Hz_%s-%s(%s)nm.txt', ...
            filesLoc, sweepInfo.type, ...
            strgfmin, strgfmax, strgfstep, ...
            strgXmin, strgXmax, strgXstep));
        end

        function calInfo = calStimGen(filesLoc, calInfo)
        %CALSTIMGEN Generates stimulus file for calibration
        %   Lo-hi freq and lo-hi amplitude sweep
        
        duration = stims.durCalc(calInfo);
        [x, t] = stims.stimGen(calInfo);
        
        calInfo.t_tot = duration.time;
        calInfo.frames = duration.frame;
        calInfo.points = length(x);

        disp('The stimulus file has been created')
        disp(append("Total time: ", string(calInfo.t_tot), " s"))
        disp(append("Frame count: ", string(calInfo.frames)))
        disp(append("Calibration: ", string(calInfo.points), " points"))
        
        figure('Name','Stimulus for Calibration','NumberTitle','off');
        plot(t, x);
        title('Stimulus for Calibration')
        xlabel('time, s')
        ylabel('voltage, V')
        
        strgfmin = string(calInfo.f_min);
        strgVmin = string(calInfo.V_min);
        strgVmax = string(calInfo.V_max);
        strgVstep = string(calInfo.V_step);
        writematrix(x', sprintf('%s\\Calibration_%sHz_%s-%s(%s)V.txt', ...
            filesLoc, strgfmin, strgVmin, strgVmax, strgVstep));
        end
    end
end

