function ppg_diary_pilot1_analysis
% ppg_diary_pilot1_analysis analyses data from the PPG Diary Pilot Study.
%
%               ppg_diary_pilot1_analysis
%
%	Inputs:
%       - A Matlab file containing the collated data, called "ppg_diary_pilot_conv_data.mat"
%               (this can be downloaded from https://doi.org/10.5281/zenodo.3268500 )
%       - The path to this file specified in the "universal_parameters"
%       function below. 
%
%	Outputs:
%       - Data processing files (in the same location as the Matlab file above)
%       - Results output verbatim in the command window
%       - Plots
%           
%   Further Information:
%       This version of ppg_diary_pilot1_analysis is provided to facilitate
%       reproduction of the PPG Diary Pilot analysis reported in:
%           Charlton P. H. et al., "Acquiring wearable photoplethysmography
%           data in daily life: the PPG Diary Pilot Study", [under review] 
%       Further information on this study can be obtained at:
%           https://peterhcharlton.github.io/ppg-diary/
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: https://peterhcharlton.github.io/ppg-diary/contributions.html
%
%   Version:
%       v.0.1 - accompanying peer review, 16th Oct 2020 by Peter H Charlton
%
%   Source:
%       This script has been adapted from scripts in the RRest toolbox,
%       available at:
%           https://github.com/peterhcharlton/RRest
%
%   Licence:
%       This program is available under the GNU public license. This
%       program is free software: you can redistribute it and/or modify 
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version. This program is distributed in
%       the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%       even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%       PARTICULAR PURPOSE.  See the GNU General Public License for more
%       details: <http://www.gnu.org/licenses/>.

%% Setup universal parameters
% ~~~ This function needs adjusting by the individual user as it contains paths
% specific to their computer ~~~
up = universal_parameters;

%% Check Matlab data file is available
if ~exist(up.paths.conv_data_file, 'file')
     fprintf('\n\n\n ~~~~~~~~~~~~\n Couldn''t find Matlab file at:\n   %s\n\n Please make sure you''ve downloaded it from:\n   https://doi.org/10.5281/zenodo.3268500\n ~~~~~~~~~~~~\n', up.paths.conv_data_file)
     return
end

%% Filter data
do_step = 0;
if ~exist(up.paths.filt_data_file, 'file') || do_step
    filter_data(up);
end

%% Identify beats in filtered data
do_step = 0;
if ~exist(up.paths.beat_data_file, 'file') || do_step
    identify_beats_in_data(up);
end

%% Export sample data
do_step = 0;
if ~exist(up.paths.sample_data_file_1_hour, 'file') || do_step
    export_sample_data(up);
end

%% Identify recording periods
do_step = 0;
if ~exist(up.paths.rec_periods_file, 'file') || do_step
    identify_rec_periods(up);
end

%% Perform signal quality analysis
do_step = 0;
if ~exist(up.paths.quality_data_file, 'file') || do_step
    assess_signal_quality(up);
end

%% Calculate trial outcomes
do_step = 1;
if ~exist(up.paths.trial_outcomes_file, 'file') || do_step
    calculate_trial_outcomes(up)
end

%% Make plots
do_step = 1;
if do_step
    make_plots(up);
end

end

function up = universal_parameters

fprintf('\n --- Setting up Universal Parameters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the root data directory (where the Matlab file is, and where the
% results will be stored)
up.paths.root_folder = '/Users/petercharlton/Desktop/temp/ppgdiary/'; % note the need for a slash at the end of the path
fprintf('\n   - The results will be stored at:\n   %s', up.paths.root_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   (don't need editing)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

up.paths.conv_data_file = [up.paths.root_folder, 'ppg_diary_pilot_conv_data.mat'];
up.paths.filt_data_file = [up.paths.root_folder, 'filt_data.mat'];
up.paths.beat_data_file = [up.paths.root_folder, 'beat_data.mat'];
up.paths.rec_periods_file = [up.paths.root_folder, 'rec_periods.mat'];
up.paths.trial_outcomes_file = [up.paths.root_folder, 'trial_outcomes.mat'];
up.paths.sample_data_file_1_hour = [up.paths.root_folder, 'PPGdiary1_1_hour_sample.mat'];
up.paths.sample_data_file_1_min = [up.paths.root_folder, 'PPGdiary1_1_min_sample.mat'];
up.paths.sample_data_file_1_pulse = [up.paths.root_folder, 'PPGdiary1_1_pulse_sample.mat'];
up.paths.quality_data_file = [up.paths.root_folder, 'quality_data.mat'];
up.paths.plot_folder = [up.paths.root_folder, 'plots', filesep];
if ~exist(up.paths.plot_folder), mkdir(up.paths.plot_folder); end

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Fpass = 0.6;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.3;   % in Hz     
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.paramSet.elim_vhf.Fpass = 20;  % in HZ
up.paramSet.elim_vhf.Fstop = 15;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.01;

% Beat detection settings
up.beat_detect.win_length = 10; % in secs
up.beat_detect.win_step_prop = 0.8;

end

function export_sample_data(up)

fprintf('\n --- Exporting sample data')

% load filtered data
load(up.paths.filt_data_file)

% identify one hour of data shortly after an annotation of sleeping
t = (1552865724-1550497798)*1000; % time in ms of sleeping since start
rel_t = [t+0.5*60*1000, t+60.5*60*1000];
rel_els = data.t_ms_rel>=rel_t(1) & data.t_ms_rel<=rel_t(2);

% extract this as a sample dataset
sample.ppg_ir = -1*data.ppg_ir(rel_els);
sample.ppg_red = -1*data.ppg_red(rel_els);
sample.fs = data.fs;
sample.description = 'Infra-red and red PPG recording from PPG Diary Pilot Study. A one-hour recording acquired during sleep.';
clear data

% save data
data = sample;
save(up.paths.sample_data_file_1_hour, 'data')

% save data (1 min sample)
no_els = 1+(sample.fs*60);
start_el = 16*60*sample.fs;
rel_els = start_el:start_el+no_els-1;
S.v = sample.ppg_ir(rel_els);
S.fs = sample.fs;
S.description = 'Infra-red PPG recording from PPG Diary Pilot Study. A one-minute recording acquired during sleep.';
save(up.paths.sample_data_file_1_min, 'S')

% save data (1 pulse wave sample)
rel_els = start_el+172:start_el+339;
S.v = sample.ppg_ir(rel_els);
S.fs = sample.fs;
S.description = 'Infra-red PPG recording from PPG Diary Pilot Study. A recording of a single pulse wave acquired during sleep.';
save(up.paths.sample_data_file_1_pulse, 'S')

end

function identify_rec_periods(up)

fprintf('\n --- Identifying recording periods')

%% Load ADL data
load(up.paths.conv_data_file, 'adl_data')

%% Identifying recording periods
% In each recording period the activity and connection status stayed
% the same. Otherwise, start a new recording period.

% get rid of actions which are actually notes
note_els = ~cellfun(@isempty, strfind(adl_data.act, 'note_'));
adl_data.act = adl_data.act(~note_els);
adl_data.t_ms_rel = adl_data.t_ms_rel(~note_els);
clear note_els

% get rid of spaces at end of actions
for s = 1 : length(adl_data.act)
    if strcmp(adl_data.act{s}(end), '_')
        adl_data.act{s} = adl_data.act{s}(1:end-1);
    end    
end

% identify recording periods
[rec_periods.deb, rec_periods.fin] = deal(nan(length(adl_data.act)-1,1));
[rec_periods.act, rec_periods.reason_disconnect] = deal(cell(length(adl_data.act)-1,1));
rec_periods.connected = false(length(adl_data.act)-1,1);
curr_connection_status = 0;
curr_activity = 'unknown';
rec_period_counter = 0;
for act_no = 1 : length(adl_data.act)
    this_txt = adl_data.act{act_no};
    this_time = adl_data.t_ms_rel(act_no);
    
    % - end of study
    if length(this_txt)>=29 && strcmp(this_txt(1:29), 'disconnected_-_study_finished')
        
            
        % - insert finish time of this period
        if rec_period_counter>0
            rec_periods.fin(rec_period_counter) = this_time;
        end
        
        break
    
    % - connections
    elseif strcmp(this_txt, 'connected')
        
        curr_disconnection_reason = '';
        
        if curr_connection_status == 0
            
            % - insert finish time of this period
            if rec_period_counter>0
                rec_periods.fin(rec_period_counter) = this_time;
            end
            
            % - move to next period
            rec_period_counter = rec_period_counter+1;
            curr_connection_status = 1;
            
            % - insert start time of this period
            rec_periods.deb(rec_period_counter) = this_time;
        else
            
            fprintf('Currently ignoring this - a connection without a disconnection beforehand')
            
        end
        
    % - disconnections
    elseif length(this_txt)>=12 && strcmp(this_txt(1:12), 'disconnected')
        
        curr_disconnection_reason = this_txt(16:end);
        
        if curr_connection_status == 1
            
            % - insert finish time of this period
            if rec_period_counter>0
                rec_periods.fin(rec_period_counter) = this_time;
            end
            
            % - move to next period
            rec_period_counter = rec_period_counter+1;
            curr_connection_status = 0;
            
            % - insert start time of this period
            rec_periods.deb(rec_period_counter) = this_time;
        else
            
            % - insert finish time of this period
            if rec_period_counter>0
                rec_periods.fin(rec_period_counter) = this_time;
            end
            
            % - move to next period
            rec_period_counter = rec_period_counter+1;
            curr_connection_status = 0;
            
            % - insert start time of this period
            rec_periods.deb(rec_period_counter) = this_time;
            
        end
    % end of activity
    elseif contains(this_txt, '_fin')
        
        rec_periods.fin(rec_period_counter) = this_time;
        
        % - move to next period
        rec_period_counter = rec_period_counter+1;
        
        % - insert start time of this period
        rec_periods.deb(rec_period_counter) = this_time;
        
        curr_activity = 'unknown';
        
    % start of activity    
    else
        
        rec_periods.fin(rec_period_counter) = this_time;
        
        % - move to next period
        rec_period_counter = rec_period_counter+1;
        
        % - insert start time of this period
        rec_periods.deb(rec_period_counter) = this_time;
                
        curr_activity = this_txt;
        
        
    end
    
    % Fill in fields
    rec_periods.reason_disconnect{rec_period_counter} = curr_disconnection_reason;
    rec_periods.connected(rec_period_counter) = curr_connection_status;
    rec_periods.act{rec_period_counter} = curr_activity;

end

% - calculate durations (in secs)
rec_periods.durn = (rec_periods.fin - rec_periods.deb)/1000;

% identify periods in which the participant noted they were not wearing the sensor 
rec_periods.noted_wearing = false(length(rec_periods.deb),1);
for s = 1 : length(rec_periods.noted_wearing)
    curr_reason_disconnect = rec_periods.reason_disconnect{s};
    % the participant noted they were wearing the sensor if it was connected
    if isempty(curr_reason_disconnect)
        rec_periods.noted_wearing(s) = true;
    end
    % the participant noted they were wearing the sensor if it was disconnected but there was no reason recorded for this
    if strcmp(curr_reason_disconnect, 'unknown')
        rec_periods.noted_wearing(s) = true;
    end
    % the participant noted they were wearing the sensor if it was disconnected but this was not intentional
    if strcmp(curr_reason_disconnect, 'connection_dropped')
        rec_periods.noted_wearing(s) = true;
    end
end

% save rec periods
save(up.paths.rec_periods_file, 'rec_periods')

end

function make_plots(up)

fprintf('\n --- Making plots')

%% Load data
fprintf('\n  - Loading data')
load(up.paths.beat_data_file);
load(up.paths.trial_outcomes_file);
load(up.paths.rec_periods_file);

%% Prepare data

% Extract heart rates
fprintf('\n  - Extracting heart rates')
hr.t = data.t_ms_rel(data.ppg_ir_p)/1000; % in secs
hr.v = round(60./diff(hr.t));
hr.t = hr.t(1:end-1);
rel_els = hr.v > 30 & hr.v < 180;
hr.t = hr.t(rel_els);
hr.v = hr.v(rel_els);

% obtain troughs
fprintf('\n  - Identifying troughs')
data.ppg_ir_tr = nan(length(data.ppg_ir_p)-2,1);
for s = 1 : length(data.ppg_ir_p)-1
    [~, temp] = min(data.ppg_ir(data.ppg_ir_p(s):data.ppg_ir_p(s+1)));
    data.ppg_ir_tr(s) = temp + data.ppg_ir_p(s) - 1;
end

%% Plot proportion of time for which signal was high quality
% obtain data
perc_tst = endpoints.prop_time_hq.acts.prop;
total_durn = endpoints.prop_time_hq.acts.t;
total_n = endpoints.prop_time_hq.acts.n;
activities = endpoints.prop_time_hq.acts.act;
activities = strrep(activities, '_', ' ');
[perc_tst, order] = sort(perc_tst, 'descend');
total_durn = total_durn(order);
total_n = total_n(order);
activities = activities(order);
include_overall = 0;
if include_overall
    perc_tst = 100*[endpoints.prop_time_hq.overall.prop; perc_tst];
    total_durn = [endpoints.prop_time_hq.overall.t; total_durn]./(60*60);
    total_n = [1; total_n];
    activities = ['Overall'; activities];
else
    perc_tst = 100*perc_tst;
    total_durn = total_durn./(60*60);
end
activities = strrep(activities, 'tv', 'watching TV');
activities = strrep(activities, 'food preparation', 'preparing food');
activities = strrep(activities, 'eat', 'eating');
activities = strrep(activities, 'computer', 'using computer');
for act_no = 1 : length(activities)
    activities{act_no} = [upper(activities{act_no}(1)), activities{act_no}(2:end)];
end

% plot bars
fprintf('\n  - Plotting proportion of time for which signal was high quality')
ftsize = 20;
figure('Position', [20,20,800,450])
subplot('Position', [0.12,0.30,0.86,0.68])
X = categorical(activities);
X = reordercats(X,activities);
hbar = bar(X,perc_tst);
ylabel({'Percentage of time for which', 'PPG signal was high quality'}, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize-2, 'YGrid', 'on')
box off
xtickangle(33)

% plot overall dashed line
hold on
xlims = xlim;
overall_perc = 100*endpoints.prop_time_hq.overall.prop;
plot(xlims, [overall_perc, overall_perc], '--r', 'LineWidth', 2);
dim = [0.7 0.62 0.1 0.1];
str = sprintf('Overall percentage: %.1f%%', overall_perc);
annotation('textbox',dim,'String',str,'FontSize', ftsize-2, 'FitBoxToText','on', 'LineStyle', 'none');

do_time_labels = 0;

if do_time_labels
    % Add duration annotations: adapted from
    % https://uk.mathworks.com/matlabcentral/answers/128396-how-do-i-label-the-bars-in-my-bar-graph-in-matlab
    
    % Get the data for all the bars that were plotted
    x = get(hbar,'XData');
    y = get(hbar,'YData');
    
    ygap = 3;  % Specify vertical gap between the bar and label
    ylimits = get(gca,'YLim');
    set(gca,'YLim',[ylimits(1),ylimits(2)+0.2*max(y)]); % Increase y limit for labels
    
    % Create labels to place over bars
    total_durn = 0.1*round(total_durn*10); % total_n; %
    labels = cellstr(num2str(total_durn));
    labels = strrep(labels, ' ', '');
    
    for i = 1:length(x) % Loop over each bar
        xpos = x(i);        % Set x position for the text label
        ypos = y(i) + ygap; % Set y position, including gap
        htext = text(xpos,ypos,labels{i},'FontSize', ftsize -4);          % Add text label
        set(htext,'VerticalAlignment','bottom',...  % Adjust properties
            'HorizontalAlignment','center')
    end
    
    
    % % title
    % dim = [.32 .86 .1 .1];
    % str = 'Proportion of time for which PPG signal was of high quality';
    % annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize+4);
end

% save
filepath = [up.paths.plot_folder, 'ppg_diary_prop_hq_acts'];
print(filepath, '-depsc')
close all

% output proportion of time for which signal was of high quality in
% selected activities:
for s = 1 : length(activities)
    fprintf('\n %s: Signal of high quality for %.1f%% of the time', activities{s}, perc_tst(s))
end

%% Plot reasons for data loss
% obtain data
perc_tpwt = endpoints.reasons_for_not_wearing.perc_tpwt;
reasons = endpoints.reasons_for_not_wearing.reason;
reasons = strrep(reasons, '_', ' ');
[perc_tpwt, order] = sort(perc_tpwt, 'descend');
reasons = reasons(order);
bad_els = strcmp(reasons, 'day off');
perc_tpwt = perc_tpwt(~bad_els);
reasons = reasons(~bad_els);
for reason_no = 1 : length(reasons)
    reasons{reason_no} = [upper(reasons{reason_no}(1)), reasons{reason_no}(2:end)];
end
reasons(strcmp(reasons, 'Social')) = {'Socialising'};

% plot
fprintf('\n  - Plotting reasons for data loss')
ftsize = 20;
figure('Position', [20,20,800,450])
subplot('Position', [0.1,0.27,0.89,0.71])
X = categorical(reasons);
X = reordercats(X,reasons);
bar(X,perc_tpwt)
xtickangle(33)
ylabel({'Percentage of', 'total possible wear time'}, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize, 'YGrid', 'on')
box off

% % title
% dim = [.42 .86 .1 .1];
% str = 'Reasons for data loss';
% annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize+4);

% save
filepath = [up.paths.plot_folder, 'ppg_diary_data_loss_reasons'];
print(filepath, '-depsc')
close all
for s = 1 : 2
    fprintf('\n Removed for %s: %.1f%% of the total possible wear time', reasons{s}, perc_tpwt(s));
end

%% Plot HRs across whole dataset
fprintf('\n  - Plotting heart rates across whole dataset')
ftsize = 16;
figure('Position', [20,20,1100,300])
subplot('Position', [0.06,0.15,0.93,0.82])
plot(hr.t/(60*60*24), medfilt1(hr.v,30), '.b');
xlim([0 28])
ylim([30 180])
set(gca, 'FontSize', ftsize)
xlabel('Time (days)', 'FontSize', ftsize)
ylabel('Heart Rate (bpm)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:28;
ax.XAxis.TickValues = 0:7:28;
grid on

% annotations
% - rest day 1
annotation('ellipse',[.23 .23 .015 .05], 'Color', 'r', 'FaceColor', 'r')
ht = text(5.3,70,'Rest day', 'Color', 'r');
set(ht,'Rotation',90,'FontSize',18)
% - rest day 2
annotation('ellipse',[.465 .23 .015 .05], 'Color', 'r', 'FaceColor', 'r')
ht = text(12.4,70,'Rest day', 'Color', 'r');
set(ht,'Rotation',90,'FontSize',18)
% - rest day 3
annotation('ellipse',[.712 .23 .015 .05], 'Color', 'r', 'FaceColor', 'r')
ht = text(19.8,70,'Rest day', 'Color', 'r');
set(ht,'Rotation',90,'FontSize',18)
% - rest day 4
annotation('ellipse',[.875 .23 .015 .05], 'Color', 'r', 'FaceColor', 'r')
ht = text(24.75,70,'Rest day', 'Color', 'r');
set(ht,'Rotation',90,'FontSize',18)

% save
filepath = [up.paths.plot_folder, 'ppg_diary_all_hr'];
print(filepath, '-depsc')
close all

%% Plot IBIs across whole dataset
fprintf('\n  - Plotting inter-beat-intervals across whole dataset')
ftsize = 16;
figure('Position', [20,20,1000,600])
subplot('Position', [0.06,0.57,0.93,0.41])
plot(hr.t/(60*60*24), medfilt1(hr.v,30), '.b');
xlim([0 28])
ylim([30 180])
set(gca, 'FontSize', ftsize)
ylabel('Heart Rate (bpm)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:28;
ax.XAxis.TickValues = 0:7:28;
grid on
subplot('Position', [0.06,0.09,0.93,0.41])
ibi.v = 1000*60./hr.v;
plot(hr.t/(60*60*24), medfilt1(ibi.v,30), '.k');
xlim([0 28])
ylim(1000*60./[180 30])
set(gca, 'FontSize', ftsize)
xlabel('Time (days)', 'FontSize', ftsize)
ylabel('Inter-beat-intervals (ms)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:28;
ax.XAxis.TickValues = 0:7:28;
grid on
filepath = [up.paths.plot_folder, 'ppg_diary_all_ibis'];
print(filepath, '-depsc')
close all


%% Plot PPG wave
fprintf('\n  - Plotting sample of high quality PPG wave')
figure('Position', [20,20,700,300])
subplot('Position', [0.10,0.22,0.87,0.75])
t = (1552865724-1550497798)*1000+2*60*60*1000-1.81*1000;
t_end = t+10000;
rel_els = find(data.t_ms_rel>t & data.t_ms_rel<t_end);
plot((data.t_ms_rel(rel_els)-t)/1000, data.ppg_ir(rel_els), 'b', 'LineWidth', 2), hold on,
plot((data.t_ms_rel(data.ppg_ir_tr)-t)/1000, data.ppg_ir(data.ppg_ir_tr), 'or', 'LineWidth', 2)
xlim(([t, t_end]-t)/1000)
ylim([-6000 4000])
set(gca, 'FontSize', ftsize, 'YTick', [])
xlabel('Time (s)', 'FontSize', ftsize)
ylabel('PPG Signal (au)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:10;
ax.XAxis.TickValues = 0:2:10;
grid on

% annotations
ah=annotation('doublearrow',[.703 .827],.37*[1,1],'Color','r', 'LineWidth', 2);
dim = [.64 .25 .2 .1];
str = 'inter-beat-interval';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize, 'Color', 'r');

% save
filepath = [up.paths.plot_folder, 'ppg_diary_ppg'];
print(filepath, '-depsc')
close all


%% Plot pulse waves during running

fig_settings.action = 'running';
fig_settings.ftsize = 22;
fig_settings.action_start = (((1552727004-1550497798)/(60*60*24))+0.002)*24*60*60*1000+(14.3*60*1000); % in ms
fig_settings.action_fin = fig_settings.action_start+(1000*60*5); % in ms
fig_settings.hr_start = 0*60; % how long to display HR before action (in secs)
fig_settings.hr_durn = 15*60; % how long to display HR after action (in secs)
fig_settings.pw_start = (fig_settings.action_fin/1000)+2.63*60; % in secs  % 2.15, 2.7 ok, 2.63 perfect
fig_settings.hr_med_filt_order = 25;
fig_settings.hr_lims = [60, 170];
fig_settings.pw_gap = 2*60; % gap between PWs (in secs)
fig_settings.action_label_hr_height = 80;

make_event_plot(hr, data, fig_settings)

% save
filepath = [up.paths.plot_folder, 'ppg_diary_running_pws'];
print(filepath, '-depsc')
close all

%% Plot heart rates during running
fprintf('\n  - Plotting sample of running')
ftsize = 24;
figure('Position', [20,20,700,300])
subplot('Position', [0.15,0.22,0.83,0.75])
t = (1552727004-1550497798)/(60*60*24);
t_start = t+0.002;
plot(((hr.t/(60*60*24))-t_start)*(24*60)-7.751, medfilt1(hr.v,30), '.b'); xlim([0 0.025]*(24*60)-7.751)
set(gca, 'FontSize', ftsize)
xlabel('Time (mins)', 'FontSize', ftsize)
ylabel('Heart Rate (bpm)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -10:30;
ax.XAxis.TickValues = -5:5:25;
grid on

% add annotations
dim = [.2 .55 .1 .1];
str = 'Rest';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
dim = [.41 .85 .1 .1];
str = 'Running';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
dim = [.62 .68 .1 .1];
str = 'Recovery';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);

% save
filepath = [up.paths.plot_folder, 'ppg_diary_running'];
print(filepath, '-depsc')
close all

%% Plot pulse waves whilst falling asleep

fig_settings.action = 'falling asleep';
fig_settings.ftsize = 24;
rel_rec_period = find(strcmp(rec_periods.act, 'sleeping'),1);
fig_settings.action_start = rec_periods.deb(rel_rec_period); % in ms
fig_settings.action_fin = fig_settings.action_start+(1000*60*4); % in ms
fig_settings.hr_start = 1.7*60; % how long to display HR before action (in secs)
fig_settings.hr_durn = 26*60; % how long to display HR after action (in secs)
fig_settings.pw_start = (fig_settings.action_start/1000)+4*60+26; % in secs
fig_settings.hr_med_filt_order = 10;
fig_settings.hr_lims = [30, 100];
fig_settings.pw_gap = 4*60; % gap between PWs (in secs)
fig_settings.action_label_hr_height = 80;

make_event_plot(hr, data, fig_settings)

% save
filepath = [up.paths.plot_folder, 'ppg_diary_falling_asleep_pws'];
print(filepath, '-depsc')
close all

% %% Plot stair climbing
% 
% fig_settings.action = 'stair climbing';
% fig_settings.ftsize = 24;
% rel_rec_period = find(strcmp(rec_periods.act, 'climbing_stairs'),1);
% fig_settings.action_start = rec_periods.deb(rel_rec_period); % in ms
% fig_settings.action_fin = rec_periods.fin(rel_rec_period); % in ms
% fig_settings.hr_start = 0; % how long to display HR before action (in secs)
% fig_settings.hr_durn = 2*60; % how long to display HR after action (in secs)
% fig_settings.pw_start = 10; % since start of action (in secs)
% fig_settings.hr_med_filt_order = 10;
% fig_settings.hr_lims = [40, 160];
% fig_settings.pw_gap = 30; % gap between PWs (in secs)
% fig_settings.action_label_hr_height = 50;
% 
% make_event_plot(hr, data, fig_settings)
% 
% % save
% filepath = [up.paths.plot_folder, 'ppg_diary_stairs'];
% print(filepath, '-depsc')
% close all

%% Plot heart rate whilst sleeping
fprintf('\n  - Plotting sample of sleeping')
figure('Position', [20,20,700,300])
subplot('Position', [0.15,0.22,0.80,0.75])
t = (1552865724-1550497798)/(60*60*24);
plot(((hr.t/(60*60*24))-t)*(24*60), medfilt1(hr.v,30), '.b'); xlim([-2 2]*60)
set(gca, 'FontSize', ftsize)
xlabel('Time (mins)', 'FontSize', ftsize)
ylabel('Heart Rate (bpm)', 'FontSize', ftsize)
box off
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -120:10:120;
ax.XAxis.TickValues = -120:60:120;
grid on

% add annotations
dim = [.37 .8 .1 .1];
str = 'Awake';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
dim = [.55 .8 .1 .1];
str = 'Asleep';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);

% save
filepath = [up.paths.plot_folder, 'ppg_diary_sleeping'];
print(filepath, '-depsc')
close all

end

function make_event_plot(hr, data, fig_settings)

fprintf('\n  - Plotting sample of %s', fig_settings.action)
figure('Position', [20,20,1200,440])
subplot('Position', [0.10,0.16,0.85,0.52])
t_start = fig_settings.action_start/(1000); % in secs
t_fin = fig_settings.action_fin/(1000); % in secs
plot((hr.t-t_start)/60, medfilt1(hr.v,fig_settings.hr_med_filt_order), '.-b'); % in mins
xlim([(-fig_settings.hr_start/60) (t_fin-t_start+fig_settings.hr_durn)/60]) % in mins
set(gca, 'FontSize', fig_settings.ftsize)
xlabel('Time (mins)', 'FontSize', fig_settings.ftsize)
ylabel('Heart Rate (bpm)', 'FontSize', fig_settings.ftsize)
box off
ax = gca;
if t_fin < 5
    ax.XAxis.TickValues = t_start:t_fin;
end
grid on
ylim(fig_settings.hr_lims)

% add annotations
% - arrow indicating time of action
[arrow_coords.x(1), arrow_coords.y(1)] = axis_coords_to_fig_coords((t_start-t_start)/60,fig_settings.action_label_hr_height);
[arrow_coords.x(2), arrow_coords.y(2)] = axis_coords_to_fig_coords((t_fin-t_start)/60,fig_settings.action_label_hr_height);
ah=annotation('doublearrow',arrow_coords.x,arrow_coords.y,'Color','k', 'LineWidth', 2);
% - label of action
[label_coords(1), label_coords(2)] = axis_coords_to_fig_coords(0.3,fig_settings.action_label_hr_height);
dim = [label_coords .1 .1];
str = fig_settings.action;
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', fig_settings.ftsize-4, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% add pulse waves
ref_ax = gca;
sig_times = fig_settings.pw_start:fig_settings.pw_gap:t_fin+fig_settings.hr_durn; % in secs
durn = 10;
ftsize_small = fig_settings.ftsize - 8;
ftsize_v_small = ftsize_small-3;
pw_plot_width = 0.08;
pw_plot_height = 0.24;
start_trs = [2,2,4,1];
for sample_no = 1 : length(sig_times)
    curr_start_ms = 1000*sig_times(sample_no);
    curr_end_ms = curr_start_ms + 1000*durn;
    rel_sig_els = find(data.t_ms_rel>=curr_start_ms & data.t_ms_rel <= curr_end_ms);
    rel_trs = data.ppg_ir_tr(data.ppg_ir_tr>=rel_sig_els(1) & data.ppg_ir_tr<=rel_sig_els(end));
%     start_tr = start_trs(sample_no);
    start_tr = 1;
    rel_sig.v = data.ppg_ir(rel_trs(start_tr):rel_trs(start_tr+1));
    rel_sig.t = data.t_ms_rel(rel_trs(start_tr):rel_trs(start_tr+1)); rel_sig.t = (rel_sig.t-rel_sig.t(1))/1000;
    %     figure, plot(rel_sig.t,rel_sig.v)
    %     close(gcf)
    
    % plot this pulse wave
    time_of_pw = (sig_times(sample_no)-t_start)/60; % in mins
    [label_coords(1), ~] = axis_coords_to_fig_coords(time_of_pw,fig_settings.hr_lims(2),ref_ax);
    label_coords(2) = 0.75;
    label_coords(1) = label_coords(1)-0.5*pw_plot_width;
    hold on, axes('Position', [label_coords,pw_plot_width,pw_plot_height])
    plot(rel_sig.t, rel_sig.v, 'b', 'LineWidth', 2)
    
    % - tidy up
    ylim([min(rel_sig.v)-0.1*range(rel_sig.v), max(rel_sig.v)+0.1*range(rel_sig.v)])
    xlim([rel_sig.t(1), rel_sig.t(end)])
    xlabel('Time (s)', 'FontSize', ftsize_small)
    if sample_no == 1
        ylabel('PPG (au)', 'FontSize', ftsize_small)
    end
    set(gca, 'FontSize', ftsize_v_small, 'YTick', [])
    if rel_sig.t(end)<1
        set(gca,'XTick', 0:0.25:rel_sig.t(end))
    else
        set(gca,'XTick', 0:0.5:rel_sig.t(end))
    end
    box off
    
    % - label of time of pulse wave
    dim = [label_coords(1)+(pw_plot_width/2) label_coords(2)+(pw_plot_height*3/4) .1 .1];
    str = [num2str(round(time_of_pw)) ' mins'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize_v_small, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

end

end

function [fig_x, fig_y] = axis_coords_to_fig_coords(ax_x,ax_y, ref_ax)

% adapted from: https://uk.mathworks.com/matlabcentral/answers/166823-converting-plot-coordinates-to-normalized-coordinates

% get axis limits in fig coords
if nargin<3
    temp = get(gca,'Position');
else
    temp = get(ref_ax,'Position');
end
axis_edges_in_fig_coords.x = [temp(1), temp(1)+temp(3)];
axis_edges_in_fig_coords.y = [temp(2), temp(2)+temp(4)];
clear temp

% get axis limits in axis coords
if nargin<3
    axis_edges_in_axis_coords.x = get(gca,'XLim');
    axis_edges_in_axis_coords.y = get(gca,'YLim');
else
    axis_edges_in_axis_coords.x = get(ref_ax,'XLim');
    axis_edges_in_axis_coords.y = get(ref_ax,'YLim');
end

% convert desired axis position into fig coordinates
fig_x = axis_edges_in_fig_coords.x(1) + ...
    ((ax_x-axis_edges_in_axis_coords.x(1))/range(axis_edges_in_axis_coords.x))* ...
    (axis_edges_in_fig_coords.x(2)-axis_edges_in_fig_coords.x(1));
fig_y = axis_edges_in_fig_coords.y(1) + ...
    ((ax_y-axis_edges_in_axis_coords.y(1))/range(axis_edges_in_axis_coords.y))* ...
    (axis_edges_in_fig_coords.y(2)-axis_edges_in_fig_coords.y(1));

end

function filter_data(up)

fprintf('\n --- Filtering Data: ')

file_data = load(up.paths.conv_data_file, 'data');

sigs = {'red', 'ir'};
for sig_no = 1 : length(sigs)
    curr_sig = sigs{sig_no};
    fprintf([curr_sig ' ']);
    eval(['s.v = file_data.data.ppg_' curr_sig ';']);
    s.fs = file_data.data.fs;
    s_evlf = elim_vlfs(s, up);
    s_filt = elim_vhfs(s_evlf, up);
    eval(['data.ppg_' curr_sig ' = round(s_filt.v*10)/10;']);
    clear s s_evlf
end

data.fs = file_data.data.fs;
data.t_ms_rel = file_data.data.t_ms_rel;

% save data
fprintf('saving');
clear file_data
save(up.paths.filt_data_file, 'data', '-v7.3')

end

function assess_signal_quality(up)

fprintf('\n --- Assessing signal quality: ')

%% Load data
fprintf('\n - Loading data')
% Load PPG data
load(up.paths.beat_data_file);
data = rmfield(data, 'ppg_red');
data = rmfield(data, 'ppg_red_p');
data = rmfield(data, 'v');
% Load ADL data
load(up.paths.conv_data_file, 'adl_data')
% Load recording periods
load(up.paths.rec_periods_file);

%% Assess quality during each recording period

% setup
options.do_plot = 0;                  %  - (1 or 0, default value of 1) A logical indicating whether or not to make plots
% options.exclude_low_quality_data %  - (default of 1) A logical indicating whether or not to exclude low quality pulse waves from the calculation of median values of pulse wave indices
% options.plot_third_deriv         %  - (default of 0) A logical indicating whether or not to plot the third derivative
% options.close_figures            %  - (default of 1) A logical indicating whether or not to close all Matlab figures before running the analysis
options.do_filter = 0;                %  - (default of 1) A logical indicating whether or not to filter pulse waves prior to analysis
% options.save_folder              %  - (default of '') A string containing the path of a directory in which to save plots
% options.save_file                %  - (default of '') A string containing the filename under which to save plots (without an extension)
% options.beat_detector            %  - (default of 'IMS') The type of beat detector to use when analysing a prolonged recording of multiple pulse waves. Definitions: IMS - [incremental merge segmentation algorithm](https://doi.org/10.1109/EMBC.2012.6346628) (implemented by M.A.F Pimentel as part of the [RRest](https://github.com/peterhcharlton/RRest) toolbox).
% options.verbose                  %  - (default of 0) A logical indicating whether or not to provide text outputs for any warnings
% options.plot_areas               %  - (default of 0) A logical indicating whether or not to plot the systolic and diastolic areas on the pulse wave
% options.plot_pw_only             %  - (default of 0) A logical indicating whether or not to plot only the pulse wave and not its derivatives
% options.normalise_pw             %  - (default of 1) A logical indicating whether or not to normalise the pulse wave to occupy a range of 0 to 1
options.calc_pw_inds = 0;

% variables to store results in
[qual.no_onsets, qual.no_hq_onsets, qual.durn_hq_onsets, qual.durn_lq_onsets] = deal( nan(length(rec_periods.deb),1) );

% cycle through periods
for period_no = 1 : length(rec_periods.deb)
    
    if rem(period_no, 10) == 0
        fprintf('\n    - Period %d of %d', period_no, length(rec_periods.deb));
    end
    
    % skip if there was no connection during this period
    if ~rec_periods.connected(period_no)
        continue
    end
    
    % identify the data during this period
    period_deb = rec_periods.deb(period_no); % in ms
    period_fin = rec_periods.fin(period_no); % in ms
    rel_els = data.t_ms_rel >= period_deb & data.t_ms_rel <= period_fin;
    
    % extract the data for this period
    period_data.v = data.ppg_ir(rel_els);
    period_data.fs = data.fs;
    
    if isempty(period_data.v) || (length(period_data.v)/period_data.fs) < 5
        continue
    end
    
    % assess signal quality
    [~, ~, pulses, ~] = PulseAnalyse(period_data, options);
    pulses.durn_samps = diff(pulses.onsets);
    if ~isempty(pulses.durn_samps)
        pulses.durn_samps = [pulses.durn_samps; pulses.durn_samps(end)];
        
        % extract results to keep
        qual.no_onsets(period_no) = length(pulses.quality);
        qual.no_hq_onsets(period_no) = mean(pulses.quality);
        qual.durn_hq_onsets(period_no) = sum(pulses.durn_samps(pulses.quality));
        qual.durn_lq_onsets(period_no) = sum(pulses.durn_samps(~pulses.quality));
    end
    clear period_* rel_els pulses
end

% collate and save results
save(up.paths.quality_data_file, 'qual');

end

function identify_beats_in_data(up)

fprintf('\n --- Identifying beats in data: ')

load(up.paths.filt_data_file);
data.ppg_red = -1*data.ppg_red;
data.ppg_ir = -1*data.ppg_ir;

sigs = {'ir','red'};
for sig_no = 1 : length(sigs)
    curr_sig = sigs{sig_no};
    fprintf([curr_sig ' ']);
    eval(['data.v = data.ppg_' curr_sig ';']);
    win_indices = find_win_indices(data, up);
    p = [];
    f = waitbar(0,'Identifying beats');
    for win_no = 1 : length(win_indices.deb)
        if win_indices.complete_log(win_no)
            if rem(win_no,1000) == 0
                waitbar(win_no/length(win_indices.deb),f);
            end
            rel_el1 = win_indices.deb(win_no);
            rel_el2 = win_indices.deb(win_no)+win_indices.no_samps_in_win;
            s.v = data.v(rel_el1:rel_el2);
            s.fs = data.fs;
            try
                win_p = abd_algorithm2(s);
            catch
                win_p = [];
            end
            win_p = win_p+rel_el1-1;
            p = [p; win_p(:)];
            % plot_els = win_p-rel_el1+1; plot(s.v), hold on, plot(plot_els, s.v(plot_els), '*r')
            clear s win_p rel_el1 rel_el2
        end
    end
    close(f)
    clear win_no
    p = unique(p);
    
    eval(['data.ppg_' curr_sig '_p = p;']);
    clear s p f
end

% save data
fprintf('saving');
save(up.paths.beat_data_file, 'data', '-v7.3')

end

function win_indices = find_win_indices(data, up)

win_indices.no_samps_in_win = round(data.fs*up.beat_detect.win_length);
win_indices.no_samps_bet_wins = round(data.fs*up.beat_detect.win_length*up.beat_detect.win_step_prop);
win_start_inds = downsample(1:length(data.t_ms_rel), win_indices.no_samps_bet_wins);
win_indices.deb = win_start_inds(1:end-1);
win_indices.complete_log = false(size(win_indices.deb));
for win_no = 1 : length(win_indices.deb)
    rel_el1 = win_indices.deb(win_no);
    rel_el2 = win_indices.deb(win_no)+win_indices.no_samps_in_win;
    curr_win_durn = (data.t_ms_rel(rel_el2) - data.t_ms_rel(rel_el1));
    if curr_win_durn == (up.beat_detect.win_length*1000)
        win_indices.complete_log(win_no) = true;
    end
    
end

end

function calculate_trial_outcomes(up)

fprintf('\n --- Calculating trial outcomes')

fprintf('\n - Loading data')

%% Load ADL data
load(up.paths.conv_data_file, 'adl_data')

%% Load quality file
load(up.paths.quality_data_file)

%% Load recording periods file
load(up.paths.rec_periods_file)

%% Calculate values (which are later used to calculate metrics)

% - no rest days
no_rest_days = 4;
% - total study time
total_study_time = (rec_periods.fin(end) - rec_periods.deb(1))/1000; % in secs
% - total time for which the participant noted they were not wearing the sensor
total_noted_not_wearing = sum(rec_periods.durn(~rec_periods.noted_wearing));
% - total signal time (time for which a PPG signal was acquired from a participant)
total_signal_time = sum(rec_periods.durn(rec_periods.connected));

% - frequency of reasons for not wearing the sensor
temp_reasons = rec_periods.reason_disconnect(~rec_periods.connected);
temp_reasons = categorical(temp_reasons);
reasons_not_worn.reason = categories(temp_reasons);
reasons_not_worn.freq = countcats(temp_reasons);
reasons_not_worn.total_durn = nan(length(reasons_not_worn.reason),1);
for s = 1 : length(reasons_not_worn.reason)
    reasons_not_worn.total_durn(s) = sum(rec_periods.durn(strcmp(rec_periods.reason_disconnect, reasons_not_worn.reason{s})));
end
clear temp_reasons s
% (used code from: https://uk.mathworks.com/matlabcentral/answers/115838-count-occurrences-of-string-in-a-single-cell-array-how-many-times-a-string-appear )

% - proportion of recorded time for which signal was of high quality (according to activity)
activities = unique(rec_periods.act);
total_durn_samps = sum(qual.durn_hq_onsets(~isnan(qual.durn_hq_onsets)))+sum(qual.durn_lq_onsets(~isnan(qual.durn_lq_onsets)));
total_durn = total_durn_samps/100; % in secs
prop_time_hq.overall.t = total_durn; % in secs
prop_time_hq.overall.prop = (sum(qual.durn_hq_onsets(~isnan(qual.durn_hq_onsets)))/100)/prop_time_hq.overall.t;
prop_time_hq.acts.t = nan(length(activities),1); % in secs
prop_time_hq.acts.prop = nan(length(activities),1);
prop_time_hq.acts.n = nan(length(activities),1);
prop_time_hq.acts.act = activities;
for act_no = 1 : length(activities)
    curr_act = activities{act_no};
    rel_els = strcmp(rec_periods.act, curr_act);
    total_durn_samps = sum(qual.durn_hq_onsets(rel_els & ~isnan(qual.durn_hq_onsets)))+sum(qual.durn_lq_onsets(rel_els & ~isnan(qual.durn_lq_onsets)));
    total_durn = total_durn_samps/100; % in secs
    prop_time_hq.acts.t(act_no) = total_durn; % in secs
    prop_time_hq.acts.n(act_no) = sum(rel_els); % in secs
    prop_time_hq.acts.prop(act_no) = (sum(qual.durn_hq_onsets(rel_els & ~isnan(qual.durn_hq_onsets)))/100)/prop_time_hq.acts.t(act_no);
end
clear act_no activities total_durn total_durn_samps

%% Calculate metrics

% - total possible wear time (4 weeks, 6 days per week, in secs)
metrics.TPWT = total_study_time - (no_rest_days*24*60*60);

% - total wear time (secs)
metrics.TWT = total_study_time - total_noted_not_wearing;

% - total signal time (time for which a PPG signal was acquired from a participant)
metrics.TST = total_signal_time;

% - Percentage of TPWT for which the participant did not wear the sensor due to each reason
reasons_not_worn.perc_tpwt = 100*reasons_not_worn.total_durn./metrics.TPWT;

%% Calculate endpoints

fprintf('\n ~~~ Endpoints ~~~')

% - total possible wear time
endpoints.TPWT = metrics.TPWT/(60*60*24);
fprintf('\n - Total Possible Wear Time: %.1f days (%.1f hours, %.1f %%)', endpoints.TPWT, endpoints.TPWT*24, 100)

% - length of time for which the participant noted they wore the sensor
endpoints.duration_TWT = metrics.TWT/(24*60*60);
endpoints.TWT_as_perc_of_TPWT = 100*metrics.TWT/metrics.TPWT;
fprintf('\n - Wear time: %.1f days (%.1f hours, %.1f %%)', endpoints.duration_TWT, endpoints.duration_TWT*24, endpoints.TWT_as_perc_of_TPWT)

% - frequency of reasons for not wearing the sensor
endpoints.reasons_for_not_wearing = reasons_not_worn;

% - total signal time
endpoints.TST = metrics.TST/(60*60);
endpoints.TST_as_perc_of_TPWT = 100*metrics.TST/metrics.TPWT;
fprintf('\n - Signal time: %.1f days (%.1f hours, %.1f %%)', endpoints.TST/24, endpoints.TST, endpoints.TST_as_perc_of_TPWT)

% - Percentage of TWT for which a PPG signal was acquired (i.e. 100 x TST / TWT)
endpoints.TST_as_perc_of_TWT = 100*metrics.TST/metrics.TWT;
fprintf('\n - Percentage of wear time for which a PPG signal was acquired: %.1f %%', endpoints.TST_as_perc_of_TWT)

% Total high quality signal time
endpoints.THQT = prop_time_hq.overall.t*prop_time_hq.overall.prop;
endpoints.THQT_as_perc_of_TPWT = 100*endpoints.THQT/metrics.TPWT;
fprintf('\n - High quality signal time: %.1f days (%.1f hours, %.1f %%)', endpoints.THQT/(60*60*24), endpoints.THQT/(60*60), endpoints.THQT_as_perc_of_TPWT)

% - Percentage of TST (near enough) for which the PPG signal was of high quality
endpoints.prop_time_hq = prop_time_hq;
fprintf('\n - Percentage of TST for which PPG was high quality: %.1f %%', prop_time_hq.overall.prop*100)

% - Number of Bluetooth disconnections
endpoints.no_disconnections = sum(strcmp(adl_data.act, 'disconnected_-_connection_dropped'));
fprintf('\n - There were %d Bluetooth disconnections', endpoints.no_disconnections)

%% Save trial outcomes
save(up.paths.trial_outcomes_file, 'endpoints');

end

function p = abd_algorithm2(data)
% ABD_ALGORITHM  identifies pulse peaks in a pulsatile signal. It is an
% adaptation of the algorithm described in:
%
%    Aboy, M. et al., An automatic beat detection algorithm for pressure
%    signals. IEEE Trans. Biomed. Eng. 2005, 52, 1662?1670,
%      http://doi.org/10.1109/TBME.2005.855725
%   
%  Inputs:
%
%    data     -  a pulsatile signal formatted as a structure, containing
%                 a vector of amplitudes, data.v, and the sampling
%                 frequency (in Hz), data.fs.
%    options  -  (optional) a structure of options, which may contain any of:
%                    options.                    - a logical (true or false)
%
%  Outputs:
%
%    p        -  a vector of pulse indices
%
%  Exemplary usage:
%
%    p = abd_algorithm(data)                              extracts pulse indices from the pulsatile signal, data.
%    p = abd_algorithm(data, options)                     uses options to specify the analysis.
%    [p, S_filt] = abd_algorithm(___)                     also outputs the filtered signals.
%
%  For further information please see the accompanying manual.
%
%  This script contains items either copied or modified from the RRest
%  toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%
% Peter H. Charlton, King's College London, August 2017

% Changes from original algorithm description:
%  1) PSD calculation method may not be exactly the same
%  2) Not conducted on windows of 10 s
%  3) Band-pass filtering may not produce exactly the right cut-offs
%  4) Wasn't sure what the lower and upper HR values were, so used 30 and 200 bpm
%  5) Changed the proportion of the lower HR value at which to draw the lower cut-off (from 50% to 80%)
%  6) Changed the percentile threshold for identifying peaks in the derivative from 90% to 75%
%  7) Haven't implemented harmonic PSD
%  8) HR estimation only considers HRs within the range of plausible HRs

% load data
if nargin<1
    load('/Users/petercharlton/Downloads/ppg_data.mat');
    data.v = -1*data.v;
end

% inputs
x = data.v;  % signal
fs = data.fs; % sampling freq
up = setup_up; % settings
w = fs*10; % window length (number of samples)
win_starts = 1:round(0.8*w):length(x);
win_starts(win_starts>=length(x)-w+1) = [];
win_starts = [win_starts, length(x)+1-w];

% before pre-processing
px = DetectMaxima(x,0);  % detect all maxima

% detect peaks in windows
all_p4 = [];
all_hr = nan(length(win_starts)-1,1);
for win_no = 1 : length(win_starts)-1
    curr_els = win_starts(win_no):win_starts(win_no)+w-1;
    curr_x = x(curr_els);
    y1 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 3*up.fh_hz); % Filter no.1
    hr = EstimateHeartRate(y1, fs, up); % Estimate HR from weakly filtered signal
    all_hr(win_no) = hr;
    if hr < 40
        a = 1;
    end
    y2 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 2.5*hr/60); % Filter no.2
    y2_deriv = EstimateDeriv(y2); % Estimate derivative from highly filtered signal
    p2 = DetectMaxima(y2_deriv,up.deriv_threshold); % Detect maxima in derivative
    % plot(x), hold on, plot(p2, x(p2), 'or')
    % plot(y2_deriv), hold on, thresh = prctile(y2_deriv, up.deriv_threshold); plot([0,length(y2_deriv)], thresh*[1,1])
    y3 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 10*hr/60);
    p3 = DetectMaxima(y3,60); % Detect maxima in moderately filtered signal
    % plot(x), hold on, plot(p3, x(p3), 'or')
    p4 = find_pulse_peaks(p2,p3);
%     plot(curr_x), hold on, plot(p4, curr_x(p4), 'or')
    all_p4 = [all_p4;win_starts(win_no)+p4-1];
end

all_p4 = unique(all_p4);

[p, fn] = IBICorrect(all_p4, px, median(all_hr), fs, up);
p = unique(p);
% plot(x), hold on, plot(p, x(p), 'or')

% plot(y1), hold on, plot(p4, y1(p4), 'or')
p = all_p4;

end

function up = setup_up

% plausible HR limits
up.fl = 30; % lower bound for HR
up.fh = 200; % upper bound for HR
up.fl_hz = up.fl/60;
up.fh_hz = up.fh/60;

% Thresholds
up.deriv_threshold = 75;            % originally 90
up.upper_hr_thresh_prop = 2.25;     % originally 1.75
up.lower_hr_thresh_prop = 0.5;     % originally 0.75

% Other parameters
up.win_size = 10; % in secs

end

function mt = DetectMaxima(sig,percentile)

% Table VI pseudocode

tr = prctile(sig, percentile);
ld = length(sig);

m = 1+find(sig(3:end) < sig(2:end-1) & sig(1:end-2) < sig(2:end-1));
mt = m(sig(m)>tr);

end

function bpf_sig = Bandpass(sig, fs, lower_cutoff, upper_cutoff)

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 1.3*lower_cutoff;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.8*lower_cutoff;   % in Hz
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 1.2*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Fstop = 0.8*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.03;

% perform BPF
s.v = sig;
s.fs = fs;
s_evlf = elim_vlfs(s, up);
s_filt = elim_vhfs(s_evlf, up);
bpf_sig = s_filt.v;

end

function s_filt = elim_vlfs(s, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vlf.Fstop up.paramSet.elim_vlf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vlf.Dstop up.paramSet.elim_vlf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

try
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
    s_filt.v = s.v-s_filt.v;
catch
    s_filt.v = s.v;
end
s_filt.fs = s.fs;

end

function s_filt = elim_vhfs(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.3355;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);
end

function hr = EstimateHeartRate(sig, fs, up)

% Estimate PSD
blackman_window = blackman(length(sig), 'periodic');
[pxx,f] = periodogram(sig, blackman_window,length(sig), fs);
% [ph, fh] = harmonicPSD(pxx,f);
ph = pxx; fh = f;

% Extract HR
rel_els = fh >= up.fl_hz & fh <= up.fh_hz;
rel_p = ph(rel_els);
rel_f = fh(rel_els);
[~,max_el] = max(rel_p);
hr = rel_f(max_el)*60;

end

function [ph, fh] = harmonicPSD(pxx,f)

% settings
n = 11;
alpha = 2;

ph = nan(length(f),1);
for freq_no = 1 : length(f)
    temp = 0;
    for k = 1 : n
        if freq_no*k > length(f)
            harmonic_p = 0;
        else
            harmonic_p = pxx(freq_no*k);
        end
        temp = temp + min([alpha*pxx(freq_no), harmonic_p]);
    end
    ph(freq_no) = temp;
end

fh = f;

end

function deriv = EstimateDeriv(sig)

% Savitzky Golay
deriv_no = 1;
win_size = 5;
deriv = savitzky_golay(sig, deriv_no, win_size);

end

function deriv = savitzky_golay(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function p4 = find_pulse_peaks(p2,p3)
p4 = nan(length(p2),1);
for k = 1 : length(p2)
    rel_el = find(p3>p2(k),1);
    if ~isempty(rel_el)
        p4(k) = p3(rel_el);
    end
end
p4 = p4(~isnan(p4));
end

function [pc, fn] = IBICorrect(p, m, hr, fs, up)

% Correct peaks' location error due to pre-processing
pc = nan(length(p),1);
for k = 1 : length(p)
    [~,rel_el] = min(abs(m-p(k)));
    pc1(k) = m(rel_el);    
end

% Correct false positives
% - identify FPs
d = diff(pc1)/fs;  % interbeat intervals in secs
fp = find_reduced_IBIs(d, median(hr), up);
% - remove FPs
pc2 = pc1(~fp);

% Correct false negatives
d = diff(pc2)/fs;  % interbeat intervals in secs
fn = find_prolonged_IBIs(d, median(hr), up);

pc = pc1;

end

function fn = find_prolonged_IBIs(IBIs, med_hr, up)

IBI_thresh = up.upper_hr_thresh_prop*60/med_hr;
fn = IBIs > IBI_thresh;

end

function fp = find_reduced_IBIs(IBIs, med_hr, up)

IBI_thresh = up.lower_hr_thresh_prop*60/med_hr;
fp = IBIs < IBI_thresh;

end
