function result = run_arousal_pipeline_v3(filepath, n3_threshold_override)
% run_arousal_pipeline_v3
% reads one night of watch data, counts the arousals the watch logged (real),
% and simulates the watch's stimulation rules sample by sample (predicted).
% inputs
%   filepath: full path to one night's csv file
%   n3_threshold_override:this night's n3 threshold from the audit (0.20 or 0.30)

%% settings
motion_threshold = 3;       % gyro |x|+|y|+|z| 
max_stim_count = 3600;    % max stim per night, in watch samples
cooldown_s = 60;      % seconds after an arousal before stim may start again
window_min_samples = 1800;    % stim only after watch sample 1800 (~30 min)
window_max_samples = 14400;   % stim before watch sample 14400 (~4 h)
stable_n_samples= 60;      % n3 must stay above threshold for the last 60 watch samples
gap_min_mult = 1.8;     % a gap of >= 1.8 normal intervals means a log line was dropped
gap_max_mult = 2.5;     % only gaps < 2.5 intervals hold samples the watch counted
long_gap_s = 300;     % a gap over 5 min = a long stretch of missing data
switch_date = datetime(2025, 3, 14);   % n3 threshold 0.20 before this date, 0.30 after

%% step 1: load the file 
% sensor_type and value as text
opts = detectImportOptions(filepath, 'FileType', 'text', 'Delimiter', '\t');
opts = setvartype(opts, {'sensor_type', 'value'}, 'string');
data = readtable(filepath, opts);

% empty text cells become "" 
data.value(ismissing(data.value)) = "";
data.sensor_type(ismissing(data.sensor_type)) = "";

% remove rows that appear more than once 
n_before = height(data);
data = unique(data, 'stable');
n_after = height(data);
if n_after < n_before
    fprintf('  Deduplicated: %d -> %d rows (%d duplicates removed)\n', n_before, n_after, n_before - n_after);
end

%% step 2: participant, condition and date from the file name
% file names look like: 1024634152_stimType1_2025-02-09.csv
[~, name, ~] = fileparts(filepath);
parts = split(name, "_");
userID = string(parts(1));
stimType= lower(string(parts(2)));      % "stimType1" -> "stimtype1"
night_date = string(parts(3));

if stimType == "nostim"
    stimType = "stimtype3";
elseif ~ismember(stimType, ["stimtype1", "stimtype2", "stimtype3"])
    stimType = "unknown";                  % the batch uses the audited condition anyway
end

%% step 3: n3 threshold for this night 
if datetime(night_date, 'InputFormat', 'yyyy-MM-dd') < switch_date
    n3_threshold = 0.20;
else
    n3_threshold = 0.30;
end

% if n3_threshold doesnt fit into the set date above
if nargin >= 2 && ~isempty(n3_threshold_override) && ~isnan(n3_threshold_override)
    n3_threshold = n3_threshold_override;
end

%% step 4: real arousals 
arousal_rows = data(data.sensor_type == "arousal", :);
n_arousals   = height(arousal_rows);

%% step 5: sleep-algorithm samples (MITSleep)
mit_rows = data(data.sensor_type == "MITSleep", :);
if height(mit_rows) < 2
    error('fewer than 2 MITSleep rows');
end
mit_rows = sortrows(mit_rows, 'timestamp');    % put samples in time order
mit_ts = double(mit_rows.timestamp);         % time of each sample (ms)
n_mit = length(mit_ts);

% remove brackets and spaces, split at commas
mit_text  = replace(mit_rows.value, ["[", "]", " "], "");
mit_parts = split(mit_text, ",");

inst_N3 = str2double(mit_parts(:, 4));        % field 4:  n3 probability at this sample
avg_N3 = str2double(mit_parts(:, 9));        % field 9:  watch's running average of field 4
real_on = mit_parts(:, 12) == "true";         % field 12: was the watch stimulating?

n_fields = size(mit_parts, 2);
buffer_n = nan(n_mit, 1);                     % field 11: samples in the watch's average buffer
if n_fields >= 11
    buffer_n = str2double(mit_parts(:, 11));
end
counter = nan(n_mit, 1);                      % field 13: watch's own stim counter
if n_fields >= 13
    counter = str2double(mit_parts(:, 13));
end

%% step 6: gyroscope (movement)
gyro_rows = data(data.sensor_type == "LSM6DSO Gyroscope", :);
gyro_ts = double(gyro_rows.timestamp);
gyro_text = replace(gyro_rows.value, ["[", "]", " "], "");
gyro_xyz= str2double(split(gyro_text, ","));
gyro_sum = abs(gyro_xyz(:, 1)) + abs(gyro_xyz(:, 2)) + abs(gyro_xyz(:, 3));   % the watch's gyroTotal

% put the gyro readings in time order
[gyro_ts, order] = sort(gyro_ts);
gyro_sum = gyro_sum(order);
n_gyro   = length(gyro_ts);

%% step 7: the watch's own sample count
% the watch counts its own samples (~1.05 s apart). if one log line was dropped,
% the watch still counted that sample, so a gap of ~2 intervals = 2 samples.
% longer gaps are times the watch produced no samples, so they count as 1.

% seconds between each sample and the one before it
dt_s = zeros(n_mit, 1);
for i = 2:n_mit
    dt_s(i) = (mit_ts(i) - mit_ts(i-1)) / 1000;
end
normal_dt = median(dt_s(2:end));               % the normal sample interval

% how many watch samples each logged sample stands for
ticks = ones(n_mit, 1);
for i = 2:n_mit
    if dt_s(i) >= gap_min_mult * normal_dt && dt_s(i) < gap_max_mult * normal_dt
        ticks(i) = round(dt_s(i) / normal_dt);
    end
end

% running watch sample count
tick_no = zeros(n_mit, 1);
tick_no(1) = ticks(1);
for i = 2:n_mit
    tick_no(i) = tick_no(i-1) + ticks(i);
end

%% step 8: first sample of each 60-sample window
% w_first(i) = earliest logged sample within the last 60 watch samples
w_first = ones(n_mit, 1);
j = 1;
for i = 1:n_mit
    while tick_no(j) <= tick_no(i) - stable_n_samples
        j = j + 1;
    end
    w_first(i) = j;
end

%% step 9: movement, checked once per sample like the watch 
% at each sample the watch looks only at its latest gyro reading
motion_flag = false(n_mit, 1);
j = 0;                                         % index of the latest gyro reading so far
for i = 1:n_mit
    % move forward to the last gyro reading at or before this sample
    while j < n_gyro && gyro_ts(j+1) <= mit_ts(i)
        j = j + 1;
    end
    if j > 0
        % several readings can share one timestamp: use the largest of them
        biggest = gyro_sum(j);
        k = j;
        while k > 1 && gyro_ts(k-1) == gyro_ts(j)
            k = k - 1;
            biggest = max(biggest, gyro_sum(k));
        end
        motion_flag(i) = biggest > motion_threshold;
    end
end

% recent movement = movement at any sample in the last 60 watch samples
recent_motion = false(n_mit, 1);
for i = 1:n_mit
    recent_motion(i) = any(motion_flag(w_first(i):i));
end

%% step 10: app restarts
% after a restart the watch's average buffer (field 11) starts again near 0.
% the watch did not stimulate again after a restart.
after_restart = false(n_mit, 1);
for i = 2:n_mit
    if buffer_n(i) - buffer_n(i-1) < -10
        after_restart(i:end) = true;
        break;
    end
end

%% step 11: simulate the watch, sample by sample 
seconds_after_start = (mit_ts - mit_ts(1)) / 1000;

stim_on      = false;                  % is the simulated watch stimulating right now?
stim_ticks   = 0;                      % watch samples of stim so far (for the nightly cap)
pred_stim    = false(n_mit, 1);        % simulated stim at each sample
pred_arousal = false(n_mit, 1);        % simulated arousal (stim stopped by n3 or movement)
last_arousal = -Inf;                   % time of the last simulated arousal (s)
fail_reason  = strings(n_mit, 1);      % why stim could not start (diagnostics only)

for i = 1:n_mit
    % this sample stops stim if n3 is low or there was recent movement
    disqualified = inst_N3(i) < n3_threshold || avg_N3(i) < n3_threshold || recent_motion(i);

    % stim is only allowed strictly between watch samples 1800 and 14400
    in_window = tick_no(i) > window_min_samples && tick_no(i) < window_max_samples;

    if ~stim_on
        % ----- stim is off: can it start? -----
        if stim_ticks >= max_stim_count
            fail_reason(i) = "max_stim_reached";
        elseif after_restart(i)
            fail_reason(i) = "after_app_restart";
        elseif ~in_window
            fail_reason(i) = "outside_time_window";
        elseif seconds_after_start(i) - last_arousal < cooldown_s
            fail_reason(i) = "within_cooldown";
        else
            % n3 must be at or above threshold across the whole window
            % (the average is checked over one extra, earlier sample)
            first   = w_first(i);
            inst_ok = all(inst_N3(first:i) >= n3_threshold);
            avg_ok  = all(avg_N3(max(first - 1, 1):i) >= n3_threshold);

            if ~(inst_ok && avg_ok)
                fail_reason(i) = "not_stable_60s";
            elseif recent_motion(i)
                fail_reason(i) = "recent_motion";
            else
                stim_on = true;                        % start stim
                pred_stim(i) = true;
                stim_ticks = stim_ticks + ticks(i);
            end
        end
    else
        % stim is on: does it stop?
        if tick_no(i) >= window_max_samples || after_restart(i)
            stim_on = false;                           % window closed or app restarted (no arousal logged)
        elseif disqualified
            stim_on = false;                           % arousal: n3 dropped or movement
            pred_arousal(i) = true;
            last_arousal = seconds_after_start(i);
        elseif stim_ticks >= max_stim_count
            stim_on = false;                           % nightly cap reached
        else
            pred_stim(i) = true;                       % keep stimulating
            stim_ticks = stim_ticks + ticks(i);
        end
    end
end

%% step 12: results 
result.userID          = userID;
result.stimType        = stimType;
result.night_date      = night_date;
result.recording_hours = (max(mit_ts) - min(mit_ts)) / (1000 * 60 * 60);

% real: what the watch logged
result.stim_on_seconds = sum(real_on);
result.n_arousals      = n_arousals;
if result.stim_on_seconds > 0
    result.arousal_rate_per_min = n_arousals / (result.stim_on_seconds / 60);
else
    result.arousal_rate_per_min = NaN;
end
result.n3_threshold = n3_threshold;
if any(real_on)
    result.min_on_inst = min(inst_N3(real_on));       % lowest n3 the watch tolerated during stim
else
    result.min_on_inst = NaN;
end

% predicted: what the simulation did
result.predicted_stim_seconds = sum(pred_stim);
result.predicted_n_arousals   = sum(pred_arousal);
if result.predicted_stim_seconds > 0
    result.predicted_arousal_rate = result.predicted_n_arousals / (result.predicted_stim_seconds / 60);
else
    result.predicted_arousal_rate = NaN;
end

% did the watch start stim right after a long logging gap? then the watch decided using
% data missing from the file, so the night is left out of the real-vs-predicted check
result.onset_after_long_gap = false;
for i = 2:n_mit
    if real_on(i) && ~real_on(i-1)                    % the watch started stimulating here
        first = max(i - stable_n_samples, 1) + 1;
        if any(dt_s(first:i) > long_gap_s)
            result.onset_after_long_gap = true;
        end
    end
end

% per-sample detail (used by the diagnostic plots; the batch ignores it)
result.detail = struct('mit_ts', mit_ts, 'tick_no', tick_no, 'inst', inst_N3, 'avg', avg_N3, ...
    'real_on', real_on, 'sim_on', pred_stim, 'sim_arousal', pred_arousal, ...
    'fail_reason', fail_reason, 'recent_motion', recent_motion, 'counter', counter, ...
    'buffer_n', buffer_n, 'arousal_ts', double(arousal_rows.timestamp));
end