function result = run_arousal_pipeline_v2(filepath, n3_threshold_override, sim_opts)

%% Arousal Rate calculation
% clear;
% filepath = 'D:\csv_extracted_v2\all_csvs\1024634152_stimType1_2025-02-09.csv';

%load the CSV
opts = detectImportOptions(filepath, 'FileType', 'text', 'Delimiter', '\t');
opts = setvartype(opts, {'sensor_type', 'value'}, 'string');    % never let readtable guess numeric
data = readtable(filepath, opts);
data.value(ismissing(data.value)) = "";          % empty fields as "" so duplicates match in unique()
data.sensor_type(ismissing(data.sensor_type)) = "";

%de duplicate time stamps
n_before = height(data);
data = unique(data, 'stable');   % remove fully-identical rows, preserve order
n_after = height(data);
if n_before ~= n_after
    fprintf('  Deduplicated: %d -> %d rows (%d duplicates removed)\n', ...
        n_before, n_after, n_before - n_after);
end

%extract info from filename
[~, name, ~] = fileparts(filepath);
parts = split(name, "_");

userID  = string(parts(1));
stimType = string(parts(2));
night_date = string(parts(3));

% normalize stimType - handle case variation and unknown labels
stimType_normalized = lower(stimType);   % "stimType1" -> "stimtype1"
if ismember(stimType_normalized, ["stimtype1", "stimtype2", "stimtype3"])
    stimType = stimType_normalized;
elseif ismember(stimType_normalized, "nostim")
    stimType = "stimtype3";
else
    stimType = "unknown";       % condition comes from the audit's stimStatus, not the filename
end

arousal_events = table();

% motion parameters
motion_window_ms = 60000;  % lookback window for recent motion (default rule)

% simulation options - read BEFORE motion is computed (defaults = current rules)
if nargin < 3 || isempty(sim_opts), sim_opts = struct(); end
if ~isfield(sim_opts, 'motion_metric'),         sim_opts.motion_metric         = "abssum";       end   % confirmed by Nathan: |x|+|y|+|z|  ("norm" = vector magnitude)
if ~isfield(sim_opts, 'motion_threshold'),      sim_opts.motion_threshold      = 3;              end   % raw gyroTotal > 3
if ~isfield(sim_opts, 'motion_stops_stim'),     sim_opts.motion_stops_stim     = true;           end
if ~isfield(sim_opts, 'motion_eval'), sim_opts.motion_eval = "epoch_latest"; end   % watch checks the latest gyro reading once per MIT sample
if ~isfield(sim_opts, 'no_stim_after_restart'), sim_opts.no_stim_after_restart = true;          end
   if ~isfield(sim_opts, 'tick_rule'), sim_opts.tick_rule = "short_only"; end   % only single dropped log lines count as watch samples
   if ~isfield(sim_opts, 'strict_window'), sim_opts.strict_window = true; end   % watch stims only for 1800 < count < 14400
if ~isfield(sim_opts, 'tick_max_mult'), sim_opts.tick_max_mult = 2.5; end   % gaps shorter than this many intervals count as missing samples  
motion_threshold = sim_opts.motion_threshold;

%% 3 finding the arousal and stim events
arousal_events = data(data.sensor_type == "arousal", :);
stim_events    = data(data.sensor_type == "stim", :);
stim_ts        = stim_events.timestamp;

%% 4 extract the p(n3) from MITsleep
%filter to MITSleep rows and put them in time order
mit_rows = data(data.sensor_type == "MITSleep", :);
if height(mit_rows) < 2, error('fewer than 2 MITSleep rows'); end
mit_rows = sortrows(mit_rows, 'timestamp');
mit_ts = double(mit_rows.timestamp);
n_mit  = length(mit_ts);

%remove the brackets from each value string
mit_clean = replace(mit_rows.value, ["[", "]", " "], "");
mit_parts = split(mit_clean, ",");

%extract the fields we need
mit_inst_N3 = str2double(mit_parts(:, 4));
mit_avg_N3  = str2double(mit_parts(:, 9));
mit_stim_on = mit_parts(:, 12) == "true";      %was stim presented yes/no

%% 5 extract the gyroscope data
gyro_rows = data(data.sensor_type == "LSM6DSO Gyroscope", :);
gyro_ts   = gyro_rows.timestamp;
gyro_clean = replace(gyro_rows.value, ["[", "]", " "], "");
gyro_parts = split(gyro_clean, ",");
gyro_xyz_rads = str2double(gyro_parts);

% CHANGED: motion measure now follows sim_opts.motion_metric (default |x|+|y|+|z|)
if string(sim_opts.motion_metric) == "norm"
    gyro_sum = sqrt(gyro_xyz_rads(:,1).^2 + gyro_xyz_rads(:,2).^2 + gyro_xyz_rads(:,3).^2);   % vector magnitude
else
    gyro_sum = abs(gyro_xyz_rads(:,1)) + abs(gyro_xyz_rads(:,2)) + abs(gyro_xyz_rads(:,3));   % watch: |x|+|y|+|z|
end

%% 6 which gyro samples exceed the motion threshold?
gyro_sum_axis_big = gyro_sum > motion_threshold;
big_motion_ts = gyro_ts(gyro_sum_axis_big);

%% 7 recent motion (default rule: any big motion in the last 60 s)
big_motion_ts = unique(double(big_motion_ts));                 % sorted
recent_motion = false(n_mit, 1);
if ~isempty(big_motion_ts)
    k = discretize(double(mit_ts), [big_motion_ts; Inf]);      % last big motion at or before each sample
    has = ~isnan(k);
    recent_motion(has) = double(mit_ts(has)) - big_motion_ts(k(has)) <= motion_window_ms;
end

%% 8 arousal calculation (using logged events)
n_arousals = height(arousal_events);
stim_on_seconds = sum(mit_stim_on); %count how many rows had stim on

if stim_on_seconds > 0
    arousal_rate_per_min = n_arousals / (stim_on_seconds / 60); %arousals per min of stim on
else
    arousal_rate_per_min = NaN;
end

%% 9 simulating stim periods - one arousal logged per transition from stim on -> disqualified

% parameters
max_stim_count     = 3600;    % max cumulative stim per night (watch samples)
post_arousal_gap_s = 60;      % cooldown after an arousal before stim can resume
window_min_samples = 1800;    % ~30 min of samples
window_max_samples = 14400;   % ~4 hr of samples
stable_n_samples   = 60;      % stable-N3 window = last 60 watch samples (~63.5 s)

% N3 threshold changed between app versions (0.20 before 2025-03-14, 0.30 after)
threshold_switch_date = datetime(2025, 3, 14);
if datetime(night_date, 'InputFormat', 'yyyy-MM-dd') < threshold_switch_date
    n3_threshold = 0.20;
else
    n3_threshold = 0.30;
end

% threshold from the audit (derived night date) overrides the filename-based rule
if nargin >= 2 && ~isempty(n3_threshold_override) && ~isnan(n3_threshold_override)
    n3_threshold = n3_threshold_override;
end

% setup
initial_ts = mit_ts(1);
seconds_after_start = (mit_ts - initial_ts) / 1000;

% watch tick count: dropped log lines still advance the watch's own counter.
% Only intervals >= 1.8x the normal one count as dropped samples, so logging
% jitter (one late line = one long + one short interval) is not double-counted.

   dt_s    = [0; diff(double(mit_ts)) / 1000];
   med_dt  = median(dt_s(2:end));
   long_iv = dt_s >= 1.8 * med_dt;                       % at least one missing log line
   switch string(sim_opts.tick_rule)
       case "current"                                    % gaps up to 5 min count as missing samples
           cnt = long_iv & dt_s <= 300;
       case "short_only"                                 % only short gaps count as missing watch samples
            cnt = long_iv & dt_s < sim_opts.tick_max_mult * med_dt;
       case "short_plus_long_gyro"                       % short gaps count; long gaps count only if the gyro kept logging
            ug = unique(double(gyro_ts));
            running = false(n_mit, 1);
            if numel(ug) > 1
                n_le    = discretize(double(mit_ts), [-Inf; ug; Inf]) - 1;
                g_in    = [0; diff(n_le)];
                g_rate  = numel(ug) / ((ug(end) - ug(1)) / 1000);
                running = g_in >= 0.5 * g_rate .* dt_s;
            end
           cnt = (long_iv & dt_s < sim_opts.tick_max_mult * med_dt) | (dt_s > 60 & running);
       otherwise                                         % "gyro_present" | "gyro_absent"
           ug = unique(double(gyro_ts));
           running = false(n_mit, 1);
           if numel(ug) > 1
               n_le    = discretize(double(mit_ts), [-Inf; ug; Inf]) - 1;   % gyro readings up to each MIT sample
               g_in    = [0; diff(n_le)];                                   % gyro readings inside each interval
               g_rate  = numel(ug) / ((ug(end) - ug(1)) / 1000);            % normal gyro readings per second
               running = g_in >= 0.5 * g_rate .* dt_s;                      % gyro kept logging through the gap
           end
           if string(sim_opts.tick_rule) == "gyro_present", cnt = long_iv & running;
           else,                                             cnt = long_iv & ~running; end
   end
   ticks = ones(n_mit, 1);
   ticks(cnt) = round(dt_s(cnt) / med_dt);
   tick_no = cumsum(ticks);

% first logged sample of the 60-tick stable window for each sample (tick_no only increases)
w_first = ones(n_mit, 1);
j = 1;
for i = 1:n_mit
    while tick_no(j) <= tick_no(i) - stable_n_samples
        j = j + 1;
    end
    w_first(i) = j;
end
stim_ticks = 0;                        % watch's own stim counter (drives the cap)

% motion checked once per MIT sample, looking back 60 watch samples (test variants)
if string(sim_opts.motion_eval) ~= "lookback_60s"
    flag = false(n_mit, 1);                                    % motion flagged at this MIT sample
    if ~isempty(gyro_ts)
        [ug, ~, ig] = unique(double(gyro_ts));                 % sorted unique gyro timestamps
        gs = accumarray(ig, gyro_sum, [], @max);               % largest reading at each timestamp
        kg = discretize(double(mit_ts), [ug; Inf]);            % latest gyro reading at or before each MIT sample
        kg(isnan(kg)) = 0;
        if string(sim_opts.motion_eval) == "epoch_latest"
            has = kg > 0;
            flag(has) = gs(kg(has)) > motion_threshold;        % only the latest reading counts
        else                                                   % "epoch_any"
            cb   = [0; cumsum(gs > motion_threshold)];
            prev = [0; kg(1:end-1)];
            flag = (cb(kg + 1) - cb(prev + 1)) > 0;            % any reading since the previous MIT sample
        end
    end
    cf = [0; cumsum(flag)];
    recent_motion = (cf((1:n_mit)' + 1) - cf(w_first)) > 0;   % flagged in any of the last 60 watch samples
end

% app restart = field 11 (avg-buffer size) drops back towards 0
if size(mit_parts, 2) >= 11, buffer_n = str2double(mit_parts(:, 11)); else, buffer_n = nan(n_mit, 1); end
after_restart = false(n_mit, 1);
restart_idx = find(diff(buffer_n) < -10, 1) + 1;              % first restart, if any
if sim_opts.no_stim_after_restart && ~isempty(restart_idx)
    after_restart(restart_idx:end) = true;                     % watch did not stim again after a restart
end

% state
stim_state = false;                    % current ON/OFF state
predicted_stim_seconds = 0;            % logged stim samples (comparable to real stim)
predicted_stim    = false(n_mit, 1);   % logical array (plotting + summing)
predicted_arousal = false(n_mit, 1);   % logical array (arousal events)
last_arousal_time = -Inf;              % time of last arousal (cooldown)
fail_reason = strings(n_mit, 1);       % diagnostics

% main loop
for i = 1:n_mit
    % is a disqualifier active at THIS sample?
    disqualified = (mit_inst_N3(i) < n3_threshold) || (mit_avg_N3(i) < n3_threshold) || ...
                   (recent_motion(i) && sim_opts.motion_stops_stim);

    if ~stim_state
        % stim is currently OFF: check whether it can START
               if stim_ticks >= max_stim_count
            fail_reason(i) = "max_stim_reached";
        elseif after_restart(i)
            fail_reason(i) = "after_app_restart";
        elseif (~sim_opts.strict_window && (tick_no(i) <  window_min_samples || tick_no(i) >  window_max_samples)) || ...
               ( sim_opts.strict_window && (tick_no(i) <= window_min_samples || tick_no(i) >= window_max_samples))
            fail_reason(i) = "outside_time_window";
        elseif seconds_after_start(i) - last_arousal_time < post_arousal_gap_s
            fail_reason(i) = "within_cooldown";
        elseif tick_no(i) < stable_n_samples
            fail_reason(i) = "not_stable_60s";
        else
            w     = w_first(i) : i;                              % inst: last 60 watch samples
            w_avg = max(w_first(i) - 1, 1) : i;                  % avg: one extra sample
            stable_N3 = all(mit_inst_N3(w) >= n3_threshold) && all(mit_avg_N3(w_avg) >= n3_threshold);

            if ~stable_N3
                fail_reason(i) = "not_stable_60s";
            elseif recent_motion(i)
                fail_reason(i) = "recent_motion";
            else
                stim_state = true;
                predicted_stim(i) = true;
                predicted_stim_seconds = predicted_stim_seconds + 1;
                stim_ticks = stim_ticks + ticks(i);
            end
        end
    else
        % stim is currently ON
        if tick_no(i) > window_max_samples || (sim_opts.strict_window && tick_no(i) >= window_max_samples) || after_restart(i)
            stim_state = false;                                  % stop, no arousal logged
        elseif disqualified
            predicted_arousal(i) = true;
            stim_state = false;
            last_arousal_time = seconds_after_start(i);
        elseif stim_ticks >= max_stim_count
            stim_state = false;
        else
            predicted_stim(i) = true;
            predicted_stim_seconds = predicted_stim_seconds + 1;
            stim_ticks = stim_ticks + ticks(i);
        end
    end
end

%% save real metrics and metadata to result struct
result.userID          = userID;
result.stimType        = stimType;
result.night_date      = night_date;
result.recording_hours = (max(mit_ts) - min(mit_ts)) / (1000*60*60);
result.stim_on_seconds = stim_on_seconds;
result.n_arousals      = n_arousals;
result.arousal_rate_per_min = arousal_rate_per_min;
result.n3_threshold    = n3_threshold;                         % which rule this night used
if any(mit_stim_on), result.min_on_inst = min(mit_inst_N3(mit_stim_on));
else, result.min_on_inst = NaN;
end

% results
result.predicted_stim_seconds = predicted_stim_seconds;
result.predicted_n_arousals   = sum(predicted_arousal);
if predicted_stim_seconds > 0
    result.predicted_arousal_rate = result.predicted_n_arousals / (predicted_stim_seconds / 60);
else
    result.predicted_arousal_rate = NaN;
end

% did any real stim onset depend on data the log doesn't contain?
% (a logging gap > 5 min inside the 60 samples before the watch started stimulating)
real_onsets = find(mit_stim_on & [false; ~mit_stim_on(1:end-1)]);
cg = cumsum(dt_s > 300);                                         % running count of long gaps
result.onset_after_long_gap = any(cg(real_onsets) - cg(max(real_onsets - stable_n_samples, 1)) > 0);

%% per-sample detail (used by diagnostic scripts; the batch ignores it)
if size(mit_parts, 2) >= 13, counter = str2double(mit_parts(:, 13)); else, counter = nan(n_mit, 1); end
% buffer_n is now computed in setup above)
result.detail = struct('mit_ts', double(mit_ts), 'tick_no', tick_no, 'inst', mit_inst_N3, 'avg', mit_avg_N3, ...
    'real_on', mit_stim_on, 'sim_on', predicted_stim, 'sim_arousal', predicted_arousal, ...
    'fail_reason', fail_reason, 'recent_motion', recent_motion, 'counter', counter, 'buffer_n', buffer_n, ...
    'arousal_ts', double(arousal_events.timestamp));