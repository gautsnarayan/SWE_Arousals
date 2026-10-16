%% main_effects_v3.m
% next-day behavior after entrainment vs sham nights (pre-registered), plus the
% exploratory comparisons audio vs sham, audio vs haptic and haptic vs sham.
% one lme per measure: DV ~ condition + (1|userID)
%
% inputs  : full_batch_h8_results.mat (pipeline batch) and the behavior csv
% outputs (in out_dir):
%   merged_all.csv            one row per night: condition, arousal data, behavior
%   main_effects_results.csv  n, F, df, p and q for every comparison and measure
%   main_effects_output.txt   everything printed to the command window
%   h1_h7_bars*.png           bar plots for each comparison


clear;

%% ---------------- Parameters ----------------
base          = 'C:\Users\sleep\Documents\Fluid_Interfaces';                 % <-- change this to your own folder
batch_dir     = fullfile(base, 'full_batch_h8_100626_new');                 % final batch
behavior_file = fullfile(base, 'year 2 study data (3).csv');
out_dir       = fullfile(base, 'main_effects');
rate_var      = "pred_rate";      % "real_rate" (watch-logged) | "pred_rate" (pipeline); identical after validation
dvs = {                            % <-- confirm against the pre-registered DV list
    'TMT_time',                       'Trail Making Test - time'
    'TMT_clicks',                     'Trail Making Test - clicks'
    'Total_LSEQ',                     'LSEQ total'
    'LSEQ_awake_following_sleep',     'LSEQ awake following sleep'
    'BRUMS_TMD',                      'BRUMS total mood disturbance'
    'BRUMS_fatigue',                  'BRUMS fatigue'
    'hr_min',                         'HR minimum'
    'sleep_time',                     'Total Sleep Time (min)'

    };

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% save everything printed below to a text file (start fresh; diary adds to the end of an old file)
log_file = fullfile(out_dir, 'main_effects_output.txt');
if isfile(log_file), delete(log_file); end
diary(log_file);

%% ---------------- Step 1: entrainment and sham nights from the final batch ----------------
S = load(fullfile(batch_dir, 'full_batch_h8_results.mat'));  
R = S.R;

% entrainment nights count only if the watch really stimulated;
% sham nights where the watch fired were already removed by the audit
fired = R.real_stim > 0;
is_entrain = ismember(R.condition, ["stimtype1", "stimtype2"]);
is_sham    = R.condition == "stimtype3";
E = R(R.error == "" & ((is_entrain & fired) | is_sham), :);
E.userID     = string(E.userID);
E.night_date = string(E.night_date);                      % derived night date (evening the night started)
E.arousal_rate = E.(rate_var);
E.stim_type  = repmat("audio", height(E), 1);
E.stim_type(E.condition == "stimtype2") = "haptic";
E.stim_type(E.condition == "stimtype3") = "sham";
% % fprintf('Entrainment nights with an arousal rate (%s): %d (audio %d, haptic %d), %d participants\n', ...
% %     rate_var, height(E), sum(E.stim_type == "audio"), sum(E.stim_type == "haptic"), numel(unique(E.userID)));

%% ---------------- Step 2: behavior table ----------------
beh = readtable(behavior_file, 'VariableNamingRule', 'modify');
beh.userID = string(beh.UserID);
if isdatetime(beh.Start_Date)
    beh_dt = beh.Start_Date;
else
    beh_dt = datetime(string(beh.Start_Date), 'InputFormat', 'MM/dd/yyyy');
end
fprintf('Behavior rows: %d | participants: %d\n', height(beh), numel(unique(beh.userID)));

% how many 0 values each measure has (a 0 can mean "not done" rather than a real score)
for k = 1:size(dvs, 1)
    if ismember(dvs{k, 1}, beh.Properties.VariableNames)
        fprintf('  %-32s %4d rows equal to 0\n', dvs{k, 1}, sum(beh.(dvs{k, 1}) == 0));
    end
end

% a tmt time or click count of 0 means the test was not done, so it is missing
tmt_cols = {'TMT_time', 'TMT_clicks'};
for k = 1:numel(tmt_cols)
    if ismember(tmt_cols{k}, beh.Properties.VariableNames)
        is_zero = beh.(tmt_cols{k}) == 0;
        fprintf('%s: %d rows with 0 set to missing\n', tmt_cols{k}, sum(is_zero));
        beh.(tmt_cols{k})(is_zero) = NaN;
    end
end

%% ---------------- Step 3: date alignment chosen empirically ----------------
offsets = [-1 0 1];                                       % behavior date shifted by this many days
ar_keys = E.userID + "_" + E.night_date;
n_match = zeros(size(offsets));
for o = 1:numel(offsets)
    bk = beh.userID + "_" + string(beh_dt + days(offsets(o)), 'yyyy-MM-dd');
    n_match(o) = numel(intersect(ar_keys, bk));
    fprintf('  behavior date %+d day: %4d matching userID+night keys\n', offsets(o), n_match(o));
end
[~, best] = max(n_match);
offset = offsets(best);
fprintf('Using behavior date %+d day (%d matches)\n', offset, n_match(best));
beh.night_date = string(beh_dt + days(offset), 'yyyy-MM-dd');
beh.key = beh.userID + "_" + beh.night_date;
E.key   = ar_keys;

% keys must be unique on both sides, or the merge pairs the wrong nights
[ub, ~, ib] = unique(beh.key);  dup_b = ub(accumarray(ib, 1) > 1);
[ue, ~, ie] = unique(E.key);    dup_e = ue(accumarray(ie, 1) > 1);
fprintf('Duplicate keys: behavior %d | arousal %d (rows with a duplicated key are dropped)\n', numel(dup_b), numel(dup_e));
beh = beh(~ismember(beh.key, dup_b), :);
E   = E(~ismember(E.key, dup_e), :);

%% ---------------- Step 4: merge ----------------
keep_e = {'key', 'userID', 'night_date', 'condition', 'stim_type', 'arousal_rate', 'real_rate', 'pred_rate', ...
          'real_stim', 'pred_stim', 'real_ar', 'pred_ar'};
dv_names = dvs(:, 1)';
missing_dvs = dv_names(~ismember(dv_names, beh.Properties.VariableNames));
if ~isempty(missing_dvs)
    fprintf('WARNING - DV columns not found in the behavior table: %s\n', strjoin(missing_dvs, ', '));
    dvs = dvs(ismember(dvs(:, 1), beh.Properties.VariableNames), :);
end
keep_b = [{'key'}, dvs(:, 1)'];
if ismember('StimNight', beh.Properties.VariableNames), keep_b = [keep_b, {'StimNight'}]; end
M = innerjoin(E(:, keep_e), beh(:, keep_b), 'Keys', 'key');
% fprintf('\nMerged entrainment nights: %d (audio %d, haptic %d), %d participants\n', height(M), ...
%     sum(M.stim_type == "audio"), sum(M.stim_type == "haptic"), numel(unique(M.userID)));

% how often the behavior file's condition label agrees with the audited condition
if ismember('StimNight', M.Properties.VariableNames)
    lab = lower(string(M.StimNight));
    fprintf('Behavior StimNight label agrees with audited condition on %d of %d nights (audited condition is used)\n', ...
        sum(lab == M.condition), height(M));
end

writetable(M, fullfile(out_dir, 'merged_all.csv'));

%% Step 5: main effects

%LSEQ will increase after nights with entrainment vs sham


% condition comes from the audit (stim_type), not the behavior file's label
M.entrainment = repmat("entrainment", height(M), 1);
M.entrainment(M.stim_type == "sham") = "non-entrainment";
M.entrainment = categorical(M.entrainment);

% Total LSEQ
lme_h1 = fitlme(M,'Total_LSEQ ~ entrainment + (1|userID)');
h1 = anova(lme_h1);
disp(h1)

% LSEQ awake following asleep
lme_h2 = fitlme(M,'LSEQ_awake_following_sleep ~ entrainment + (1|userID)');
h2 = anova(lme_h2);
disp(h2)

% BRUMS_fatigue
lme_h3 = fitlme(M,'BRUMS_fatigue ~ entrainment + (1|userID)');
h3 = anova(lme_h3);
disp(h3)

% BRUMS_TMD
lme_h4 = fitlme(M,'BRUMS_TMD ~ entrainment + (1|userID)');
h4 = anova(lme_h4);
disp(h4)

% TMT_time
lme_h6 = fitlme(M,'TMT_time ~ entrainment + (1|userID)');
h6 = anova(lme_h6);
disp(h6)

% TMT_clicks
lme_h7 = fitlme(M,'TMT_clicks ~ entrainment + (1|userID)');
h7 = anova(lme_h7);
disp(h7)

%HR min
lme_h9 = fitlme(M,'hr_min ~ entrainment + (1|userID)');
h9 = anova(lme_h9);
disp(h9)

%TST
lme_h10 = fitlme(M,'sleep_time ~ entrainment + (1|userID)');
h10 = anova(lme_h10);
disp(h10)


%% Step 6: bar graphs of each measure, sham vs entrainment
% a participant can have more than one night in a condition, so first average
% each participant's nights within a condition, then take the mean and sem
% across participants (this way no one counts twice in the error bars)

dvs     = {'Total_LSEQ', 'LSEQ_awake_following_sleep', 'BRUMS_fatigue', 'BRUMS_TMD', 'TMT_time', 'TMT_clicks', 'hr_min', 'sleep_time'};
labels  = {'Total LSEQ', 'LSEQ awake following sleep', 'BRUMS fatigue', 'BRUMS TMD', 'TMT time (s)', 'TMT clicks', 'HR min (beats)', 'TST (hrs)'};
anovas  = {h1, h2, h3, h4, h6, h7, h9, h10};                % anova tables from step 5, same order as dvs
conds   = ["non-entrainment", "entrainment"];      % left bar, right bar
bar_col = [0.60 0.60 0.60;                         % sham: grey
           0.40 0.20 0.60];                        % entrainment: purple

users = unique(M.userID);

figure('Color', 'w', 'Position', [100 100 1200 650]);
for d = 1:numel(dvs)
    y = M.(dvs{d});

    % one mean per participant per condition
    subj_mean = nan(numel(users), 2);
    for s = 1:numel(users)
        for c = 1:2
            rows = M.userID == users(s) & M.entrainment == conds(c) & ~isnan(y);
            if any(rows)
                subj_mean(s, c) = mean(y(rows));
            end
        end
    end

    % mean and sem across participants (participants missing a condition are skipped for that bar)
    group_mean = nan(1, 2);
    group_sem  = nan(1, 2);
    for c = 1:2
        vals = subj_mean(~isnan(subj_mean(:, c)), c);
        group_mean(c) = mean(vals);
        group_sem(c)  = std(vals) / sqrt(length(vals));
    end

    % draw the bars
    subplot(2, 4, d); hold on;
    for c = 1:2
        bar(c, group_mean(c), 0.6, 'FaceColor', bar_col(c, :), 'EdgeColor', 'k', 'LineWidth', 1);
    end
    errorbar(1:2, group_mean, group_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Sham', 'Entrainment'}, 'FontSize', 11, 'Box', 'off');
    xlim([0.4 2.6]);
    ylabel(labels{d});

    % f and p for the entrainment term from the lme
    tbl = anovas{d};
    row = strcmp(tbl.Term, 'entrainment');
    title({labels{d}, sprintf('F(%d,%d) = %.2f, p = %.4f', tbl.DF1(row), tbl.DF2(row), tbl.FStat(row), tbl.pValue(row))}, ...
        'FontWeight', 'normal');
end
sgtitle('Next-day behavior after entrainment vs sham nights (mean \pm SEM across participants)');

exportgraphics(gcf, fullfile(out_dir, 'h1_h7_bars.png'), 'Resolution', 300);

%% comparison 2: sham vs auditory

C2 = M(M.stim_type ~= 'haptic', :);

C2.entrainment = categorical(C2.entrainment);

% Total LSEQ
lme_h1_c2 = fitlme(C2,'Total_LSEQ ~ entrainment + (1|userID)');
h1_c2 = anova(lme_h1_c2);
disp(h1_c2)

% LSEQ awake following asleep
lme_h2_c2 = fitlme(C2,'LSEQ_awake_following_sleep ~ entrainment + (1|userID)');
h2_c2 = anova(lme_h2_c2);
disp(h2_c2)

% BRUMS_fatigue
lme_h3_c2 = fitlme(C2,'BRUMS_fatigue ~ entrainment + (1|userID)');
h3_c2 = anova(lme_h3_c2);
disp(h3_c2)

% BRUMS_TMD
lme_h4_c2 = fitlme(C2,'BRUMS_TMD ~ entrainment + (1|userID)');
h4_c2 = anova(lme_h4_c2);
disp(h4_c2)

% TMT_time
lme_h6_c2 = fitlme(C2,'TMT_time ~ entrainment + (1|userID)');
h6_c2 = anova(lme_h6_c2);
disp(h6_c2)

% TMT_clicks
lme_h7_c2 = fitlme(C2,'TMT_clicks ~ entrainment + (1|userID)');
h7_c2 = anova(lme_h7_c2);
disp(h7_c2)

%HR min
lme_h9_c2 = fitlme(C2,'hr_min ~ entrainment + (1|userID)');
h9_c2 = anova(lme_h9_c2);
disp(h9_c2)

%TST
lme_h10_c2 = fitlme(C2,'sleep_time ~ entrainment + (1|userID)');
h10_c2 = anova(lme_h10_c2);
disp(h10_c2)

%%  plotting comparison 2

dvs     = {'Total_LSEQ', 'LSEQ_awake_following_sleep', 'BRUMS_fatigue', 'BRUMS_TMD', 'TMT_time', 'TMT_clicks', 'hr_min', 'sleep_time'};
labels  = {'Total LSEQ', 'LSEQ awake following sleep', 'BRUMS fatigue', 'BRUMS TMD', 'TMT time (s)', 'TMT clicks', 'HR min (beats)', 'TST (hrs)'};
anovas_c2  = {h1_c2, h2_c2, h3_c2, h4_c2, h6_c2, h7_c2, h9_c2, h10_c2};                % anova tables from step 5, same order as dvs
conds   = ["sham", "audio"];      % left bar, right bar
bar_col = [0.60 0.60 0.60;                         % sham: grey
           0 0 1];                        % entrainment: blue

users = unique(M.userID);

figure('Color', 'w', 'Position', [100 100 1200 650]);
for d = 1:numel(dvs)
    y = C2.(dvs{d});

    % one mean per participant per condition
    subj_mean = nan(numel(users), 2);
    for s = 1:numel(users)
        for c = 1:2
            rows = C2.userID == users(s) & C2.stim_type == conds(c) & ~isnan(y);
            if any(rows)
                subj_mean(s, c) = mean(y(rows));
            end
        end
    end

    % mean and sem across participants (participants missing a condition are skipped for that bar)
    group_mean = nan(1, 2);
    group_sem  = nan(1, 2);
    for c = 1:2
        vals = subj_mean(~isnan(subj_mean(:, c)), c);
        group_mean(c) = mean(vals);
        group_sem(c)  = std(vals) / sqrt(length(vals));
    end

    % draw the bars
    subplot(2, 4, d); hold on;
    for c = 1:2
        bar(c, group_mean(c), 0.6, 'FaceColor', bar_col(c, :), 'EdgeColor', 'k', 'LineWidth', 1);
    end
    errorbar(1:2, group_mean, group_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Sham', 'Audio'}, 'FontSize', 11, 'Box', 'off');
    xlim([0.4 2.6]);
    ylabel(labels{d});

    % f and p for the entrainment term from the lme
    tbl = anovas_c2{d};
    row = strcmp(tbl.Term, 'entrainment');
    title({labels{d}, sprintf('F(%d,%d) = %.2f, p = %.4f', tbl.DF1(row), tbl.DF2(row), tbl.FStat(row), tbl.pValue(row))}, ...
        'FontWeight', 'normal');
end
sgtitle('Next-day behavior after audio vs sham nights (mean \pm SEM across participants)');

exportgraphics(gcf, fullfile(out_dir, 'h1_h7_bars_C2.png'), 'Resolution', 300);
%% comparison 3: haptic vs auditory

C3 = M(M.stim_type ~= 'sham', :);

C3.stim_type = categorical(C3.stim_type);

% Total LSEQ
lme_h1_c3 = fitlme(C3,'Total_LSEQ ~ stim_type + (1|userID)');
h1_c3 = anova(lme_h1_c3);
disp(h1_c3)

% LSEQ awake following asleep
lme_h2_c3 = fitlme(C3,'LSEQ_awake_following_sleep ~ stim_type + (1|userID)');
h2_c3 = anova(lme_h2_c3);
disp(h2_c3)

% BRUMS_fatigue
lme_h3_c3 = fitlme(C3,'BRUMS_fatigue ~ stim_type + (1|userID)');
h3_c3 = anova(lme_h3_c3);
disp(h3_c3)

% BRUMS_TMD
lme_h4_c3 = fitlme(C3,'BRUMS_TMD ~ stim_type + (1|userID)');
h4_c3 = anova(lme_h4_c3);
disp(h4_c3)

% TMT_time
lme_h6_c3 = fitlme(C3,'TMT_time ~ stim_type + (1|userID)');
h6_c3 = anova(lme_h6_c3);
disp(h6_c3)

% TMT_clicks
lme_h7_c3 = fitlme(C3,'TMT_clicks ~ stim_type + (1|userID)');
h7_c3 = anova(lme_h7_c3);
disp(h7_c3)

%HR min
lme_h9_c3 = fitlme(C3,'hr_min ~ stim_type + (1|userID)');
h9_c3 = anova(lme_h9_c3);
disp(h9_c3)

%TST
lme_h10_c3 = fitlme(C3,'sleep_time ~ stim_type + (1|userID)');
h10_c3 = anova(lme_h10_c3);
disp(h10_c3)

%%  plotting comparison 3: sound vs haptic

dvs     = {'Total_LSEQ', 'LSEQ_awake_following_sleep', 'BRUMS_fatigue', 'BRUMS_TMD', 'TMT_time', 'TMT_clicks', 'hr_min', 'sleep_time'};
labels  = {'Total LSEQ', 'LSEQ awake following sleep', 'BRUMS fatigue', 'BRUMS TMD', 'TMT time (s)', 'TMT clicks', 'HR min (beats)', 'TST (hrs)'};
anovas_c3  = {h1_c3, h2_c3, h3_c3, h4_c3, h6_c3, h7_c3, h9_c3, h10_c3};                % anova tables from step 5, same order as dvs
conds   = ["audio", "haptic"];      % left bar, right bar
bar_col = [0 0 1;                         % audio: blue
           1 0 0];                        % haptic: red

users = unique(M.userID);

figure('Color', 'w', 'Position', [100 100 1200 650]);
for d = 1:numel(dvs)
    y = C3.(dvs{d});

    % one mean per participant per condition
    subj_mean = nan(numel(users), 2);
    for s = 1:numel(users)
        for c = 1:2
            rows = C3.userID == users(s) & C3.stim_type == conds(c) & ~isnan(y);
            if any(rows)
                subj_mean(s, c) = mean(y(rows));
            end
        end
    end

    % mean and sem across participants (participants missing a condition are skipped for that bar)
    group_mean = nan(1, 2);
    group_sem  = nan(1, 2);
    for c = 1:2
        vals = subj_mean(~isnan(subj_mean(:, c)), c);
        group_mean(c) = mean(vals);
        group_sem(c)  = std(vals) / sqrt(length(vals));
    end

    % draw the bars
    subplot(2, 4, d); hold on;
    for c = 1:2
        bar(c, group_mean(c), 0.6, 'FaceColor', bar_col(c, :), 'EdgeColor', 'k', 'LineWidth', 1);
    end
    errorbar(1:2, group_mean, group_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Audio', 'Haptic'}, 'FontSize', 11, 'Box', 'off');
    xlim([0.4 2.6]);
    ylabel(labels{d});

    % f and p for the entrainment term from the lme
    tbl = anovas_c3{d};
    row = strcmp(tbl.Term, 'stim_type');
    title({labels{d}, sprintf('F(%d,%d) = %.2f, p = %.4f', tbl.DF1(row), tbl.DF2(row), tbl.FStat(row), tbl.pValue(row))}, ...
        'FontWeight', 'normal');
end
sgtitle('Next-day behavior after audio vs haptic nights (mean \pm SEM across participants)');

exportgraphics(gcf, fullfile(out_dir, 'h1_h7_bars_C3.png'), 'Resolution', 300);

%% comparison 4: haptic vs sham


C4 = M(M.stim_type ~= 'audio', :);

C4.entrainment = categorical(C4.entrainment);
C4.stim_type = categorical(C4.stim_type);

% Total LSEQ
lme_h1_c4 = fitlme(C4,'Total_LSEQ ~ entrainment + (1|userID)');
h1_c4 = anova(lme_h1_c4);
disp(h1_c4)

% LSEQ awake following asleep
lme_h2_c4 = fitlme(C4,'LSEQ_awake_following_sleep ~ entrainment + (1|userID)');
h2_c4 = anova(lme_h2_c4);
disp(h2_c4)

% BRUMS_fatigue
lme_h3_c4 = fitlme(C4,'BRUMS_fatigue ~ entrainment + (1|userID)');
h3_c4 = anova(lme_h3_c4);
disp(h3_c4)

% BRUMS_TMD
lme_h4_c4 = fitlme(C4,'BRUMS_TMD ~ entrainment + (1|userID)');
h4_c4 = anova(lme_h4_c4);
disp(h4_c4)

% TMT_time
lme_h6_c4 = fitlme(C4,'TMT_time ~ entrainment + (1|userID)');
h6_c4 = anova(lme_h6_c4);
disp(h6_c4)

% TMT_clicks
lme_h7_c4 = fitlme(C4,'TMT_clicks ~ entrainment + (1|userID)');
h7_c4 = anova(lme_h7_c4);
disp(h7_c4)

%HR min
lme_h9_c4 = fitlme(C4,'hr_min ~ entrainment + (1|userID)');
h9_c4 = anova(lme_h9_c4);
disp(h9_c4)

%TST
lme_h10_c4 = fitlme(C4,'sleep_time ~ entrainment + (1|userID)');
h10_c4 = anova(lme_h10_c4);
disp(h10_c4)

%%  plotting comparison 4


dvs     = {'Total_LSEQ', 'LSEQ_awake_following_sleep', 'BRUMS_fatigue', 'BRUMS_TMD', 'TMT_time', 'TMT_clicks', 'hr_min', 'sleep_time'};
labels  = {'Total LSEQ', 'LSEQ awake following sleep', 'BRUMS fatigue', 'BRUMS TMD', 'TMT time (s)', 'TMT clicks', 'HR min (beats)', 'TST (hrs)'};
anovas_c4  = {h1_c4, h2_c4, h3_c4, h4_c4, h6_c4, h7_c4, h9_c4, h10_c4};                % anova tables from step 5, same order as dvs
conds   = ["sham", "haptic"];      % left bar, right bar
bar_col = [0.60 0.60 0.60;                         % sham: grey
           1 0 0];                        % entrainment: red

users = unique(M.userID);

figure('Color', 'w', 'Position', [100 100 1200 650]);
for d = 1:numel(dvs)
    y = C4.(dvs{d});

    % one mean per participant per condition
    subj_mean = nan(numel(users), 2);
    for s = 1:numel(users)
        for c = 1:2
            rows = C4.userID == users(s) & C4.stim_type == conds(c) & ~isnan(y);
            if any(rows)
                subj_mean(s, c) = mean(y(rows));
            end
        end
    end

    % mean and sem across participants (participants missing a condition are skipped for that bar)
    group_mean = nan(1, 2);
    group_sem  = nan(1, 2);
    for c = 1:2
        vals = subj_mean(~isnan(subj_mean(:, c)), c);
        group_mean(c) = mean(vals);
        group_sem(c)  = std(vals) / sqrt(length(vals));
    end

    % draw the bars
    subplot(2, 4, d); hold on;
    for c = 1:2
        bar(c, group_mean(c), 0.6, 'FaceColor', bar_col(c, :), 'EdgeColor', 'k', 'LineWidth', 1);
    end
    errorbar(1:2, group_mean, group_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Sham', 'Haptic'}, 'FontSize', 11, 'Box', 'off');
    xlim([0.4 2.6]);
    ylabel(labels{d});

    % f and p for the entrainment term from the lme
    tbl = anovas_c4{d};
    row = strcmp(tbl.Term, 'entrainment');
    title({labels{d}, sprintf('F(%d,%d) = %.2f, p = %.4f', tbl.DF1(row), tbl.DF2(row), tbl.FStat(row), tbl.pValue(row))}, ...
        'FontWeight', 'normal');
end
sgtitle('Next-day behavior after haptic vs sham nights (mean \pm SEM across participants)');

exportgraphics(gcf, fullfile(out_dir, 'h1_h7_bars_C4.png'), 'Resolution', 300);

%% ---------------- FDR: benjamini-hochberg across the measures, per comparison ----------------
% each comparison is its own family of tests (one p value per measure)
all_anovas = {anovas, anovas_c2, anovas_c3, anovas_c4};
all_data   = {M, C2, C3, C4};                                              % nights used in each comparison
comp_names = {'entrainment vs sham', 'audio vs sham', 'audio vs haptic', 'haptic vs sham'};
term_names = {'entrainment', 'entrainment', 'stim_type', 'entrainment'};   % tested term in each comparison

% results table: one row per comparison x measure, filled in the loop below
res_comp = strings(0, 1);  res_dv = strings(0, 1);  res_n = [];
res_F = [];  res_df1 = [];  res_df2 = [];  res_p = [];  res_q = [];

for c = 1:numel(all_anovas)
    n_dv = numel(dvs);

    % p value of the tested term for each measure
    p = nan(n_dv, 1);
    for d = 1:n_dv
        tbl = all_anovas{c}{d};
        p(d) = tbl.pValue(strcmp(tbl.Term, term_names{c}));
    end

    % sort p from smallest to largest; q = p * m / rank
    [p_sorted, order] = sort(p);
    q_sorted = p_sorted .* n_dv ./ (1:n_dv)';

    % a q value can't be larger than the q of the next-larger p (work down from the largest)
    for i = n_dv-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    q_sorted = min(q_sorted, 1);

    % put the q values back in the original order
    q = zeros(n_dv, 1);
    q(order) = q_sorted;

    fprintf('\n=== %s: BH-FDR across %d measures ===\n', comp_names{c}, n_dv);
    for d = 1:n_dv
        fprintf('  %-32s p = %.4f   q = %.4f\n', labels{d}, p(d), q(d));

        % add this comparison x measure to the results table
        tbl = all_anovas{c}{d};
        row = strcmp(tbl.Term, term_names{c});
        res_comp(end+1, 1) = string(comp_names{c});
        res_dv(end+1, 1)   = string(dvs{d});
        res_n(end+1, 1)    = sum(~isnan(all_data{c}.(dvs{d})));             % nights in this model
        res_F(end+1, 1)    = tbl.FStat(row);
        res_df1(end+1, 1)  = tbl.DF1(row);
        res_df2(end+1, 1)  = tbl.DF2(row);
        res_p(end+1, 1)    = p(d);
        res_q(end+1, 1)    = q(d);
    end
end

results = table(res_comp, res_dv, res_n, res_F, res_df1, res_df2, res_p, res_q, ...
    'VariableNames', {'comparison', 'measure', 'n_nights', 'F', 'df1', 'df2', 'p', 'q'});
writetable(results, fullfile(out_dir, 'main_effects_results.csv'));
fprintf('\nresults saved to %s\n', fullfile(out_dir, 'main_effects_results.csv'));
diary off;