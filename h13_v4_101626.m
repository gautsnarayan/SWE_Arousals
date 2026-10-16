%% h13_all_v6.m
% everything for H13 
%
% ANALYSIS 1 - H13 
%   on ENTRAINMENT nights only, does the arousal rate, or its interaction with
%   entrainment type (audio vs haptic), predict next-morning outcomes?
%   arousal rate = watch-logged (real) rate
%
% ANALYSIS 2 - entrainment vs sham (exploratory, not pre-registered)
%   same question, but comparing entrainment nights (audio + haptic) with sham nights.
%   arousal rate = predicted rate from the pipeline. on sham nights the watch did not
%   stimulate, so the predicted rate is the rate the watch would have logged: how often
%   n3 dropped or the person moved during the stretches the watch would have stimulated.
%   the interaction asks: is the link between arousal rate and next-morning behavior
%   different when stimulation was actually delivered?
%
% both analyses use the same model, one per outcome measure:
%   DV ~ group * arousal_rate + (1|userID)
%   - group: audio vs haptic (analysis 1) or entrainment vs sham (analysis 2),
%     effects-coded, so the arousal_rate term is the slope averaged over the two groups
%   - REML fit (as lmer in R); Satterthwaite df for the p values (as pre-registered)
%   - tested terms: group:arousal_rate (main focus) and arousal_rate
%   - Benjamini-Hochberg FDR across the 7 measures, separately for each term
%
% merge: each night is matched to its own next-morning behavior record on userID + date.
% the date alignment (same day / +-1 day) is chosen from the data (most matches).
% TMT time / clicks of 0 mean the test was not done and are set to missing.
%
% how to run
%   1. set base below to your own folder
%   2. run with min_stim_s = 0   -> main result (all nights)
%   3. run with min_stim_s = 300 -> sensitivity check (nights with < 5 min of stim
%      dropped, because the arousal rate is unstable when stim is very short)
%


clear;

%% ---------------- Parameters ----------------
base          = 'C:\Users\sleep\Documents\Fluid_Interfaces';                 % <-- change this to your own folder
batch_dir     = fullfile(base, 'full_batch_h8_100626_new');                 % final batch
behavior_file = fullfile(base, 'year 2 study data (3).csv');
out_dir       = fullfile(base, 'h13');
min_stim_s    = 300;                % 0 = all nights (main result) | 300 = sensitivity: drop nights with < 5 min stim
dvs = {                            % pre-registered H13 outcome measures
    'TMT_time',                       'Trail Making Test - time'
    'TMT_clicks',                     'Trail Making Test - clicks'
    'Total_LSEQ',                     'LSEQ total'
    'LSEQ_awake_following_sleep',     'LSEQ awake following sleep'
    'LSEQ_behavior_following_waking', 'LSEQ behaviour following waking'
    'BRUMS_TMD',                      'BRUMS total mood disturbance'
    'BRUMS_fatigue',                  'BRUMS fatigue'
    };

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% file name suffix: none for the main run, e.g. "_minstim300" for the sensitivity run
if min_stim_s > 0
    tag = sprintf('_minstim%d', min_stim_s);
else
    tag = '';
end

% save everything printed below to a text file (start fresh; diary adds to the end of an old file)
log_file = fullfile(out_dir, ['h13_output' tag '.txt']);
if isfile(log_file), delete(log_file); end
diary(log_file);

if min_stim_s > 0
    fprintf('*** SENSITIVITY RUN: nights with < %d s of stim excluded ***\n\n', min_stim_s);
else
    fprintf('*** MAIN RUN: all nights ***\n\n');
end

%% ---------------- Step 1: nights from the final batch ----------------
S = load(fullfile(batch_dir, 'full_batch_h8_results.mat'));  R = S.R;
R.userID     = string(R.userID);
R.night_date = string(R.night_date);                      % derived night date (evening the night started)

ok      = R.error == "";
is_entr = ok & ismember(R.condition, ["stimtype1", "stimtype2"]);
is_sham = ok & R.condition == "stimtype3";

% sample flow for analysis 1: one row per step = {what is left, nights, participants}
flow = cell(0, 3);
flow(end+1, :) = {'audited nights (all conditions)', sum(ok), numel(unique(R.userID(ok)))};
flow(end+1, :) = {'entrainment nights only (H13 excludes sham)', sum(is_entr), numel(unique(R.userID(is_entr)))};

% analysis 1: entrainment nights with a watch-logged arousal rate
E1 = R(is_entr & ~isnan(R.real_rate), :);
E1.arousal_rate = E1.real_rate;
E1.stim_s       = E1.real_stim;                           % stim seconds behind that rate
E1.group        = repmat("audio", height(E1), 1);
E1.group(E1.condition == "stimtype2") = "haptic";
flow(end+1, :) = {'with an arousal rate (some stim)', height(E1), numel(unique(E1.userID))};
fprintf('Analysis 1: %d entrainment nights (audio %d, haptic %d), %d participants\n', ...
    height(E1), sum(E1.group == "audio"), sum(E1.group == "haptic"), numel(unique(E1.userID)));

% analysis 2: entrainment and sham nights with a predicted arousal rate
% (a sham night has no predicted rate if the simulated watch would never have started stim)
E2 = R((is_entr | is_sham) & ~isnan(R.pred_rate), :);
E2.arousal_rate = E2.pred_rate;
E2.stim_s       = E2.pred_stim;
E2.group        = repmat("entrainment", height(E2), 1);
E2.group(E2.condition == "stimtype3") = "sham";
fprintf('Analysis 2: %d entrainment and %d sham nights with a predicted rate, %d participants\n', ...
    sum(E2.group == "entrainment"), sum(E2.group == "sham"), numel(unique(E2.userID)));
fprintf('            (%d sham nights left out: the watch would never have started stim)\n', ...
    sum(is_sham) - sum(E2.group == "sham"));

%% ---------------- Step 2: behavior table ----------------
beh = readtable(behavior_file, 'VariableNamingRule', 'modify');
beh.userID = string(beh.UserID);
if isdatetime(beh.Start_Date)
    beh_dt = beh.Start_Date;
else
    beh_dt = datetime(string(beh.Start_Date), 'InputFormat', 'MM/dd/yyyy');
end
fprintf('\nBehavior rows: %d | participants: %d\n', height(beh), numel(unique(beh.userID)));

% a tmt time or click count of 0 means the test was not done, so it is missing
tmt_cols = {'TMT_time', 'TMT_clicks'};
for k = 1:numel(tmt_cols)
    if ismember(tmt_cols{k}, beh.Properties.VariableNames)
        is_zero = beh.(tmt_cols{k}) == 0;
        fprintf('%s: %d rows with 0 set to missing\n', tmt_cols{k}, sum(is_zero));
        beh.(tmt_cols{k})(is_zero) = NaN;
    end
end

% keep only the outcome measures that exist in the behavior table
dv_names = dvs(:, 1)';
missing_dvs = dv_names(~ismember(dv_names, beh.Properties.VariableNames));
if ~isempty(missing_dvs)
    fprintf('WARNING - DV columns not found in the behavior table: %s\n', strjoin(missing_dvs, ', '));
    dvs = dvs(ismember(dvs(:, 1), beh.Properties.VariableNames), :);
end
n_dv = size(dvs, 1);

%% ---------------- Step 3: date alignment chosen from the data ----------------
E1.key = E1.userID + "_" + E1.night_date;
E2.key = E2.userID + "_" + E2.night_date;

% try shifting the behavior date by -1, 0 and +1 day; use the shift with the most matches
offsets = [-1 0 1];
n_match = zeros(size(offsets));
for o = 1:numel(offsets)
    bk = beh.userID + "_" + string(beh_dt + days(offsets(o)), 'yyyy-MM-dd');
    n_match(o) = numel(intersect(E1.key, bk));
    fprintf('  behavior date %+d day: %4d matching userID+night keys\n', offsets(o), n_match(o));
end
[~, best] = max(n_match);
offset = offsets(best);
fprintf('Using behavior date %+d day (%d matches)\n', offset, n_match(best));
beh.night_date = string(beh_dt + days(offset), 'yyyy-MM-dd');
beh.key = beh.userID + "_" + beh.night_date;

% a key that appears more than once can't be matched safely, so those rows are dropped
counts = groupcounts(beh, 'key');  dup_b = counts.key(counts.GroupCount > 1);
counts = groupcounts(E1, 'key');   dup_1 = counts.key(counts.GroupCount > 1);
counts = groupcounts(E2, 'key');   dup_2 = counts.key(counts.GroupCount > 1);
fprintf('Duplicate keys: behavior %d | analysis 1 nights %d | analysis 2 nights %d (dropped)\n', ...
    numel(dup_b), numel(dup_1), numel(dup_2));

% sample flow for analysis 1: why entrainment nights don't reach the merge
s = E1(~ismember(E1.key, dup_1), :);
flow(end+1, :) = {'night not duplicated', height(s), numel(unique(s.userID))};
s = s(ismember(s.key, beh.key), :);
flow(end+1, :) = {'has a next-morning behavior record', height(s), numel(unique(s.userID))};
s = s(~ismember(s.key, dup_b), :);
flow(end+1, :) = {'behavior record not duplicated (= merged)', height(s), numel(unique(s.userID))};

beh = beh(~ismember(beh.key, dup_b), :);
E1  = E1(~ismember(E1.key, dup_1), :);
E2  = E2(~ismember(E2.key, dup_2), :);

%% ---------------- Step 4: merge nights with behavior ----------------
keep_e = {'key', 'userID', 'night_date', 'condition', 'group', 'arousal_rate', 'stim_s', ...
          'real_rate', 'pred_rate', 'real_stim', 'pred_stim', 'real_ar', 'pred_ar'};
keep_b = [{'key'}, dvs(:, 1)'];
if ismember('StimNight', beh.Properties.VariableNames), keep_b = [keep_b, {'StimNight'}]; end

M1 = innerjoin(E1(:, keep_e), beh(:, keep_b), 'Keys', 'key');
M2 = innerjoin(E2(:, keep_e), beh(:, keep_b), 'Keys', 'key');
fprintf('\nMerged, analysis 1: %d nights (audio %d, haptic %d), %d participants\n', height(M1), ...
    sum(M1.group == "audio"), sum(M1.group == "haptic"), numel(unique(M1.userID)));
fprintf('Merged, analysis 2: %d nights (entrainment %d, sham %d), %d participants\n', height(M2), ...
    sum(M2.group == "entrainment"), sum(M2.group == "sham"), numel(unique(M2.userID)));

% how often the behavior file's condition label agrees with the audited condition
if ismember('StimNight', M2.Properties.VariableNames)
    lab = lower(string(M2.StimNight));
    fprintf('Behavior StimNight label agrees with audited condition on %d of %d nights (audited condition is used)\n', ...
        sum(lab == M2.condition), height(M2));
end

% sensitivity rule: drop nights with little stim (the arousal rate is unstable when stim is short)
if min_stim_s > 0
    keep1 = M1.stim_s >= min_stim_s;
    keep2 = M2.stim_s >= min_stim_s;
    fprintf('SENSITIVITY: min stim %d s drops %d of %d nights in analysis 1 and %d of %d in analysis 2\n', ...
        min_stim_s, sum(~keep1), height(M1), sum(~keep2), height(M2));
    M1 = M1(keep1, :);
    M2 = M2(keep2, :);
    flow(end+1, :) = {sprintf('at least %d s of stim (sensitivity)', min_stim_s), height(M1), numel(unique(M1.userID))};
end

writetable(M1, fullfile(out_dir, ['h13_merged_entrainment' tag '.csv']));
writetable(M2, fullfile(out_dir, ['h13_merged_entrain_vs_sham' tag '.csv']));

%% ---------------- Step 5: one LME per measure, for each analysis ----------------
analysis_names = {'ANALYSIS 1 - H13: audio vs haptic, entrainment nights, watch-logged arousal rate', ...
                  'ANALYSIS 2 - exploratory: entrainment vs sham, predicted arousal rate'};
group_levels   = {["audio", "haptic"], ["entrainment", "sham"]};
data_sets      = {M1, M2};
results        = cell(1, 2);

for a = 1:2
    lv = group_levels{a};
    D_all = data_sets{a};
    D_all.group  = categorical(D_all.group, lv);
    D_all.userID = categorical(D_all.userID);
    data_sets{a} = D_all;                                 % keep this version for the plots

    T = table('Size', [n_dv, 14], ...
        'VariableTypes', [{'string', 'string', 'double', 'double'}, repmat({'double'}, 1, 10)], ...
        'VariableNames', {'dv', 'label', 'n', 'n_part', 'b_rate', 'F_rate', 'df2_rate', 'p_rate', ...
                          'b_int', 'F_int', 'df2_int', 'p_int', 'slope_group1', 'slope_group2'});

    fprintf('\n\n######## %s ########\n', analysis_names{a});
    for k = 1:n_dv
        dv = dvs{k, 1};
        D  = D_all(~isnan(D_all.(dv)) & ~isnan(D_all.arousal_rate), :);
        lme = fitlme(D, sprintf('%s ~ group * arousal_rate + (1|userID)', dv), ...
            'FitMethod', 'REML', 'DummyVarCoding', 'effects');
        an = anova(lme, 'DFMethod', 'satterthwaite');     % Satterthwaite F tests per term

        fprintf('\n==================== %s (%s) ====================\n', dvs{k, 2}, dv);
        fprintf('Model: %s\n', char(lme.Formula));
        fprintf('n = %d nights, %d participants | REML | effects-coded group (%s vs %s)\n', ...
            height(D), numel(unique(D.userID)), lv(1), lv(2));
        fprintf('\nFit statistics:\n');
        disp(lme.ModelCriterion);

        fprintf('Fixed effects (Satterthwaite df):\n');
        [~, ~, fe] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
        disp(fe);

        fprintf('Random effects:\n');
        [~, ~, re] = covarianceParameters(lme);
        disp(re{1});                                      % participant intercept SD
        fprintf('Residual SD: %.4f\n\n', sqrt(lme.MSE));

        fprintf('ANOVA (type III, Satterthwaite df):\n');
        disp(an);

        % the two tested terms in the anova table
        terms = string(an.Term);
        ir = terms == "arousal_rate";
        ii = contains(terms, "group") & contains(terms, "arousal_rate");

        % and their coefficients
        [beta, nm] = fixedEffects(lme);
        nm = string(nm.Name);
        is_int = contains(nm, "group") & contains(nm, "arousal_rate");
        b_rate = beta(nm == "arousal_rate");
        b_int  = beta(is_int);

        % slope in each group. with effects coding the group column is +1 for the
        % level in the coefficient name and -1 for the other level
        if contains(nm(is_int), lv(1))
            slope_1 = b_rate + b_int;
            slope_2 = b_rate - b_int;
        else
            slope_1 = b_rate - b_int;
            slope_2 = b_rate + b_int;
        end

        T(k, :) = {string(dv), string(dvs{k, 2}), height(D), numel(unique(D.userID)), ...
            b_rate, an.FStat(ir), an.DF2(ir), an.pValue(ir), ...
            b_int,  an.FStat(ii), an.DF2(ii), an.pValue(ii), slope_1, slope_2};
    end

    % say which group is which in the saved table
    T.group1 = repmat(lv(1), n_dv, 1);
    T.group2 = repmat(lv(2), n_dv, 1);
    results{a} = T;
end

%% ---------------- Step 6: Benjamini-Hochberg FDR across the measures, per term ----------------
% in each analysis, two families of tests: the interaction p values and the arousal_rate p values
for a = 1:2
    T = results{a};
    for f = 1:2
        if f == 1
            p = T.p_int;
        else
            p = T.p_rate;
        end
        m = numel(p);

        % sort p from smallest to largest; q = p * m / rank
        [p_sorted, order] = sort(p);
        q_sorted = p_sorted .* m ./ (1:m)';

        % a q value can't be larger than the q of the next-larger p (work down from the largest)
        for i = m-1:-1:1
            q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
        end
        q_sorted = min(q_sorted, 1);

        % put the q values back in the original measure order
        q = zeros(m, 1);
        q(order) = q_sorted;

        if f == 1
            T.q_int = q;
        else
            T.q_rate = q;
        end
    end
    results{a} = T;
end

%% ---------------- Step 7: report ----------------
% sample flow for analysis 1: what is left after each step, and what that step removed
T = results{1};
fprintf('\n\n=== SAMPLE FLOW (analysis 1): from audited nights to the H13 models ===\n');
fprintf('%-46s %7s %7s | %8s %8s\n', 'step (what is left)', 'nights', 'people', '-nights', '-people');
for i = 1:size(flow, 1)
    if i == 1
        fprintf('%-46s %7d %7d |\n', flow{i, 1}, flow{i, 2}, flow{i, 3});
    else
        fprintf('%-46s %7d %7d | %8d %8d\n', flow{i, 1}, flow{i, 2}, flow{i, 3}, ...
            flow{i - 1, 2} - flow{i, 2}, flow{i - 1, 3} - flow{i, 3});
    end
end
fprintf('then, per model, nights with that measure missing are left out:\n');
for k = 1:n_dv
    fprintf('  %-44s %7d %7d | %8d %8d\n', T.label(k), T.n(k), T.n_part(k), ...
        height(M1) - T.n(k), numel(unique(M1.userID)) - T.n_part(k));
end
fprintf('(-people counts participants who lost ALL their remaining nights at that step)\n');

% summary table for each analysis
for a = 1:2
    T  = results{a};
    lv = group_levels{a};
    fprintf('\n\n=== %s ===\n', analysis_names{a});
    fprintf('DV ~ group * arousal_rate + (1|userID), REML, Satterthwaite df\n');
    if min_stim_s > 0
        fprintf('SENSITIVITY: nights with < %d s of stim excluded\n', min_stim_s);
    end
    fprintf('BH-FDR across the %d measures, separately for the interaction and for arousal rate\n\n', n_dv);
    fprintf('%-32s %4s %4s | %-53s | %-32s\n', '', '', '', ...
        sprintf('INTERACTION: %s vs %s x arousal rate', lv(1), lv(2)), 'AROUSAL RATE (avg of both groups)');
    fprintf('%-32s %4s %4s | %11s %11s %6s %6s %7s %7s | %8s %6s %6s %7s %7s\n', 'DV', 'n', 'ppl', ...
        lv(1), lv(2), 'F', 'df', 'p', 'q', 'slope', 'F', 'df', 'p', 'q');
    for k = 1:n_dv
        fprintf('%-32s %4d %4d | %+11.3f %+11.3f %6.2f %6.1f %7.3f %7.3f | %+8.3f %6.2f %6.1f %7.3f %7.3f\n', ...
            T.label(k), T.n(k), T.n_part(k), ...
            T.slope_group1(k), T.slope_group2(k), T.F_int(k), T.df2_int(k), T.p_int(k), T.q_int(k), ...
            T.b_rate(k), T.F_rate(k), T.df2_rate(k), T.p_rate(k), T.q_rate(k));
    end
    fprintf('\n%s / %s = change in the DV per +1 arousal/min on %s / %s nights.\n', lv(1), lv(2), lv(1), lv(2));
    fprintf('interaction p = do those two slopes differ? F(1, df), Satterthwaite df.\n');
    fprintf('arousal rate  = the slope averaged over both groups (effects coding).\n');
    fprintf('q             = Benjamini-Hochberg FDR across the %d measures, separately for each term.\n', n_dv);
end

% save the results
writetable(results{1}, fullfile(out_dir, ['h13_results' tag '.csv']));
writetable(results{2}, fullfile(out_dir, ['h13_entrain_vs_sham_results' tag '.csv']));
T1 = results{1};  T2 = results{2};
save(fullfile(out_dir, ['h13_results' tag '.mat']), 'T1', 'T2', 'M1', 'M2', 'flow', 'offset', 'min_stim_s');

%% ---------------- Step 8: scatter plots, one figure per analysis ----------------
% each panel: one measure vs arousal rate, one point per night, coloured by group,
% with the lme p and q (interaction and arousal rate) and the spearman correlation
group_names = {{'Audio', 'Haptic'}, {'Entrainment', 'Sham'}};
group_cols  = {[0.10 0.65 0.30; 0.00 0.45 0.85], ...     % audio green, haptic blue
               [1 0 1; 0 0 0]};        % entrainment purple, sham black
x_labels    = {'Arousal rate (/min of stim, watch-logged)', 'Predicted arousal rate (/min of stim)'};
fig_titles  = {'H13: next-morning behaviour vs arousal rate on entrainment nights, audio vs haptic', ...
               'Exploratory: next-morning behaviour vs predicted arousal rate, entrainment vs sham'};
file_names  = {'h13_behaviour_vs_arousal', 'h13_entrain_vs_sham_vs_arousal'};

if min_stim_s > 0
    nights_text = sprintf('Nights shown: stim >= %d s', min_stim_s);
else
    nights_text = 'Nights shown: all nights';
end

for a = 1:2
    P   = data_sets{a};
    T   = results{a};
    lv  = group_levels{a};
    gn  = group_names{a};
    col = group_cols{a};
    in_g1 = P.group == lv(1);

    figure('Color', 'w', 'Position', [40 60 1600 820]);
    for k = 1:n_dv
        subplot(2, 4, k); hold on;
        y = P.(dvs{k, 1});
        x = P.arousal_rate;
        good = ~isnan(x) & ~isnan(y);

        % points
        scatter(x(good & in_g1),  y(good & in_g1),  18, col(1, :), 'filled', 'MarkerFaceAlpha', 0.45);
        scatter(x(good & ~in_g1), y(good & ~in_g1), 18, col(2, :), 'filled', 'MarkerFaceAlpha', 0.45);

        % least-squares line for each group (drawn over that group's range)
        for g = 1:2
            if g == 1
                m = good & in_g1;
            else
                m = good & ~in_g1;
            end
            if sum(m) > 2
                pf = polyfit(x(m), y(m), 1);
                xx = linspace(min(x(m)), max(x(m)), 100);
                plot(xx, polyval(pf, xx), '-', 'Color', col(g, :), 'LineWidth', 1.8);
            end
        end

        % descriptive spearman correlation: all nights, then each group
        [rho, p]    = corr(x(good), y(good), 'Type', 'Spearman');
        [rho_1, ~]  = corr(x(good & in_g1),  y(good & in_g1),  'Type', 'Spearman');
        [rho_2, ~]  = corr(x(good & ~in_g1), y(good & ~in_g1), 'Type', 'Spearman');
        rho_line = sprintf('\\rho = %.2f (p = %.3f) | %s %.2f, %s %.2f', rho, p, ...
            lower(gn{1}), rho_1, lower(gn{2}), rho_2);

        % lme p and q for this measure (rows of T are in the same order as dvs)
        lme_line = sprintf('Interaction p = %.3f (q = %.3f) | Rate p = %.3f (q = %.3f)', ...
            T.p_int(k), T.q_int(k), T.p_rate(k), T.q_rate(k));

        title({dvs{k, 2}, lme_line, rho_line}, 'FontSize', 8.5);
        xlabel(x_labels{a});
        ylabel(dvs{k, 2});
        box off;
    end

    % legend and key in the empty 8th panel
    subplot(2, 4, 8); axis off; hold on;
    scatter(nan, nan, 40, col(1, :), 'filled');
    scatter(nan, nan, 40, col(2, :), 'filled');
    legend(gn, 'Location', 'north', 'FontSize', 11); legend boxoff;
    % key = {'\bfHow to read each panel\rm', ...
    %        'Interaction p: LME, does the arousal-rate slope', ...
    %        sprintf('   differ between %s and %s nights?', lower(gn{1}), lower(gn{2})), ...
    %        'Rate p: LME, arousal-rate slope averaged', ...
    %        sprintf('   over %s and %s nights', lower(gn{1}), lower(gn{2})), ...
    %        '   Model: DV ~ group x arousal rate', ...
    %        '   + (1 | participant), REML, Satterthwaite df', ...
    %        sprintf('q: Benjamini-Hochberg FDR across the %d', n_dv), ...
    %        '   measures, separately for each term', ...
    %        sprintf('\\rho: Spearman, all nights then %s, %s;', lower(gn{1}), lower(gn{2})), ...
    %        '   descriptive only (participants repeat)', ...
    %        'Lines: least-squares fit per group', ...
    %        nights_text};
    % text(0, 0.72, key, 'FontSize', 8.5, 'Units', 'normalized', 'VerticalAlignment', 'top');

    sgtitle(sprintf('%s (n = %d nights, %d participants)', fig_titles{a}, height(P), numel(unique(P.userID))));
    % exportgraphics(gcf, fullfile(out_dir, [file_names{a} tag '.png']), 'Resolution', 300);
end

diary off;