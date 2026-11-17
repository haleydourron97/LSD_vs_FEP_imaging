
%% === Run Cross-Validation Binary PLS Classifier (LSD vs Placebo, Within-Subject) ===
% ===============================================================
% LSD vs Placebo: Binary PLS Classifier with FOLD-WISE Covariate Control
% WITHIN-SUBJECT design — each subject contributes an LSD scan AND a
% Placebo scan. Both scans of a subject are ALWAYS kept in the same fold.
% ===============================================================
%
% Why this differs from the FEP vs HC version:
%   FEP vs HC is BETWEEN-subject (one scan per person), so rows are
%   partitioned directly, stratified by GROUP.
%
%   Here, the LSD and Placebo scans of the same subject share
%   subject-specific FC structure unrelated to the drug contrast. If one
%   scan were in training and the other in test, the model could partly
%   learn "this is subject X" instead of "this is LSD". So the partition is
%   built on unique SUBJECTS (not rows), stratified by SITE (Basel / LPM),
%   and both scans follow their subject. LSD/Placebo balance is automatic
%   because every subject contributes one of each.
%
% Structure matches the FEP vs HC script:
%   - The partition is built ONCE per repeat (rng(r)) and shared by all
%     edge sets, so every model in a repeat uses exactly the same splits.
%   - Fold-level models (BETA, B_cov, x_mean, x_std) are stored on every
%     fold of every repeat and averaged across all K*num_repeats folds at
%     the end to produce ONE averaged model per edge set.
%   - Fold-wise covariate residualization, fold-wise z-scoring and PLS
%     training are identical to the FEP vs HC script.
%
% Configuration
num_repeats = 100;
K = 10;
limit_pls_components = true;
shuffle_labels = false;  % true = within-subject label swap (null model)
verbose_folds = false;   % true = print every fold (326 sets x 100 repeats x 10 folds = 326,000 lines)

% Number of edge/network models
num_edge_sets = length(all_net_edges);
assert(numel(model_labels) == num_edge_sets, ...
    'model_labels (%d) and all_net_edges (%d) differ in length', numel(model_labels), num_edge_sets);
num_subs = length(label_vector); % total rows (subjects x 2 sessions)
assert(num_subs == size(sfc_mat, 3) && num_subs == size(covariates, 1) && num_subs == numel(subject_id_vector), ...
    'label_vector (%d), sfc_mat (%d), covariates (%d) and subject_id_vector (%d) differ', ...
    num_subs, size(sfc_mat, 3), size(covariates, 1), numel(subject_id_vector));

% Unique subjects, fixed order (site-prefixed IDs, e.g. Basel_sub-004,
% LPM_sub-004 — globally unique across sites)
unique_subjects = unique(subject_id_vector, 'stable');
num_unique_subs = length(unique_subjects);
fprintf('%d scans across %d subjects (within-subject): LSD=%d | Placebo=%d\n', ...
    num_subs, num_unique_subs, sum(label_vector == 1), sum(label_vector == 0));

% -------------------------------------------------------
% SITE for each unique subject, from the ID prefix (Basel_ / LPM_)
% -------------------------------------------------------
site_labels_unique = strings(num_unique_subs, 1);
site_labels_unique(startsWith(unique_subjects, "Basel_")) = "Basel";
site_labels_unique(startsWith(unique_subjects, "LPM_"))   = "LPM";
if any(site_labels_unique == "")
    error('Unrecognized site prefix for subject(s): %s', strjoin(unique_subjects(site_labels_unique == ""), ', '));
end
fprintf('Site breakdown: Basel=%d | LPM=%d (total %d)\n', ...
    sum(site_labels_unique == "Basel"), sum(site_labels_unique == "LPM"), num_unique_subs);

%% Initialize performance containers
train_accuracies_all = cell(num_edge_sets, 1);
test_accuracies_all  = cell(num_edge_sets, 1);
per_class_acc_all    = cell(num_edge_sets, 1);

% Scan-level held-out predictions across repeats (both LSD and Placebo rows)
heldout_class_all = cell(num_edge_sets, 1);
heldout_score_all = cell(num_edge_sets, 1);

% Fold-level model containers, to be averaged after all repeats/folds
fold_BETA_all = cell(num_edge_sets, 1);
fold_Bcov_all = cell(num_edge_sets, 1);
fold_mean_all = cell(num_edge_sets, 1);
fold_std_all  = cell(num_edge_sets, 1);

for inet = 1:num_edge_sets
    train_accuracies_all{inet} = [];
    test_accuracies_all{inet}  = [];
    per_class_acc_all{inet}    = [];

    heldout_class_all{inet} = NaN(num_subs, num_repeats);
    heldout_score_all{inet} = NaN(num_subs, num_repeats);

    fold_BETA_all{inet} = [];
    fold_Bcov_all{inet} = [];
    fold_mean_all{inet} = [];
    fold_std_all{inet}  = [];
end

% ===============================================================
% Outer loop over CV repetitions
% ===============================================================
for r = 1:num_repeats
    fprintf('Repetition %d of %d...\n', r, num_repeats);
    rng(r);   % reproducible partitions (and shuffles) per repeat

    % -------------------------------------------------------
    % PAIRED, SITE-STRATIFIED SUBJECT-LEVEL K-FOLD PARTITION
    %
    % Built ONCE per repeat and shared by all edge sets. Partition is on
    % unique SUBJECTS, stratified by SITE; both scans of a subject are then
    % expanded into the same fold below.
    % -------------------------------------------------------
    cv_partition = cvpartition(site_labels_unique, 'KFold', K, 'Stratify', true);

    % Row indices per fold (same for every edge set in this repeat)
    train_idx_all = cell(K, 1);
    test_idx_all  = cell(K, 1);
    for k = 1:K
        train_subj_ids = unique_subjects(training(cv_partition, k));
        test_subj_ids  = unique_subjects(test(cv_partition, k));

        % Sanity check: no subject in both train and test
        overlap = intersect(train_subj_ids, test_subj_ids);
        if ~isempty(overlap)
            error('Subject leakage detected between train/test: %s', strjoin(overlap, ', '));
        end

        % Expand to ROW indices: BOTH sessions per subject
        train_idx_all{k} = ismember(subject_id_vector, train_subj_ids);
        test_idx_all{k}  = ismember(subject_id_vector, test_subj_ids);
    end

    % Loop over each edge set / network model
    for inet = 1:num_edge_sets

        % -----------------------------------------------------------
        % Build RAW feature matrix for this edge set.
        %   rows    = scans (subjects x 2 sessions)
        %   columns = connectivity edges
        % NOT residualized yet — residualization happens inside each fold.
        % -----------------------------------------------------------
        edges_tmp = all_net_edges{inet};
        num_edges = size(edges_tmp, 1);

        X_raw = NaN(num_subs, num_edges);
        for i = 1:num_edges
            X_raw(:, i) = squeeze(sfc_mat(edges_tmp(i, 1), edges_tmp(i, 2), :));
        end

        fold_train_acc     = zeros(K, 1);
        fold_test_acc      = zeros(K, 1);
        per_fold_class_acc = zeros(K, 2); % columns: [Placebo, LSD]

        % Exactly one held-out prediction per scan for this repeat
        heldout_class_this_repeat = NaN(num_subs, 1);
        heldout_score_this_repeat = NaN(num_subs, 1);

        % ===========================================================
        % Inner loop over CV folds
        % ===========================================================
        for k = 1:K

            train_idx = train_idx_all{k};
            test_idx  = test_idx_all{k};

            % -------------------------------------------------------
            % Split RAW connectivity features, labels and covariates
            % (intercept added separately for train and test)
            % -------------------------------------------------------
            X_train_raw = X_raw(train_idx, :);
            X_test_raw  = X_raw(test_idx, :);

            y_train = label_vector(train_idx);
            y_test  = label_vector(test_idx);

            C_train = [ones(sum(train_idx), 1), covariates(train_idx, :)];
            C_test  = [ones(sum(test_idx), 1),  covariates(test_idx, :)];

            % -------------------------------------------------------
            % FOLD-WISE RESIDUALIZATION: NO DATA LEAKAGE
            % Betas fit on TRAINING scans only, applied to train and test.
            % -------------------------------------------------------
            B_cov = C_train \ X_train_raw;
            x_train = X_train_raw - C_train * B_cov;
            x_test  = X_test_raw  - C_test  * B_cov;

            % Optional permutation test (WITHIN-subject):
            % swap LSD/Placebo labels for a random half of the TRAINING
            % subjects, so each subject still has one LSD + one Placebo
            % label (a free randperm across rows would break the pairing).
            if shuffle_labels
                [~, ~, g] = unique(subject_id_vector(train_idx));
                flip = rand(max(g), 1) < 0.5;
                y_train = double(xor(y_train, flip(g)));
            end

            %% ===============================================================
% Average fold-level models across ALL folds and repeats
% (K * num_repeats models per edge set)
%
% Averaging strictly train-only, fold-level models keeps this consistent
% with the FEP/HC pipeline, and means the saved model was never fit on
% 100% of the data.
% ===============================================================
output_dir   = basel_output_dir;          % from the reload block
model_labels = cellstr(model_labels);     % works whether stored as string or cell

avg_BETA_final = cell(num_edge_sets, 1);
avg_B_cov_all  = cell(num_edge_sets, 1);  % [Int, Age, Male, FD] x n_edges
avg_x_mean_all = cell(num_edge_sets, 1);
avg_x_std_all  = cell(num_edge_sets, 1);

for inet = 1:num_edge_sets
    if isempty(fold_BETA_all{inet})
        warning('No fold-level models found for edge set %03d - skipping.', inet);
        continue;
    end

    avg_BETA_final{inet} = mean(fold_BETA_all{inet}, 3);
    avg_B_cov_all{inet}  = mean(fold_Bcov_all{inet}, 3);
    avg_x_mean_all{inet} = mean(fold_mean_all{inet}, 1);
    % Pool SDs via variance rather than averaging SDs directly
    avg_x_std_all{inet}  = sqrt(mean(fold_std_all{inet}.^2, 1));
end
fprintf('Averaged models computed for %d edge sets (from %d fold-level models each)\n', ...
    num_edge_sets, size(fold_BETA_all{1}, 3));

model_path = fullfile(output_dir, 'AveragedPLSModels_LSDvsPlacebo.mat');
save(model_path, 'avg_BETA_final', 'avg_B_cov_all', 'avg_x_mean_all', 'avg_x_std_all', 'model_labels');
fprintf('Saved averaged models to: %s\n', model_path);


%% === Save all subject-level held-out predictions (LSD vs Placebo) ===
% For symptom predictions if we end up doing those.
% (Long-format CSV skipped: ~5M rows; this .mat holds the same information.)
heldout_save_path = fullfile(output_dir, 'BASEL_HeldOutPredictions_AllEdgeSets_LSDvsPlacebo.mat');
save(heldout_save_path, 'heldout_class_all', 'heldout_score_all', ...
    'label_vector', 'subject_id_vector', 'model_labels', '-v7.3');
fprintf('Saved subject-level held-out predictions to: %s\n', heldout_save_path);


%% === Subject-level summary CSV: average across repeats ===
subject_id_cell = cellstr(subject_id_vector);
num_subs        = numel(subject_id_cell);

summary_subjects   = {};
summary_labels     = [];
summary_edge_idx   = [];
summary_edge_names = {};
mean_pred_class    = [];
mean_pred_score    = [];
sd_pred_score      = [];
prop_lsd_class     = [];

for inet = 1:num_edge_sets
    class_mat = heldout_class_all{inet};   % [num_subs x num_repeats]
    score_mat = heldout_score_all{inet};   % [num_subs x num_repeats]

    summary_subjects   = [summary_subjects; subject_id_cell(:)];
    summary_labels     = [summary_labels; label_vector(:)];
    summary_edge_idx   = [summary_edge_idx; repmat(inet, num_subs, 1)];
    summary_edge_names = [summary_edge_names; repmat(model_labels(inet), num_subs, 1)];

    mean_pred_class = [mean_pred_class; mean(class_mat, 2, 'omitnan')];
    mean_pred_score = [mean_pred_score; mean(score_mat, 2, 'omitnan')];
    sd_pred_score   = [sd_pred_score; std(score_mat, 0, 2, 'omitnan')];
    prop_lsd_class  = [prop_lsd_class; mean(class_mat == 1, 2, 'omitnan')];   % 1 = LSD
end

summary_groups = repmat({''}, numel(summary_labels), 1);
summary_groups(summary_labels == 0) = {'Placebo'};
summary_groups(summary_labels == 1) = {'LSD'};

T_summary = table(summary_subjects, summary_labels, summary_groups, summary_edge_idx, ...
    summary_edge_names, mean_pred_class, prop_lsd_class, mean_pred_score, sd_pred_score, ...
    'VariableNames', {'SubjectID', 'TrueLabel', 'TrueGroup', 'EdgeSetIndex', 'EdgeSetLabel', ...
                      'MeanHeldOutPredClass', 'PropRepeatsClassifiedLSD', ...
                      'MeanHeldOutPredScore', 'SDHeldOutPredScore'});

summary_csv_path = fullfile(output_dir, 'BASEL_HeldOutPredictions_SubjectSummary_LSDvsPlacebo.csv');
writetable(T_summary, summary_csv_path);
fprintf('Saved subject-summary CSV: %s\n', summary_csv_path);


%% =========================================================================
% PER-ITERATION METRICS (one row per model x repeat) — LSD vs Placebo
% -------------------------------------------------------------------------
% Collapses each repeat's K folds into one repeat-level value; d' per repeat
% from the scan-level held-out predictions.
% per_class_acc_all columns are [Placebo, LSD] (class 0, class 1).
% Hit = LSD called LSD; False alarm = Placebo called LSD.
% =========================================================================
num_models = numel(model_labels);
n_repeats_observed = numel(test_accuracies_all{1}) / K;
n_iter_rows = num_models * n_repeats_observed;

Model       = cell(n_iter_rows, 1);
EdgeSet     = NaN(n_iter_rows, 1);
Repeat      = NaN(n_iter_rows, 1);
TestAcc     = NaN(n_iter_rows, 1);
Placebo_Acc = NaN(n_iter_rows, 1);
LSD_Acc     = NaN(n_iter_rows, 1);
Dprime      = NaN(n_iter_rows, 1);

row_idx = 1;
for inet = 1:num_models
    % --- collapse K-fold accuracy down to one value per repeat ---
    test_vals = test_accuracies_all{inet};             % [K*num_repeats x 1]
    n_rep     = numel(test_vals) / K;
    test_repeat_means = mean(reshape(test_vals, K, n_rep), 1)';

    class_vals = per_class_acc_all{inet};              % [K*num_repeats x 2]
    placebo_repeat_means = mean(reshape(class_vals(:,1), K, n_rep), 1)';
    lsd_repeat_means     = mean(reshape(class_vals(:,2), K, n_rep), 1)';

    % --- d' per repeat, from scan-level held-out predictions ---
    pred_mat = heldout_class_all{inet};                % [num_subs x num_repeats]
    dprime_vals = NaN(n_rep, 1);
    for r = 1:n_rep
        pred = pred_mat(:, r);
        hits   = sum(pred == 1 & label_vector == 1);   % LSD called LSD
        misses = sum(pred == 0 & label_vector == 1);   % LSD missed
        fas    = sum(pred == 1 & label_vector == 0);   % Placebo called LSD
        crs    = sum(pred == 0 & label_vector == 0);   % Placebo called Placebo

        % log-linear correction avoids +/-Inf when hit/FA rate = 0 or 1
        hit_rate = (hits + 0.5) / (hits + misses + 1);
        fa_rate  = (fas  + 0.5) / (fas + crs + 1);
        dprime_vals(r) = norminv(hit_rate) - norminv(fa_rate);
    end

    rows = row_idx:(row_idx + n_rep - 1);
    Model(rows)       = repmat(model_labels(inet), n_rep, 1);
    EdgeSet(rows)     = inet;
    Repeat(rows)      = (1:n_rep)';
    TestAcc(rows)     = test_repeat_means;
    Placebo_Acc(rows) = placebo_repeat_means;
    LSD_Acc(rows)     = lsd_repeat_means;
    Dprime(rows)      = dprime_vals;
    row_idx = row_idx + n_rep;
end

per_iteration_summary = table(EdgeSet, Model, Repeat, TestAcc, Placebo_Acc, LSD_Acc, Dprime);
per_iteration_summary{:, {'TestAcc','Placebo_Acc','LSD_Acc'}} = ...
    per_iteration_summary{:, {'TestAcc','Placebo_Acc','LSD_Acc'}} * 100;  % accuracy as %, d' unchanged

per_iteration_path = fullfile(output_dir, 'BASEL_FinalSummary_LSDvsPlacebo_PerIteration_WithDprime.csv');
writetable(per_iteration_summary, per_iteration_path);
fprintf('Per-iteration summary (%d rows) saved to: %s\n', n_iter_rows, per_iteration_path);


%% =========================================================================
% PER-MODEL SUMMARY (mean / SD / 95% CI across repeats) — LSD vs Placebo
% -------------------------------------------------------------------------
% CI is descriptive only: it reflects variability across CV partitions,
% not independent-sample uncertainty (repeats are not independent).
% =========================================================================
metrics = {'TestAcc', 'Placebo_Acc', 'LSD_Acc', 'Dprime'};

[~, first_idx] = unique(per_iteration_summary.EdgeSet, 'stable');
EdgeSet_s = per_iteration_summary.EdgeSet(first_idx);
Model_s   = per_iteration_summary.Model(first_idx);
num_models_summary = numel(EdgeSet_s);

for m = 1:numel(metrics)
    mean_col.(metrics{m}) = NaN(num_models_summary, 1);
    sd_col.(metrics{m})   = NaN(num_models_summary, 1);
    ci_lo.(metrics{m})    = NaN(num_models_summary, 1);
    ci_hi.(metrics{m})    = NaN(num_models_summary, 1);
end

for i = 1:num_models_summary
    rows  = per_iteration_summary.EdgeSet == EdgeSet_s(i);
    n_rep = sum(rows);

    for m = 1:numel(metrics)
        vals = per_iteration_summary.(metrics{m})(rows);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 0, 'omitnan');
        if n_rep > 1
            half_width = tinv(0.975, n_rep - 1) * sd / sqrt(n_rep);
        else
            half_width = NaN;
        end
        mean_col.(metrics{m})(i) = mu;
        sd_col.(metrics{m})(i)   = sd;
        ci_lo.(metrics{m})(i)    = mu - half_width;
        ci_hi.(metrics{m})(i)    = mu + half_width;
    end
end

per_model_summary = table(EdgeSet_s, Model_s, ...
    mean_col.TestAcc,     sd_col.TestAcc,     ci_lo.TestAcc,     ci_hi.TestAcc, ...
    mean_col.Placebo_Acc, sd_col.Placebo_Acc, ci_lo.Placebo_Acc, ci_hi.Placebo_Acc, ...
    mean_col.LSD_Acc,     sd_col.LSD_Acc,     ci_lo.LSD_Acc,     ci_hi.LSD_Acc, ...
    mean_col.Dprime,      sd_col.Dprime,      ci_lo.Dprime,      ci_hi.Dprime, ...
    'VariableNames', { ...
        'EdgeSet', 'Model', ...
        'TestAcc_Mean',     'TestAcc_SD',     'TestAcc_CI_Lo',     'TestAcc_CI_Hi', ...
        'Placebo_Acc_Mean', 'Placebo_Acc_SD', 'Placebo_Acc_CI_Lo', 'Placebo_Acc_CI_Hi', ...
        'LSD_Acc_Mean',     'LSD_Acc_SD',     'LSD_Acc_CI_Lo',     'LSD_Acc_CI_Hi', ...
        'Dprime_Mean',      'Dprime_SD',      'Dprime_CI_Lo',      'Dprime_CI_Hi'});

per_model_path = fullfile(output_dir, 'BASEL_FinalSummary_LSDvsPlacebo_PerModel_MeanAcrossRepeats.csv');
writetable(per_model_summary, per_model_path);
fprintf('Per-model summary (mean/SD/CI across %d repeats) saved to: %s\n', ...
    n_repeats_observed, per_model_path);

            % -------------------------------------------------------
            % Z-SCORING: ALSO FIT ON TRAINING DATA ONLY
            % -------------------------------------------------------
            x_mean = mean(x_train, 1);
            x_std  = std(x_train, 0, 1);
            x_std(x_std == 0) = 1; % avoid division by zero for constant edges

            x_train = (x_train - x_mean) ./ x_std;
            x_test  = (x_test  - x_mean) ./ x_std;

            % -------------------------------------------------------
            % Train PLS classifier on TRAINING subjects only
            % -------------------------------------------------------
            if limit_pls_components
                ncomp = min([size(x_train, 1)-1, size(x_train, 2), 3]);
                [~, ~, ~, ~, BETA] = plsregress(x_train, y_train, ncomp);
            else
                [~, ~, ~, ~, BETA] = plsregress(x_train, y_train);
            end

            % -------------------------------------------------------
            % Predict held-out TEST subjects (both of their scans)
            % -------------------------------------------------------
            y_pred_test = [ones(size(x_test, 1), 1), x_test] * BETA;
            pred_class_test = double(y_pred_test >= 0.5);

            heldout_score_this_repeat(test_idx) = y_pred_test;
            heldout_class_this_repeat(test_idx) = pred_class_test;

            % Test / training accuracy
            fold_test_acc(k) = mean(pred_class_test == y_test);
            y_pred_train = [ones(size(x_train, 1), 1), x_train] * BETA;
            fold_train_acc(k) = mean(double(y_pred_train >= 0.5) == y_train);

            % Per-class test accuracy
            for cls = 0:1
                total = sum(y_test == cls);
                correct = sum((pred_class_test == cls) & (y_test == cls));
                per_fold_class_acc(k, cls+1) = correct / max(total, 1);
            end

            if verbose_folds, fprintf(' Edge %d (%s) | Fold %d — TrainAcc: %.3f | TestAcc: %.3f | Placebo: %.3f, LSD: %.3f\n', ...
                inet, model_labels{inet}, k, fold_train_acc(k), fold_test_acc(k), ...
                per_fold_class_acc(k,1), per_fold_class_acc(k,2)); end

            % -------------------------------------------------------
            % Store this fold's model for later averaging
            % -------------------------------------------------------
            fold_BETA_all{inet} = cat(3, fold_BETA_all{inet}, BETA);   % [p+1 x 1 x n_folds_total]
            fold_Bcov_all{inet} = cat(3, fold_Bcov_all{inet}, B_cov);  % [n_cov+1 x p x n_folds_total]
            fold_mean_all{inet} = [fold_mean_all{inet}; x_mean];       % [n_folds_total x p]
            fold_std_all{inet}  = [fold_std_all{inet};  x_std];        % [n_folds_total x p]
        end

        % Sanity check: every scan should have one held-out prediction
        if any(isnan(heldout_class_this_repeat)) || any(isnan(heldout_score_this_repeat))
            warning('Missing held-out predictions in Repeat %d, Edge Set %d', r, inet);
        end

        % Store this repeat's scan-level predictions
        heldout_class_all{inet}(:, r) = heldout_class_this_repeat;
        heldout_score_all{inet}(:, r) = heldout_score_this_repeat;

        % Append fold-wise performance
        train_accuracies_all{inet} = [train_accuracies_all{inet}; fold_train_acc];
        test_accuracies_all{inet}  = [test_accuracies_all{inet}; fold_test_acc];
        per_class_acc_all{inet}    = [per_class_acc_all{inet}; per_fold_class_acc];
    end
end

%% =========================================================================
% PERMUTATION TEST — LSD vs Placebo (Basel), within-subject label swaps
% -------------------------------------------------------------------------
% Run AFTER the reload block (needs sfc_mat, covariates, label_vector,
% subject_id_vector, all_net_edges, model_labels, basel_output_dir).
%
% Design:
%   - ONE fixed 10-fold partition (subject-level, site-stratified, seeded).
%     Fold-wise residualization + z-scoring do not use labels, so they are
%     computed ONCE per edge set and reused for every permutation; only the
%     PLS fit is redone. This makes 5,000 permutations feasible.
%   - Permutation 0 = TRUE labels on the same partition (the observed value
%     the null is compared against).
%   - Permutations 1..n_perm: for each subject, LSD/Placebo labels are
%     swapped with probability 0.5 (keeps the paired design intact).
%     The SAME label permutations are used for every edge set.
%   - p = (#null >= observed + 1) / (n_perm + 1), one-sided; BH-FDR across
%     all edge sets.
%   - Checkpoints after every edge set, so a crash / forced restart resumes
%     where it stopped.
% =========================================================================
n_perm       = 5000;
K            = 10;
cv_seed      = 2026;     % fixed partition
perm_seed    = 42;       % fixed label permutations
max_ncomp    = 3;
output_dir   = basel_output_dir;
model_labels = cellstr(model_labels);
checkpoint_path = fullfile(output_dir, 'BASEL_PermutationTest_LSDvsPlacebo_checkpoint.mat');

num_edge_sets = numel(all_net_edges);
num_subs      = numel(label_vector);

%% === Fixed subject-level, site-stratified partition ===
unique_subjects = unique(subject_id_vector, 'stable');
num_unique_subs = numel(unique_subjects);
site_labels_unique = strings(num_unique_subs, 1);
site_labels_unique(startsWith(unique_subjects, "Basel_")) = "Basel";
site_labels_unique(startsWith(unique_subjects, "LPM_"))   = "LPM";

rng(cv_seed);
cv_partition = cvpartition(site_labels_unique, 'KFold', K, 'Stratify', true);
train_idx_all = cell(K, 1);
test_idx_all  = cell(K, 1);
for k = 1:K
    train_idx_all{k} = ismember(subject_id_vector, unique_subjects(training(cv_partition, k)));
    test_idx_all{k}  = ismember(subject_id_vector, unique_subjects(test(cv_partition, k)));
end

%% === Fixed within-subject label permutations ===
[~, subj_row] = ismember(subject_id_vector, unique_subjects);   % subject index per scan
rng(perm_seed);
flips = rand(num_unique_subs, n_perm) < 0.5;                     % true = swap this subject
Y_all = [label_vector(:), double(xor(label_vector(:), flips(subj_row, :)))];  % [scans x (1+n_perm)], col 1 = true labels

%% === Output containers (or resume from checkpoint) ===
if exist(checkpoint_path, 'file')
    load(checkpoint_path, 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', 'done');
    fprintf('Resuming from checkpoint: %d of %d edge sets already done.\n', sum(done), num_edge_sets);
else
    obs_acc     = nan(num_edge_sets, 1);
    obs_dprime  = nan(num_edge_sets, 1);
    null_acc    = nan(num_edge_sets, n_perm, 'single');
    null_dprime = nan(num_edge_sets, n_perm, 'single');
    done        = false(num_edge_sets, 1);
end

z_rate = @(k, n) norminv((k + 0.5) ./ (n + 1));   % log-linear correction
n1 = sum(label_vector == 1);   % LSD scans (unchanged by within-subject swaps)
n0 = sum(label_vector == 0);   % Placebo scans

%% === Main loop: edge sets ===
t_start = tic;
for inet = 1:num_edge_sets
    if done(inet), continue; end
    t_set = tic;

    % --- Raw features for this edge set ---
    edges_tmp = all_net_edges{inet};
    X_raw = NaN(num_subs, size(edges_tmp, 1));
    for i = 1:size(edges_tmp, 1)
        X_raw(:, i) = squeeze(sfc_mat(edges_tmp(i,1), edges_tmp(i,2), :));
    end

    % --- Fold-wise residualization + z-scoring (label-free, done ONCE) ---
    x_train_all = cell(K, 1);
    x_test_all  = cell(K, 1);
    for k = 1:K
        tr = train_idx_all{k};  te = test_idx_all{k};
        C_train = [ones(sum(tr), 1), covariates(tr, :)];
        C_test  = [ones(sum(te), 1), covariates(te, :)];
        B_cov   = C_train \ X_raw(tr, :);
        x_tr    = X_raw(tr, :) - C_train * B_cov;
        x_te    = X_raw(te, :) - C_test  * B_cov;
        mu = mean(x_tr, 1);  sd = std(x_tr, 0, 1);  sd(sd == 0) = 1;
        x_train_all{k} = (x_tr - mu) ./ sd;
        x_test_all{k}  = (x_te - mu) ./ sd;
    end

    % --- Permutation 0 (true labels) + n_perm null permutations ---
    acc_tmp = nan(1, n_perm + 1);
    dp_tmp  = nan(1, n_perm + 1);
    parfor p = 1:(n_perm + 1)          % runs as a normal for-loop without the Parallel Toolbox
        y = Y_all(:, p);
        pred = nan(num_subs, 1);
        for k = 1:K
            tr = train_idx_all{k};  te = test_idx_all{k};
            ncomp = min([sum(tr)-1, size(x_train_all{k}, 2), max_ncomp]);
            [~, ~, ~, ~, BETA] = plsregress(x_train_all{k}, y(tr), ncomp);
            pred(te) = double([ones(sum(te), 1), x_test_all{k}] * BETA >= 0.5);
        end
        acc_tmp(p) = mean(pred == y);
        dp_tmp(p)  = z_rate(sum(pred == 1 & y == 1), n1) - z_rate(sum(pred == 1 & y == 0), n0);
    end

    obs_acc(inet)        = acc_tmp(1);
    obs_dprime(inet)     = dp_tmp(1);
    null_acc(inet, :)    = acc_tmp(2:end);
    null_dprime(inet, :) = dp_tmp(2:end);
    done(inet) = true;

    % --- Checkpoint + progress ---
    save(checkpoint_path, 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', 'done', '-v7.3');
    elapsed = toc(t_start);
    n_done_now = sum(done);
    fprintf('[%03d/%03d] %-45s d''=%5.2f acc=%5.1f%% | %.0fs this set | ~%.1f h left\n', ...
        inet, num_edge_sets, model_labels{inet}, obs_dprime(inet), 100*obs_acc(inet), ...
        toc(t_set), elapsed / max(n_done_now, 1) * (num_edge_sets - n_done_now) / 3600);
end

%% === p-values + BH-FDR across edge sets ===
p_acc    = (sum(null_acc    >= obs_acc,    2) + 1) / (n_perm + 1);
p_dprime = (sum(null_dprime >= obs_dprime, 2) + 1) / (n_perm + 1);

q_acc = nan(num_edge_sets, 1);  q_dprime = nan(num_edge_sets, 1);
[ps, idx] = sort(p_acc);     q_acc(idx)    = min(1, cummin(ps .* num_edge_sets ./ (1:num_edge_sets)', 'reverse'));
[ps, idx] = sort(p_dprime);  q_dprime(idx) = min(1, cummin(ps .* num_edge_sets ./ (1:num_edge_sets)', 'reverse'));

%% === Save ===
perm_results = table((1:num_edge_sets)', model_labels(:), ...
    100*obs_acc, p_acc, q_acc, obs_dprime, p_dprime, q_dprime, ...
    100*mean(null_acc, 2), mean(null_dprime, 2), ...
    'VariableNames', {'EdgeSet', 'Model', ...
    'Acc_FixedPartition', 'Acc_perm_p', 'Acc_perm_q', ...
    'Dprime_FixedPartition', 'Dprime_perm_p', 'Dprime_perm_q', ...
    'NullAcc_Mean', 'NullDprime_Mean'});

writetable(perm_results, fullfile(output_dir, 'BASEL_PermutationTest_LSDvsPlacebo_Summary.csv'));
save(fullfile(output_dir, 'BASEL_PermutationTest_LSDvsPlacebo.mat'), ...
    'perm_results', 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', ...
    'p_acc', 'p_dprime', 'q_acc', 'q_dprime', 'flips', 'n_perm', 'K', ...
    'cv_seed', 'perm_seed', 'train_idx_all', 'test_idx_all', 'model_labels', '-v7.3');

fprintf('Done. %d of %d edge sets significant (d'', FDR q < .05)\n', sum(q_dprime < .05), num_edge_sets);