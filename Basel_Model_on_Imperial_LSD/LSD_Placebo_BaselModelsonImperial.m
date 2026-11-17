% Applying Basel Model (within-subject) to Imperial data (within-subject)

%% Making matrices for control groups 
% Define directories for LSD datasets
Imperial_LSD_dir = '/Volumes/toshiba/LSD_November11/Imperial_Only/Imperial/NILEARN_Processed_Pearson_Z_transformed_w_Motion';
covariates_dir = '/Volumes/toshiba/LSD_November11/Imperial_Only/Imperial/Covariates';
output_dir_Imperial = '/Volumes/toshiba/LSD_November11/BaselOnImperial';

%% Setup Paths
excluded_ids_Imperial = {'018', '019', '007', '014'}; % Subjects to excluded 
Imperial_LSD_cov_file = fullfile(covariates_dir, 'Imperial_LSD_Covariates_RUN1w_AGE_GENDER.csv'); % LSD covariate files
Imperial_Placebo_cov_file = fullfile(covariates_dir, 'Imperial_Placebo_Covariates_RUN1w_AGE_GENDER.csv'); % Placebo covariate files


%% =========================== LSD LOADING =========================== %%
Imperial_LSD_cov_table = readtable(Imperial_LSD_cov_file);
imperial_lsd_files = dir(fullfile(Imperial_LSD_dir, 'sub-*_ses-LSD_run-1_combinedsubcorticalOnly_TS_denoised_NO_GSR__MOTION_NILEARN_z_transformed_matrix.csv'));
all_ids_imperial_lsd = regexp({imperial_lsd_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_imperial_lsd = unique(cellfun(@(x) x{1}, all_ids_imperial_lsd, 'UniformOutput', false));

LSD_UsedImperial = sort(setdiff(all_ids_imperial_lsd, excluded_ids_Imperial));
fprintf('Imperial LSD: Excluding subjects %s\n', strjoin(excluded_ids_Imperial, ', '));
num_Imperial_LSD = length(LSD_UsedImperial);

fc_matrices_Imperial_LSD = cell(num_Imperial_LSD, 1);
subject_ids_Imperial_LSD = strings(num_Imperial_LSD, 1);
mean_fd_Imperial_LSD = nan(num_Imperial_LSD, 1);
age_Imperial_LSD = nan(num_Imperial_LSD, 1);
male_Imperial_LSD = nan(num_Imperial_LSD, 1);

for i = 1:num_Imperial_LSD
    sub_id = LSD_UsedImperial{i};
    pattern = sprintf('sub-%s_ses-LSD_run-1_combinedsubcorticalOnly_TS_denoised_NO_GSR__MOTION_NILEARN_z_transformed_matrix.csv', sub_id);
    file = dir(fullfile(Imperial_LSD_dir, pattern));
    if isempty(file)
        warning('File not found for Imperial subject %s', sub_id);
        continue;
    end
    file_path = fullfile(Imperial_LSD_dir, file(1).name);
    fc_matrices_Imperial_LSD{i} = readmatrix(file_path);
    subject_ids_Imperial_LSD(i) = string(['Imperial_sub-' sub_id]);

    cov_idx = find(strcmp(Imperial_LSD_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    mean_fd_Imperial_LSD(i) = Imperial_LSD_cov_table.Mean_framewise_displacement(cov_idx);
    age_Imperial_LSD(i) = Imperial_LSD_cov_table.Age(cov_idx);
    male_Imperial_LSD(i) = Imperial_LSD_cov_table.Male(cov_idx);
end

save(fullfile(output_dir_Imperial, 'Imperial_LSD_Selections.mat'), ...
    'LSD_UsedImperial', 'fc_matrices_Imperial_LSD', 'subject_ids_Imperial_LSD', ...
    'mean_fd_Imperial_LSD', 'age_Imperial_LSD', 'male_Imperial_LSD');

fprintf('Loaded %d Imperial LSD matrices and matched covariates.\n', num_Imperial_LSD);


%% =========================== PLACEBO LOADING =========================== %%
Imperial_Placebo_cov_table = readtable(Imperial_Placebo_cov_file);

imperial_placebo_files = dir(fullfile(Imperial_LSD_dir, 'sub-*_ses-plcb_run-1_combinedsubcorticalOnly_TS_denoised_NO_GSR__MOTION_NILEARN_z_transformed_matrix.csv'));
all_ids_imperial_plcb = regexp({imperial_placebo_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_imperial_plcb = unique(cellfun(@(x) x{1}, all_ids_imperial_plcb, 'UniformOutput', false));

Placebo_UsedImperial = sort(setdiff(all_ids_imperial_plcb, excluded_ids_Imperial));
%fprintf('Imperial Placebo: Excluding subjects %s\n', strjoin(excluded_ids_Imperial, ', '));
num_Imperial_Placebo = length(Placebo_UsedImperial);

fc_matrices_Imperial_Placebo = cell(num_Imperial_Placebo, 1);
subject_ids_Imperial_Placebo = strings(num_Imperial_Placebo, 1);
mean_fd_Imperial_Placebo = nan(num_Imperial_Placebo, 1);
age_Imperial_Placebo = nan(num_Imperial_Placebo, 1);
male_Imperial_Placebo = nan(num_Imperial_Placebo, 1);

for i = 1:num_Imperial_Placebo
    sub_id = Placebo_UsedImperial{i};
    pattern = sprintf('sub-%s_ses-plcb_run-1_combinedsubcorticalOnly_TS_denoised_NO_GSR__MOTION_NILEARN_z_transformed_matrix.csv', sub_id);
    file = dir(fullfile(Imperial_LSD_dir, pattern));
    if isempty(file)
        warning('File not found for Imperial subject %s', sub_id);
        continue;
    end
    file_path = fullfile(Imperial_LSD_dir, file(1).name);
    fc_matrices_Imperial_Placebo{i} = readmatrix(file_path);
    subject_ids_Imperial_Placebo(i) = string(['Imperial_sub-' sub_id]);

    cov_idx = find(strcmp(Imperial_Placebo_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    mean_fd_Imperial_Placebo(i) = Imperial_Placebo_cov_table.Mean_framewise_displacement(cov_idx);
    age_Imperial_Placebo(i) = Imperial_Placebo_cov_table.Age(cov_idx);
    male_Imperial_Placebo(i) = Imperial_Placebo_cov_table.Male(cov_idx);
end

save(fullfile(output_dir_Imperial, 'Imperial_Placebo_Selections.mat'), ...
    'Placebo_UsedImperial', 'fc_matrices_Imperial_Placebo', 'subject_ids_Imperial_Placebo', ...
    'mean_fd_Imperial_Placebo', 'age_Imperial_Placebo', 'male_Imperial_Placebo');

fprintf('Loaded %d Imperial placebo matrices and matched covariates.\n', num_Imperial_Placebo);


%% === Final Output for LSD vs Placebo  ===
% Combine FC matrices
fc_matrices_all = [fc_matrices_Imperial_LSD; ...
                   fc_matrices_Imperial_Placebo];
% Combine subject IDs
subject_ids = [subject_ids_Imperial_LSD; ...
               subject_ids_Imperial_Placebo];
% Combine ALL covariates: Age, Male, Mean FD (same order)
age_all = [age_Imperial_LSD; ...
           age_Imperial_Placebo];
male_all = [male_Imperial_LSD; ...
            male_Imperial_Placebo];
mean_fd_all = [mean_fd_Imperial_LSD; ...
               mean_fd_Imperial_Placebo];
% Subject ID vector for within-subject pairing in CV folds
subject_id_vector = subject_ids;

%% === Define Binary Class Labels for LSD and Placebo ===
num_LSD     = num_Imperial_LSD;
num_Placebo = num_Imperial_Placebo;

% Total subject count
num_subs = num_LSD + num_Placebo;

% label count must match loaded FC matrices
if num_subs ~= length(fc_matrices_all)
    error('Mismatch in subject count: expected %d, found %d matrices', num_subs, length(fc_matrices_all));
end

% Create binary class label vector: 1 = LSD, 0 = Placebo 
label_vector = [ ...
    ones(num_LSD, 1); ...     % LSD = 1
    zeros(num_Placebo, 1) ... % Placebo = 0
];

fprintf('Created binary label vector [%d subjects] (LSD=1, Placebo=0)\n', num_subs);
fprintf('   LSD: %d | Placebo: %d\n', num_LSD, num_Placebo);

%% === Additional sanity checks and final save ===
if num_subs ~= length(subject_ids) || num_subs ~= length(age_all) ...
        || num_subs ~= length(male_all) ...
        || num_subs ~= length(mean_fd_all) || num_subs ~= length(label_vector)
    error('Mismatch in number of subjects across matrices, covariates, labels, or IDs');
end

% Within-subject pairing: row i of the LSD block and row i of the Placebo
% block must be the same subject
assert(isequal(LSD_UsedImperial(:), Placebo_UsedImperial(:)), 'Imperial LSD and Placebo subject lists differ');

% NaN covariates would make residuals NaN
assert(~any(isnan([age_all; male_all; mean_fd_all])), 'NaNs present in Imperial covariates');

% Build final FC array: ROI x ROI x subjects
sfc_mat = cat(3, fc_matrices_all{:});

% Build covariates matrix: columns = [Age, Male, Mean_FD]
covariates = [age_all, male_all, mean_fd_all];

% Save outputs
save(fullfile(output_dir_Imperial, 'Imperial_sfc_mat_LSD_vs_Placebo.mat'), 'sfc_mat');
save(fullfile(output_dir_Imperial, 'Imperial_covariates_LSD_vs_Placebo.mat'), 'covariates');
save(fullfile(output_dir_Imperial, 'Imperial_label_vector_LSD_vs_Placebo.mat'), 'label_vector');
save(fullfile(output_dir_Imperial, 'Imperial_subject_id_vector_LSD_vs_Placebo.mat'), 'subject_id_vector');
writematrix(subject_ids, fullfile(output_dir_Imperial, 'Imperial_subject_ids_LSD_vs_Placebo.txt'));

fprintf('Saved sFC matrix, full covariates (Age/Male/FD), labels, and subject IDs for Imperial LSD vs Placebo classification.\n');


%% === Define Paths: Edge Sets + Model Labels (Schaefer100 + Tian38) ===
edge_dir  = '/Volumes/toshiba/LSD_November11/Within_Between_allEdge/Schaefer100_Tian38';
edge_file = fullfile(edge_dir, 'edge_sets_and_labels.mat');

%% BASEL AVERAGED MODEL
% Folder with AveragedPLSModels_LSDvsPlacebo_1comp/_2comp/_3comp.mat
% (used by the 1-2-3 component application script)
model_folder = '/Volumes/toshiba/LSD_November11/LSD_BASELS/OutputBaselOnly/schaeferSubcortical_Nov17';
% sfc_mat (Imperial, ROI x ROI x scans) was built above
num_rois = size(sfc_mat, 1);

%% === Load Edge Sets + Model Labels (built and saved TOGETHER) ===
if ~exist(edge_file, 'file'), error('Missing file: %s', edge_file); end
E = load(edge_file, 'edge_sets', 'model_labels', 'model_type', ...
    'model_i', 'model_j', 'group_labels', 'roi_group');

all_net_edges = E.edge_sets;          % 326 x 1: full, 25 within, 300 between
model_labels  = E.model_labels;       % 326 x 1, same order as all_net_edges
model_type    = E.model_type;         % 'full' / 'within' / 'between'
group_labels  = E.group_labels;       % the 25 group names
num_nets      = numel(group_labels);

%% === Sanity Checks ===
% 1. Atlas matches the Imperial FC matrices (combinedsubcorticalOnly = 138 ROIs)
assert(num_rois == numel(E.roi_group), ...
    'Imperial sfc_mat is %dx%d but edge sets are for %d ROIs (Schaefer100+Tian38 expected)', ...
    num_rois, num_rois, numel(E.roi_group));

% 2. One label per edge set, and no edge index out of range
assert(numel(all_net_edges) == numel(model_labels), 'Mismatch between model labels and edge sets!');
assert(max(cellfun(@(e) max(e(:)), all_net_edges)) <= num_rois, 'Edge index exceeds number of ROIs');

% 3. Every set contains only edges from the groups its label names
for m = 2:numel(all_net_edges)
    e    = all_net_edges{m};
    pair = sort([E.roi_group(e(:,1)), E.roi_group(e(:,2))], 2);
    assert(all(pair(:,1) == E.model_i(m) & pair(:,2) == E.model_j(m)), ...
        'Edge set %d (%s) contains edges from the wrong groups', m, model_labels{m});
end

fprintf('Loaded %d edge sets (%d full, %d within, %d between) for %d ROIs\n', ...
    numel(all_net_edges), sum(strcmp(model_type,'full')), sum(strcmp(model_type,'within')), ...
    sum(strcmp(model_type,'between')), num_rois);
fprintf('%d model labels match %d edge sets (verified edge-by-edge).\n', ...
    numel(model_labels), numel(all_net_edges));

%% =========================================================================
% APPLY BASEL LSD vs PLACEBO CLASSIFIERS (1, 2 AND 3 PLS COMPONENTS) TO Imperial
% -------------------------------------------------------------------------
% Run AFTER the Imperial loading block (needs sfc_mat, covariates, label_vector,
% num_subs, num_LSD, num_Placebo, LSD_UsedImperial, Placebo_UsedImperial,
% all_net_edges, model_labels, output_dir_Imperial).
%
% For each component count in ncomp_list:
%   1. Load AveragedPLSModels_LSDvsPlacebo_<n>comp.mat (from the Basel CV run)
%   2. Residualize Imperial (resid_mode, same as before)
%   3. Apply the model with Basel z-scoring AND within-Imperial z-scoring
%   4. Permutation test (same 5,000 within-subject swaps for every ncomp,
%      every edge set and both z-scoring modes)
% Outputs are tagged _1comp / _2comp / _3comp; a comparison table across
% component counts is written at the end (one per z-scoring mode).
% =========================================================================

%% === Imperial variable names ===
% The loading block may still use the older RCH names; map them here so
% everything below (and every output file) says Imperial.
if ~exist('output_dir_Imperial', 'var'),  output_dir_Imperial  = output_dir_RCH;  end
if ~exist('LSD_UsedImperial', 'var'),     LSD_UsedImperial     = LSD_UsedRCH;     end
if ~exist('Placebo_UsedImperial', 'var'), Placebo_UsedImperial = Placebo_UsedRCH; end

%% === Settings ===
% Folder holding AveragedPLSModels_LSDvsPlacebo_1comp.mat / _2comp / _3comp
% (the Basel CV output folder, basel_output_dir in the Basel script)
model_folder = '/Volumes/toshiba/LSD_November11/LSD_BASELS/OutputBaselOnly/schaeferSubcortical_Nov17';
ncomp_list   = [1 2 3];
resid_mode   = 'basel';        % 'basel' (primary) | 'hybrid' | 'local' — see notes below
zscore_modes = {'basel_zscore', 'within_zscore'};
n_perm       = 5000;

model_labels  = cellstr(model_labels);
num_edge_sets = numel(all_net_edges);
assert(numel(model_labels) == num_edge_sets, 'model_labels and all_net_edges differ in length');

% Check all model files exist before starting
for c = 1:numel(ncomp_list)
    f = fullfile(model_folder, sprintf('AveragedPLSModels_LSDvsPlacebo_%dcomp.mat', ncomp_list(c)));
    if ~exist(f, 'file'), error('Missing Basel model file: %s', f); end
end

%% === PERMUTATION SETUP (within-subject label swaps) — built ONCE ===
% The Basel model is fixed, so only the Imperial labels are permuted.
% Residualization and z-scoring do not use labels, so they are not redone.
% Each permutation swaps LSD/Placebo for each subject with probability 0.5,
% keeping the paired design intact. The SAME permutations are used for
% every component count, every edge set AND both z-scoring modes.

% Pairing check: row i of the LSD block and row i of the Placebo block must be
% the same subject (both lists are sorted, so this holds if the sets match)
assert(isequal(LSD_UsedImperial(:), Placebo_UsedImperial(:)), ...
    'LSD and Placebo subject lists differ — within-subject permutation needs complete pairs');

n1       = num_LSD;                         % LSD scans
n0       = num_Placebo;                     % Placebo scans
n_pairs  = n1;
pair_idx = [(1:n_pairs)'; (1:n_pairs)'];    % subject index for every row

rng(42);                                    % reproducible permutations
flips = rand(n_pairs, n_perm) < 0.5;        % true = swap this subject's labels
perm_labels = double(xor(label_vector, flips(pair_idx, :)));   % [N_scans x n_perm]

% d' with log-linear correction (avoids Inf when a rate is 0% or 100%)
z_rate = @(k, n) norminv((k + 0.5) ./ (n + 1));
p_fun  = @(null, obs) (sum(null >= obs) + 1) / (n_perm + 1);   % one-sided

% Design matrix: one row per Imperial scan, columns = [Intercept, Age, Male, FD]
% (the Imperial scan's own covariate values)
cov_X = [ones(num_subs, 1), covariates];

res = struct();   % res.comp1.basel_zscore, res.comp2.within_zscore, ...

%% =========================================================================
% OUTER LOOP: NUMBER OF PLS COMPONENTS
% =========================================================================
for icn = 1:numel(ncomp_list)
    ncomp     = ncomp_list(icn);
    ncomp_tag = sprintf('_%dcomp', ncomp);
    comp_key  = sprintf('comp%d', ncomp);
    fprintf('\n#################### %d PLS COMPONENT(S) ####################\n', ncomp);

    %% === Load Basel Averaged Classifier for Transfer ===
    model_file = fullfile(model_folder, ['AveragedPLSModels_LSDvsPlacebo' ncomp_tag '.mat']);
    M = load(model_file);   % into a struct, so nothing carries over from the previous ncomp
    avg_BETA_final = M.avg_BETA_final;
    avg_B_cov_all  = M.avg_B_cov_all;
    avg_x_mean_all = M.avg_x_mean_all;
    avg_x_std_all  = M.avg_x_std_all;
    if isfield(M, 'max_ncomp')
        assert(M.max_ncomp == ncomp, '%s says %d components, expected %d', model_file, M.max_ncomp, ncomp);
    end
    clear M
    fprintf('Loaded averaged Basel LSD vs Placebo model (%d comp) from: %s\n', ncomp, model_file);

    % The Basel model must have been trained on the same edge sets, in the same order
    assert(numel(avg_BETA_final) == num_edge_sets, ...
        'Basel model has %d edge sets but %d are loaded — retrain Basel on the Schaefer100+Tian38 edge sets', ...
        numel(avg_BETA_final), num_edge_sets);
    for m = 1:num_edge_sets
        n_edges = size(all_net_edges{m}, 1);
        assert(size(avg_BETA_final{m}, 1) == n_edges + 1 && numel(avg_x_mean_all{m}) == n_edges, ...
            'Edge set %d (%s): Basel model has %d weights, expected %d (edges + intercept)', ...
            m, model_labels{m}, size(avg_BETA_final{m}, 1), n_edges + 1);
    end

    % Basel betas must line up with the covariate columns: rows of
    % avg_B_cov_all{inet} are [Intercept; b_Age; b_Male; b_FD], one column per edge
    assert(size(avg_B_cov_all{1}, 1) == size(cov_X, 2), ...
        'Basel betas have %d rows but Imperial covariates give %d columns', size(avg_B_cov_all{1}, 1), size(cov_X, 2));

    %% === Residualize Imperial FC Data ===
    % ---------------------------------------------------------------------
    % predicted part = Intercept + b_Age*Age + b_Male*Male + b_FD*FD
    % residual       = actual edge value - predicted part
    % Covariate VALUES always come from each Imperial scan; resid_mode only
    % changes where the BETAS come from:
    %   'basel'  = PRIMARY. Intercept AND slopes from Basel training data
    %              (avg_B_cov_all, averaged over all CV folds). "Frozen".
    %   'hybrid' = Basel slopes, intercept replaced by Imperial's own edge mean
    %              (removes any site-wide offset; Reviewer comment 1).
    %   'local'  = SENSITIVITY. Betas refitted on Imperial itself.
    % The residualization/z-scoring averages do not depend on ncomp (same
    % folds for every component count), so 'basel' residuals are identical
    % across ncomp; they are recomputed here only to keep each pass self-contained.
    % ---------------------------------------------------------------------
    X_resid_all = cell(num_edge_sets, 1);
    for inet = 1:num_edge_sets
        e = all_net_edges{inet};
        X_raw = NaN(num_subs, size(e, 1));
        for i = 1:size(e, 1)
            X_raw(:, i) = squeeze(sfc_mat(e(i,1), e(i,2), :));
        end

        switch resid_mode
            case 'basel'
                X_resid_all{inet} = X_raw - cov_X * avg_B_cov_all{inet};
            case 'hybrid'
                R = X_raw - cov_X * avg_B_cov_all{inet};
                X_resid_all{inet} = R - mean(R, 1);
            case 'local'
                X_resid_all{inet} = X_raw - cov_X * (cov_X \ X_raw);
            otherwise
                error('Unknown resid_mode: %s', resid_mode);
        end
    end
    fprintf('Residualized Imperial data (%s mode) for %d edge sets.\n', resid_mode, num_edge_sets);
    save(fullfile(output_dir_Imperial, ['Imperial_X_resid_all_' resid_mode ncomp_tag '.mat']), ...
        'X_resid_all', 'resid_mode', 'ncomp', '-v7.3');

    %% === Apply Basel-Trained Classifier to Imperial (both z-scoring modes) ===
    %   'basel_zscore'  = PRIMARY. Basel TRAINING mean/SD (frozen).
    %   'within_zscore' = SENSITIVITY. Imperial's own mean/SD across all Imperial scans
    %                     (labels not used); removes site-wide offset/scale.
    for zm = 1:numel(zscore_modes)
        zmode = zscore_modes{zm};
        fprintf('\n=========== %d COMP | Z-SCORING: %s ===========\n', ncomp, upper(zmode));

        acc_obs   = nan(num_edge_sets, 1);  p_acc    = nan(num_edge_sets, 1);
        auc_obs   = nan(num_edge_sets, 1);  p_auc    = nan(num_edge_sets, 1);
        dprime_obs= nan(num_edge_sets, 1);  p_dprime = nan(num_edge_sets, 1);
        crit_obs  = nan(num_edge_sets, 1);
        pair_diff = nan(num_edge_sets, 1);  p_pair   = nan(num_edge_sets, 1);
        prop_pos  = nan(num_edge_sets, 1);
        LSD_hit   = nan(num_edge_sets, 1);  Placebo_FA = nan(num_edge_sets, 1);
        predicted_probs_all  = cell(num_edge_sets, 1);   % continuous scores
        predicted_labels_all = cell(num_edge_sets, 1);   % 1 = LSD, 0 = Placebo

        for inet = 1:num_edge_sets

            % === 1. Residualized Imperial features [scans x edges] ===
            X = X_resid_all{inet};

            % === 2. Z-score: EXTERNAL (Basel) or WITHIN (Imperial) parameters ===
            switch zmode
                case 'basel_zscore'
                    mu    = avg_x_mean_all{inet};
                    sigma = avg_x_std_all{inet};
                case 'within_zscore'
                    mu    = mean(X, 1);
                    sigma = std(X, 0, 1);
            end
            sigma(sigma == 0) = 1;
            X_zscored = (X - mu) ./ sigma;

            % === 3. Predict with the Basel model for this ncomp ===
            y_pred  = [ones(size(X_zscored,1), 1), X_zscored] * avg_BETA_final{inet};
            y_class = double(y_pred >= 0.5);   % 1 = LSD, 0 = Placebo

            % === 4. Hits / false alarms vs ground truth ===
            hits = sum(y_class == 1 & label_vector == 1);   % LSD called LSD
            fas  = sum(y_class == 1 & label_vector == 0);   % Placebo called LSD
            acc  = mean(y_class == label_vector);

            % === 5. Observed metrics + permutation nulls ===
            null_acc = mean(y_class == perm_labels, 1);

            % AUC via ranks — threshold-free, not affected by a global offset
            r = tiedrank(y_pred);
            auc_fun  = @(L) (sum(r .* L, 1) - n1*(n1+1)/2) / (n1*n0);
            auc      = auc_fun(label_vector);
            null_auc = auc_fun(perm_labels);

            % d' (separability) and criterion c (bias)
            hit_z  = z_rate(hits, n1);
            fa_z   = z_rate(fas,  n0);
            dprime = hit_z - fa_z;
            crit   = -0.5 * (hit_z + fa_z);
            null_dp = z_rate(sum(y_class == 1 & perm_labels == 1, 1), n1) ...
                    - z_rate(sum(y_class == 1 & perm_labels == 0, 1), n0);

            % Within-subject: LSD score minus Placebo score per subject
            d_pair    = y_pred(1:n1) - y_pred(n1+1:end);
            null_diff = mean(d_pair .* (1 - 2*flips), 1);   % label swap = sign flip

            % === 6. Store ===
            LSD_hit(inet)    = hits;    Placebo_FA(inet) = fas;
            acc_obs(inet)    = acc;     p_acc(inet)    = p_fun(null_acc, acc);
            auc_obs(inet)    = auc;     p_auc(inet)    = p_fun(null_auc, auc);
            dprime_obs(inet) = dprime;  p_dprime(inet) = p_fun(null_dp, dprime);
            crit_obs(inet)   = crit;
            pair_diff(inet)  = mean(d_pair);
            prop_pos(inet)   = mean(d_pair > 0);
            p_pair(inet)     = p_fun(null_diff, mean(d_pair));
            predicted_probs_all{inet}  = y_pred;
            predicted_labels_all{inet} = y_class;

            % === 7. Print (one line per edge set) ===
            fprintf('[%dc %s] %03d %-30s Acc %5.1f%% | Hit %d/%d FA %d/%d | d''=%5.2f p=%.4f | AUC=%.2f p=%.4f | LSD>Plcb %3.0f%% p=%.4f\n', ...
                ncomp, zmode, inet, model_labels{inet}, 100*acc, hits, n1, fas, n0, ...
                dprime, p_dprime(inet), auc, p_auc(inet), 100*prop_pos(inet), p_pair(inet));
        end

        % === FDR (Benjamini-Hochberg) across edge sets, per metric ===
        p_list = {p_acc, p_auc, p_dprime, p_pair};
        q_list = cell(size(p_list));
        for k = 1:numel(p_list)
            pv = p_list{k}(:);  m = numel(pv);
            [ps, idx] = sort(pv);
            q = nan(m, 1);  q(idx) = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));
            q_list{k} = q;
        end

        % === Summary table for THIS ncomp x z-scoring mode ===
        T = table((1:num_edge_sets)', model_labels(:), repmat(ncomp, num_edge_sets, 1), ...
            repmat(string(resid_mode), num_edge_sets, 1), repmat(string(zmode), num_edge_sets, 1), ...
            LSD_hit, repmat(n1, num_edge_sets, 1), Placebo_FA, repmat(n0, num_edge_sets, 1), ...
            100*acc_obs, p_acc, q_list{1}, auc_obs, p_auc, q_list{2}, ...
            dprime_obs, p_dprime, q_list{3}, crit_obs, ...
            pair_diff, 100*prop_pos, p_pair, q_list{4}, repmat(n_perm, num_edge_sets, 1), ...
            'VariableNames', {'EdgeSet', 'Model', 'NComp', 'ResidMode', 'ZscoreMode', ...
            'LSD_Hits', 'LSD_N', 'Placebo_FalseAlarms', 'Placebo_N', ...
            'Accuracy_Pct', 'Acc_perm_p', 'Acc_perm_q', 'AUC', 'AUC_perm_p', 'AUC_perm_q', ...
            'Dprime', 'Dprime_perm_p', 'Dprime_perm_q', 'Criterion_c', ...
            'WithinSubj_MeanDiff_LSDminusPlcb', 'WithinSubj_Pct_LSDhigher', ...
            'WithinSubj_perm_p', 'WithinSubj_perm_q', 'N_permutations'});

        csv_path = fullfile(output_dir_Imperial, sprintf('BASEL_on_Imperial_resid-%s_%s%s_Summary.csv', resid_mode, zmode, ncomp_tag));
        writetable(T, csv_path);
        fprintf('Saved %d comp %s summary: %s\n', ncomp, zmode, csv_path);

        res.(comp_key).(zmode).table      = T;
        res.(comp_key).(zmode).pred_probs = predicted_probs_all;
        res.(comp_key).(zmode).pred_class = predicted_labels_all;
    end

    % Predicted labels/scores for this ncomp (for the FEP-like / HC-like and
    % violin scripts), primary z-scoring mode
    predicted_labels_all = res.(comp_key).basel_zscore.pred_class;
    predicted_probs_all  = res.(comp_key).basel_zscore.pred_probs;
    save(fullfile(output_dir_Imperial, ['BASEL_on_Imperial_resid-' resid_mode '_predictions' ncomp_tag '.mat']), ...
        'predicted_labels_all', 'predicted_probs_all', 'label_vector', 'model_labels', 'ncomp', 'resid_mode');
end

%% =========================================================================
% COMPARISON ACROSS COMPONENT COUNTS (one wide table per z-scoring mode)
% =========================================================================
ncomp_str = strjoin(arrayfun(@num2str, ncomp_list, 'UniformOutput', false), '-');
for zm = 1:numel(zscore_modes)
    zmode = zscore_modes{zm};
    C = table((1:num_edge_sets)', model_labels(:), 'VariableNames', {'EdgeSet', 'Model'});
    fprintf('\n=== %s: Basel model applied to Imperial, by number of components ===\n', zmode);
    for icn = 1:numel(ncomp_list)
        t  = sprintf('_%dcomp', ncomp_list(icn));
        Tk = res.(sprintf('comp%d', ncomp_list(icn))).(zmode).table;
        C.(['Accuracy_Pct' t])  = Tk.Accuracy_Pct;
        C.(['Dprime' t])        = Tk.Dprime;
        C.(['Dprime_perm_q' t]) = Tk.Dprime_perm_q;
        C.(['AUC' t])           = Tk.AUC;
        C.(['AUC_perm_q' t])    = Tk.AUC_perm_q;
        C.(['WithinSubj_perm_q' t]) = Tk.WithinSubj_perm_q;
        fprintf('%d comp: full-connectome d'' = %.2f, AUC = %.2f | d'' FDR q<.05: %d/%d | AUC q<.05: %d/%d\n', ...
            ncomp_list(icn), Tk.Dprime(1), Tk.AUC(1), sum(Tk.Dprime_perm_q < .05), num_edge_sets, ...
            sum(Tk.AUC_perm_q < .05), num_edge_sets);
    end
    cmp_path = fullfile(output_dir_Imperial, sprintf('BASEL_on_Imperial_resid-%s_%s_CompareNcomp_%scomp.csv', resid_mode, zmode, ncomp_str));
    writetable(C, cmp_path);
    fprintf('Saved comparison: %s\n', cmp_path);
end

save(fullfile(output_dir_Imperial, sprintf('BASEL_on_Imperial_resid-%s_AllNcomp_%scomp.mat', resid_mode, ncomp_str)), ...
    'res', 'ncomp_list', 'resid_mode', 'zscore_modes', 'flips', 'n_perm', 'model_labels', '-v7.3');
fprintf('\nDone. All results in: %s\n', output_dir_Imperial);