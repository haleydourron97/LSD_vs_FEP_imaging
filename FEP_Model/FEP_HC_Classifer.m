%% ============================================================
% FEP vs HC Classification using PLS-DA 
% ------------------------------------------------------------
%   1)LOADING DATA: This script loads functional connectivity (FC) matrices and
%   demographic covariates for first-episode psychosis (FEP) and
%   healthy control (HC) participants. 
%    
%   2) LOADING EDGES & LABELS: The edges and labels for the full connectome, between-network, and
%   within-network are loaded. 
%  
%   3)CLASSIFER: The partial least squares discriminant classifer is run to classify
%   FEP or HC. The averaged model is saved. The model is ran with 1, 2, 3
%   components as a sensitivity test for 
%
%   4)PERMUATION TESTING: Permuatation testing is done by shuffling the labels of FEP and HC 
%    5,000 times to create a null distribution. 
%   
%   5)P-Values & Correction: We then compare the results with and without false discovery rate
%   correction, with permutation p-values generate from the permuatation-based
%   null distribution. 

%% 1) Loading DATA %%

% === Define Directories ===
FEP_HC_dir    = '/Volumes/toshiba/LSD_November11/FEPModels/FEP_HC_Data_fMRI/NILEARN_Processed_Pearson_Z_transformed_w_Motion/combinedsubcorticalOnly_NO_GSR';
FEP_dir       = '/Volumes/toshiba/LSD_November11/FEPModels/FEP_Data_fMRI/NILEARN_Processed_Pearson_Z_transformed_w_Motion/combinedsubcorticalOnly_NO_GSR';
covariate_input_dir = '/Volumes/toshiba/LSD_November11/FEPModels/Covariates';
output_dir    = '/Volumes/toshiba/LSD_November11/FEPModels/Output';


fc_matrices = {};             % Combined FEP + HC matrices
subject_ids = strings(144,1); % Store subject ID order: 72 FEP + 72 HC % changed because one FEP had NAs for some subcortical regions

%% === Load FEP Covariates ===
cov_table_FEP = readtable(fullfile(covariate_input_dir, 'FEP.csv'));

% Get FEP file list
FEP_files = dir(fullfile(FEP_dir, 'sub-*_combinedsubcorticalOnly_TS_denoised_NO_GSR_*.csv'));

% Extract numeric subject IDs from filenames
all_subnums = arrayfun(@(f) str2double(regexp(f.name, '(?<=^sub-)\d+', 'match', 'once')), FEP_files);

% Keep only subjects listed in the covariate table (it already contains only subjects who meet the inclusion criteria)
valid_idx   = ismember(all_subnums, cov_table_FEP.subid);
FEP_files   = FEP_files(valid_idx);
all_subnums = all_subnums(valid_idx);

n_files = length(FEP_files);

% Log selected filenames
log_file_FEP = fullfile(covariate_input_dir, 'selected_FEP.txt');
fid = fopen(log_file_FEP, 'w');

selected_ids_FEP = zeros(n_files, 1); % numeric subid values

for i = 1:n_files
    file_path = fullfile(FEP_dir, FEP_files(i).name);
    fprintf(fid, '%s\n', FEP_files(i).name);

    matrix_data = readmatrix(file_path);
    fc_matrices{end+1} = matrix_data;

    selected_ids_FEP(i) = all_subnums(i);          % numeric ID (e.g., 26)
    subject_ids(i) = "sub-" + all_subnums(i);      % string ID (e.g., "sub-26")
end

fclose(fid);

fprintf('Loaded %d FEP files (subjects in covariate table) and saved list to %s\n', n_files, log_file_FEP);

% Match subid from file to subid in covariate table
[is_match, cov_idx] = ismember(selected_ids_FEP, cov_table_FEP.subid);

% Extract matched covariates
Age   = cov_table_FEP.Age(cov_idx);
Male  = cov_table_FEP.Male(cov_idx);
FD    = cov_table_FEP.Mean_framewise_displacement(cov_idx);

% Combine covariates into a matrix
covariates_FEP = [Age, Male, FD];

% Define number of FEP subjects for use in HC block
num_FEP = n_files;

%% === Load HC Covariates ===
cov_table_HC = readtable(fullfile(covariate_input_dir, 'FEP_HC.csv'));
% Get HC file list
FEP_HC_files = dir(fullfile(FEP_HC_dir, 'sub-*_combinedsubcorticalOnly_TS_denoised_NO_GSR_*.csv'));
% Extract 'sub-XX' prefixes from filenames
file_prefixes = arrayfun(@(f) regexp(f.name, 'sub-\d{2,3}', 'match', 'once'), ...
    FEP_HC_files, 'UniformOutput', false);
file_prefixes = string(file_prefixes);
% Keep only subjects listed in the covariate table (it already contains only subjects who meet the inclusion criteria)
valid_idx     = ismember(double(extractAfter(file_prefixes, "sub-")), cov_table_HC.subid);
FEP_HC_files  = FEP_HC_files(valid_idx);
file_prefixes = file_prefixes(valid_idx);

% Randomly select 72 from remaining
rng(42);
selected_idx_HC      = randperm(length(FEP_HC_files), 72);
selected_files_HC    = FEP_HC_files(selected_idx_HC);
selected_prefixes_HC = file_prefixes(selected_idx_HC);

% Log + load
log_file_HC = fullfile(covariate_input_dir, 'selected_FEP_HC.txt');
fid = fopen(log_file_HC, 'w');
for i = 1:72
    file_path = fullfile(FEP_HC_dir, selected_files_HC(i).name);
    fprintf(fid, '%s\n', selected_files_HC(i).name);
    matrix_data = readmatrix(file_path);
    fc_matrices{end+1} = matrix_data;
    % Insert into preallocated subject_ids array (kept as "sub-XX" string)
    subject_ids(num_FEP + i) = selected_prefixes_HC(i);
end
fclose(fid);
fprintf('Loaded 72 HC files (subjects in covariate table) and saved list to %s\n', log_file_HC);

num_HC   = numel(selected_files_HC); % Number of HCs

% -----------------------------------------------------------------
% Match to covariates using numeric subid.
% Filenames still use "sub-XX" format, so extract the numeric part
% from selected_prefixes_HC before matching against cov_table_HC.subid.
% -----------------------------------------------------------------
selected_ids_HC = double(extractAfter(selected_prefixes_HC, "sub-"));

[is_match, cov_idx] = ismember(selected_ids_HC, cov_table_HC.subid);

% Extract covariates
Age   = cov_table_HC.Age(cov_idx);
Male  = cov_table_HC.Male(cov_idx);
FD    = cov_table_HC.Mean_framewise_displacement(cov_idx);

% Combine into a matrix
covariates_HC = [Age, Male, FD];

%% === Define Binary Class Labels for FEP and HC ===
% Total subject count
num_subs = num_FEP + num_HC;

% Sanity check: label count must match loaded FC matrices
if num_subs ~= length(fc_matrices)
    error('Mismatch in subject count: expected %d, found %d matrices', num_subs, length(fc_matrices));
end

% Create binary class label vector: 1 = FEP, 0 = HC
label_vector = [ ...
    ones(num_FEP, 1); ...  % FEP = 1
    zeros(num_HC, 1) ...   % HC  = 0
];
fprintf('Created binary label vector [%d subjects] (FEP=1, HC=0)\n', num_subs);
fprintf('   FEP: %d | HC: %d\n', num_FEP, num_HC);

%% === Combine covariates across groups ===
covariates_all = [covariates_FEP; covariates_HC];
age_all     = covariates_all(:, 1);
male_all    = covariates_all(:, 2);
mean_fd_all = covariates_all(:, 3);

% Numeric subject IDs in the same order (FEP first, then HC)
subject_id_vector = [selected_ids_FEP; selected_ids_HC];

%% === Additional sanity checks and final save ===
if num_subs ~= length(subject_ids) || num_subs ~= length(age_all) ...
        || num_subs ~= length(male_all) ...
        || num_subs ~= length(mean_fd_all) || num_subs ~= length(label_vector) ...
        || num_subs ~= length(subject_id_vector)
    error('Mismatch in number of subjects across matrices, covariates, labels, or IDs');
end
if any(isnan(age_all)) || any(isnan(male_all)) || any(isnan(mean_fd_all))
    warning('NaNs present in covariates — check for missing subjects before residualizing');
end
% Build final FC array: ROI x ROI x subjects
sfc_mat = cat(3, fc_matrices{:});
num_rois = size(sfc_mat, 1);
if ~isequal(size(sfc_mat), [138 138 num_subs])
    warning('Unexpected sfc_mat size: %s (expected [138 138 %d])', mat2str(size(sfc_mat)), num_subs);
end
% Build covariates matrix: columns = [Age, Male, Mean_FD]
covariates = [age_all, male_all, mean_fd_all];
% Save outputs
tag = 'FEP_DATASETS_combinedsubcorticalOnly_NO_GSR';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
save(fullfile(output_dir, [tag '_sfc_mat_FEP_vs_HC.mat']), 'sfc_mat');
save(fullfile(output_dir, [tag '_covariates_FEP_vs_HC.mat']), 'covariates');
save(fullfile(output_dir, [tag '_label_vector_FEP_vs_HC.mat']), 'label_vector');
save(fullfile(output_dir, [tag '_subject_id_vector_FEP_vs_HC.mat']), 'subject_id_vector');
writematrix(subject_ids, fullfile(output_dir, [tag '_subject_ids_FEP_vs_HC.txt']));
fprintf('Saved sFC matrix, full covariates (Age/Male/FD), labels, and subject IDs for FEP vs HC classification.\n');

%% === Define Paths for PLS-DA Classifer & Edges ===
edge_dir       = '/Volumes/toshiba/LSD_November11/Within_Between_allEdge/Schaefer100_Tian38';
edge_file      = fullfile(edge_dir, 'edge_sets_and_labels.mat');

% === Edge Sets + Model Labels  ===
if ~exist(edge_file, 'file'), error('Missing file: %s', edge_file); end
E = load(edge_file, 'edge_sets', 'model_labels', 'model_type', ...
                    'model_i', 'model_j', 'group_labels', 'roi_group');

all_net_edges = E.edge_sets;          % 326 x 1: full, 25 within, 300 between
model_labels  = E.model_labels;       % 326 x 1, same order as all_net_edges
model_type    = E.model_type;         % 'full' / 'within' / 'between'
group_labels  = E.group_labels;       % the 25 group names
num_nets      = numel(group_labels);

%% === Sanity Checks ===
% 1. Subjects: matrices, labels, covariates and IDs all line up
assert(num_subs == numel(label_vector) && num_subs == size(covariates, 1) ...
    && num_subs == numel(subject_id_vector), ...
    'Subject count mismatch: sfc_mat %d, labels %d, covariates %d, IDs %d', ...
    num_subs, numel(label_vector), size(covariates, 1), numel(subject_id_vector));
assert(all(ismember(label_vector, [0 1])), 'label_vector must contain only 0 (HC) and 1 (FEP)');
assert(isequal(label_vector(:), [ones(num_FEP,1); zeros(num_HC,1)]), ...
    'label_vector is not ordered FEP first, then HC (expected by the loading block)');
if any(isnan(covariates(:)))
    warning('NaNs present in covariates - check for missing subjects before residualizing');
end

% 2. Atlas matches the FC matrices
assert(num_rois == numel(E.roi_group), ...
    'sfc_mat is %dx%d but edge sets are for %d ROIs (Schaefer100+Tian38 expected)', ...
    num_rois, num_rois, numel(E.roi_group));

% 3. One label per edge set, and no edge index out of range
assert(numel(all_net_edges) == numel(model_labels), 'Mismatch between model labels and edge sets!');
assert(max(cellfun(@(e) max(e(:)), all_net_edges)) <= num_rois, 'Edge index exceeds number of ROIs');

% 4. Every set contains only edges from the groups its label names
%    (this is the check that would have caught the old ordering bug)
for m = 2:numel(all_net_edges)
    e    = all_net_edges{m};
    pair = sort([E.roi_group(e(:,1)), E.roi_group(e(:,2))], 2);
    assert(all(pair(:,1) == E.model_i(m) & pair(:,2) == E.model_j(m)), ...
        'Edge set %d (%s) contains edges from the wrong groups', m, model_labels{m});
end

% 5. Edge values used as features are finite (diagonal / bad values excluded)
full_e = all_net_edges{1};
lin    = sub2ind([num_rois num_rois], full_e(:,1), full_e(:,2));
vals   = reshape(sfc_mat, num_rois*num_rois, num_subs);
vals   = vals(lin, :);
assert(all(isfinite(vals(:))), 'Non-finite FC values found in the edges used as features');

fprintf('Loaded sfc_mat [%d x %d x %d]: FEP = %d | HC = %d\n', num_rois, num_rois, num_subs, num_FEP, num_HC);
fprintf('Loaded %d edge sets (%d full, %d within, %d between) for %d ROIs\n', ...
    numel(all_net_edges), sum(strcmp(model_type,'full')), sum(strcmp(model_type,'within')), ...
    sum(strcmp(model_type,'between')), num_rois);
fprintf('%d model labels match %d edge sets (verified edge-by-edge).\n', ...
    numel(model_labels), numel(all_net_edges));

%% === Run Cross-Validation Binary PLS Classifier (FEP vs HC, Between-Subject) ===
% ===============================================================
% FEP vs HC: Binary PLS Classifier with FOLD-WISE Covariate Control
% BETWEEN-SUBJECT design — each subject contributes exactly ONE scan
% (either FEP or HC), so a row-wise K-fold partition stratified by
% group is appropriate.
% Run once each for 1, 2 and 3 PLS components (ncomp_list).
% ===============================================================
%
%   Here every row is a different person, so rows can be partitioned
%   directly. Folds are stratified by GROUP (FEP/HC) so each fold keeps
%   roughly the same FEP:HC ratio as the full sample.
%
% Structure matches the LSD vs Placebo script:
%   - The partition is built ONCE per repeat (rng(r)) and shared by all
%     edge sets, so every model in a repeat uses exactly the same splits
%     (differences between edge sets reflect the features, not the split).
%     The partitions are also identical across component counts.
%   - Fold-level models (BETA, B_cov, x_mean, x_std) are stored on every
%     fold of every repeat and averaged across all K*num_repeats folds at
%     the end to produce ONE averaged model per edge set. No model is ever
%     fit on 100% of the data.
%   - Fold-wise covariate residualization, fold-wise z-scoring and PLS
%     training are identical to the LSD vs Placebo script.

% Configuration
num_repeats = 100;
K = 10;
ncomp_list = [1 2 3];    % comparing results by number of components
verbose_folds = false;   % true = print every fold (326 sets x 100 repeats x 10 folds = 326,000 lines per ncomp)

% output_dir is defined in the loading block (FEPModels/Output)
model_labels = cellstr(model_labels);

% Number of edge/network models
num_edge_sets = length(all_net_edges);
assert(numel(model_labels) == num_edge_sets, ...
    'model_labels (%d) and all_net_edges (%d) differ in length', numel(model_labels), num_edge_sets);
num_subs = length(label_vector); % total rows (one scan per subject)
assert(num_subs == size(sfc_mat, 3) && num_subs == size(covariates, 1) && num_subs == numel(subject_id_vector), ...
    'label_vector (%d), sfc_mat (%d), covariates (%d) and subject_id_vector (%d) differ', ...
    num_subs, size(sfc_mat, 3), size(covariates, 1), numel(subject_id_vector));
fprintf('%d subjects (between-subject): FEP=%d | HC=%d\n', ...
    num_subs, sum(label_vector == 1), sum(label_vector == 0));

% Collected across component counts for the comparison table at the end
dprime_by_ncomp  = NaN(num_edge_sets, numel(ncomp_list));
acc_by_ncomp     = NaN(num_edge_sets, numel(ncomp_list));
dprime_sd_ncomp  = NaN(num_edge_sets, numel(ncomp_list));

% ===============================================================
% LOOP OVER NUMBER OF PLS COMPONENTS
% Outputs are tagged _1comp / _2comp / _3comp so nothing is overwritten.
% ===============================================================
for ic = 1:numel(ncomp_list)
    max_ncomp = ncomp_list(ic);
    ncomp_tag = sprintf('_%dcomp', max_ncomp);
    fprintf('\n===== Running with max %d PLS component(s) =====\n', max_ncomp);
    t_ncomp = tic;

    % Initialize performance containers
    train_accuracies_all = cell(num_edge_sets, 1);
    test_accuracies_all  = cell(num_edge_sets, 1);
    per_class_acc_all    = cell(num_edge_sets, 1);

    % Subject-level held-out predictions across repeats
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
        fprintf('[%dcomp] Repetition %d of %d...\n', max_ncomp, r, num_repeats);
        rng(r);   % reproducible partitions per repeat (identical across ncomp)

        % -------------------------------------------------------
        % GROUP-STRATIFIED SUBJECT-LEVEL K-FOLD PARTITION
        %
        % Built ONCE per repeat and shared by all edge sets. One row = one
        % subject, so partitioning rows = partitioning subjects. Stratified
        % by GROUP so each fold has ~the same FEP:HC ratio.
        % -------------------------------------------------------
        cv_partition = cvpartition(label_vector, 'KFold', K, 'Stratify', true);

        % Row indices per fold (same for every edge set in this repeat)
        train_idx_all = cell(K, 1);
        test_idx_all  = cell(K, 1);
        for k = 1:K
            train_idx_all{k} = training(cv_partition, k);
            test_idx_all{k}  = test(cv_partition, k);

            % Sanity check: no subject in both train and test
            if any(train_idx_all{k} & test_idx_all{k})
                error('Subject leakage detected between train/test in fold %d', k);
            end
        end

        % Loop over each edge set / network model
        for inet = 1:num_edge_sets

            % -----------------------------------------------------------
            % Build feature matrix for this edge set.
            %   rows    = subjects
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
            per_fold_class_acc = zeros(K, 2); % columns: [HC, FEP]

            % Exactly one held-out prediction per subject for this repeat
            heldout_class_this_repeat = NaN(num_subs, 1);
            heldout_score_this_repeat = NaN(num_subs, 1);

            % ===========================================================
            % loop over CV folds
            % ===========================================================
            for k = 1:K

                train_idx = train_idx_all{k};
                test_idx  = test_idx_all{k};

                % -------------------------------------------------------
                % Split connectivity features, labels and covariates
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
                % Betas fit on TRAINING subjects only, applied to train and test.
                % -------------------------------------------------------
                B_cov = C_train \ X_train_raw;
                x_train = X_train_raw - C_train * B_cov;
                x_test  = X_test_raw  - C_test  * B_cov;

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
                % (capped for edge sets with fewer edges than max_ncomp)
                % -------------------------------------------------------
                ncomp = min([size(x_train, 1)-1, size(x_train, 2), max_ncomp]);
                [~, ~, ~, ~, BETA] = plsregress(x_train, y_train, ncomp);

                % -------------------------------------------------------
                % Predict held-out TEST subjects
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

                if verbose_folds, fprintf(' [%dcomp] Edge %d (%s) | Fold %d — TrainAcc: %.3f | TestAcc: %.3f | HC: %.3f, FEP: %.3f\n', ...
                    max_ncomp, inet, model_labels{inet}, k, fold_train_acc(k), fold_test_acc(k), ...
                    per_fold_class_acc(k,1), per_fold_class_acc(k,2)); end

                % -------------------------------------------------------
                % Store this fold's model for later averaging
                % -------------------------------------------------------
                fold_BETA_all{inet} = cat(3, fold_BETA_all{inet}, BETA);   % [p+1 x 1 x n_folds_total]
                fold_Bcov_all{inet} = cat(3, fold_Bcov_all{inet}, B_cov);  % [n_cov+1 x p x n_folds_total]
                fold_mean_all{inet} = [fold_mean_all{inet}; x_mean];       % [n_folds_total x p]
                fold_std_all{inet}  = [fold_std_all{inet};  x_std];        % [n_folds_total x p]
            end

            % Sanity check: every subject should have one held-out prediction
            if any(isnan(heldout_class_this_repeat)) || any(isnan(heldout_score_this_repeat))
                warning('Missing held-out predictions in Repeat %d, Edge Set %d (%dcomp)', r, inet, max_ncomp);
            end

            % Store this repeat's subject-level predictions
            heldout_class_all{inet}(:, r) = heldout_class_this_repeat;
            heldout_score_all{inet}(:, r) = heldout_score_this_repeat;

            % Append fold-wise performance
            train_accuracies_all{inet} = [train_accuracies_all{inet}; fold_train_acc];
            test_accuracies_all{inet}  = [test_accuracies_all{inet}; fold_test_acc];
            per_class_acc_all{inet}    = [per_class_acc_all{inet}; per_fold_class_acc];
        end
    end

    % === Check the CV actually finished before averaging ===
    n_models = cellfun(@(b) size(b, 3), fold_BETA_all);
    assert(all(n_models == K * num_repeats), ...
        'CV did not finish (%dcomp): %d of %d edge sets have all %d fold models (set 1 has %d).', ...
        max_ncomp, sum(n_models == K*num_repeats), numel(n_models), K*num_repeats, n_models(1));
    fprintf('CV complete (%dcomp): %d edge sets x %d fold-level models each.\n', ...
        max_ncomp, num_edge_sets, K*num_repeats);


    % ===============================================================
    % Average fold-level models across ALL folds and repeats
    % (K * num_repeats models per edge set)
    %
    % Averaging strictly train-only, fold-level models keeps this consistent
    % with the LSD/Placebo pipeline, and means the saved model was never fit
    % on 100% of the data.
    % ===============================================================
    avg_BETA_final = cell(num_edge_sets, 1);
    avg_B_cov_all  = cell(num_edge_sets, 1);  % [Int, Age, Male, FD] x n_edges
    avg_x_mean_all = cell(num_edge_sets, 1);
    avg_x_std_all  = cell(num_edge_sets, 1);

    for inet = 1:num_edge_sets
        avg_BETA_final{inet} = mean(fold_BETA_all{inet}, 3);
        avg_B_cov_all{inet}  = mean(fold_Bcov_all{inet}, 3);
        avg_x_mean_all{inet} = mean(fold_mean_all{inet}, 1);
        % Pool SDs via variance rather than averaging SDs directly
        avg_x_std_all{inet}  = sqrt(mean(fold_std_all{inet}.^2, 1));
    end
    fprintf('Averaged models computed for %d edge sets (from %d fold-level models each)\n', ...
        num_edge_sets, size(fold_BETA_all{1}, 3));

    model_path = fullfile(output_dir, ['AveragedPLSModels_FEPvsHC' ncomp_tag '.mat']);
    save(model_path, 'avg_BETA_final', 'avg_B_cov_all', 'avg_x_mean_all', 'avg_x_std_all', ...
        'model_labels', 'max_ncomp');
    fprintf('Saved averaged models to: %s\n', model_path);


    % === Save all subject-level held-out predictions (FEP vs HC) ===
    heldout_save_path = fullfile(output_dir, ['FEP_HeldOutPredictions_AllEdgeSets_FEPvsHC' ncomp_tag '.mat']);
    save(heldout_save_path, 'heldout_class_all', 'heldout_score_all', ...
        'label_vector', 'subject_id_vector', 'model_labels', 'max_ncomp', '-v7.3');
    fprintf('Saved subject-level held-out predictions to: %s\n', heldout_save_path);


    % === Subject-level summary CSV: average across repeats ===
    % subject_id_vector is NUMERIC here (e.g. 26), so build "sub-26" strings
    % (cellstr on a numeric vector would error).
    subject_id_cell = cellstr("sub-" + string(subject_id_vector(:)));

    summary_subjects   = {};
    summary_labels     = [];
    summary_edge_idx   = [];
    summary_edge_names = {};
    mean_pred_class    = [];
    mean_pred_score    = [];
    sd_pred_score      = [];
    prop_fep_class     = [];

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
        prop_fep_class  = [prop_fep_class; mean(class_mat == 1, 2, 'omitnan')];   % 1 = FEP
    end

    summary_groups = repmat({''}, numel(summary_labels), 1);
    summary_groups(summary_labels == 0) = {'HC'};
    summary_groups(summary_labels == 1) = {'FEP'};

    T_summary = table(summary_subjects, summary_labels, summary_groups, summary_edge_idx, ...
        summary_edge_names, mean_pred_class, prop_fep_class, mean_pred_score, sd_pred_score, ...
        'VariableNames', {'SubjectID', 'TrueLabel', 'TrueGroup', 'EdgeSetIndex', 'EdgeSetLabel', ...
                          'MeanHeldOutPredClass', 'PropRepeatsClassifiedFEP', ...
                          'MeanHeldOutPredScore', 'SDHeldOutPredScore'});

    summary_csv_path = fullfile(output_dir, ['FEP_HeldOutPredictions_SubjectSummary_FEPvsHC' ncomp_tag '.csv']);
    writetable(T_summary, summary_csv_path);
    fprintf('Saved subject-summary CSV: %s\n', summary_csv_path);


    % =========================================================================
    % PER-ITERATION METRICS (one row per model x repeat) — FEP vs HC
    % -------------------------------------------------------------------------
    % Collapses each repeat's K folds into one repeat-level value; d' per repeat
    % from the subject-level held-out predictions.
    % per_class_acc_all columns are [HC, FEP] (class 0, class 1).
    % Hit = FEP called FEP; False alarm = HC called FEP.
    % =========================================================================
    num_models = numel(model_labels);
    n_repeats_observed = numel(test_accuracies_all{1}) / K;
    n_iter_rows = num_models * n_repeats_observed;

    Model   = cell(n_iter_rows, 1);
    EdgeSet = NaN(n_iter_rows, 1);
    Repeat  = NaN(n_iter_rows, 1);
    TestAcc = NaN(n_iter_rows, 1);
    HC_Acc  = NaN(n_iter_rows, 1);
    FEP_Acc = NaN(n_iter_rows, 1);
    Dprime  = NaN(n_iter_rows, 1);

    row_idx = 1;
    for inet = 1:num_models
        % --- collapse K-fold accuracy down to one value per repeat ---
        test_vals = test_accuracies_all{inet};             % [K*num_repeats x 1]
        n_rep     = numel(test_vals) / K;
        test_repeat_means = mean(reshape(test_vals, K, n_rep), 1)';

        class_vals = per_class_acc_all{inet};              % [K*num_repeats x 2]
        hc_repeat_means  = mean(reshape(class_vals(:,1), K, n_rep), 1)';
        fep_repeat_means = mean(reshape(class_vals(:,2), K, n_rep), 1)';

        % --- d' per repeat, from subject-level held-out predictions ---
        pred_mat = heldout_class_all{inet};                % [num_subs x num_repeats]
        dprime_vals = NaN(n_rep, 1);
        for r = 1:n_rep
            pred = pred_mat(:, r);
            hits   = sum(pred == 1 & label_vector == 1);   % FEP called FEP
            misses = sum(pred == 0 & label_vector == 1);   % FEP called HC
            fas    = sum(pred == 1 & label_vector == 0);   % HC called FEP
            crs    = sum(pred == 0 & label_vector == 0);   % HC called HC

            % log-linear correction avoids +/-Inf when hit/FA rate = 0 or 1
            hit_rate = (hits + 0.5) / (hits + misses + 1);
            fa_rate  = (fas  + 0.5) / (fas + crs + 1);
            dprime_vals(r) = norminv(hit_rate) - norminv(fa_rate);
        end

        rows = row_idx:(row_idx + n_rep - 1);
        Model(rows)   = repmat(model_labels(inet), n_rep, 1);
        EdgeSet(rows) = inet;
        Repeat(rows)  = (1:n_rep)';
        TestAcc(rows) = test_repeat_means;
        HC_Acc(rows)  = hc_repeat_means;
        FEP_Acc(rows) = fep_repeat_means;
        Dprime(rows)  = dprime_vals;
        row_idx = row_idx + n_rep;
    end

    per_iteration_summary = table(EdgeSet, Model, Repeat, TestAcc, HC_Acc, FEP_Acc, Dprime);
    per_iteration_summary{:, {'TestAcc','HC_Acc','FEP_Acc'}} = ...
        per_iteration_summary{:, {'TestAcc','HC_Acc','FEP_Acc'}} * 100;

    per_iteration_path = fullfile(output_dir, ['FEP_FinalSummary_FEPvsHC_PerIteration_WithDprime' ncomp_tag '.csv']);
    writetable(per_iteration_summary, per_iteration_path);
    fprintf('Per-iteration summary (%d rows) saved to: %s\n', n_iter_rows, per_iteration_path);


    % =========================================================================
    % PER-MODEL SUMMARY (mean / SD / 95% CI across repeats) — FEP vs HC
    % -------------------------------------------------------------------------
    % CI is descriptive only: it reflects variability across CV partitions,
    % not independent-sample uncertainty (repeats are not independent).
    % =========================================================================
    metrics = {'TestAcc', 'HC_Acc', 'FEP_Acc', 'Dprime'};

    [~, first_idx] = unique(per_iteration_summary.EdgeSet, 'stable');
    EdgeSet_s = per_iteration_summary.EdgeSet(first_idx);
    Model_s   = per_iteration_summary.Model(first_idx);
    num_models_summary = numel(EdgeSet_s);

    clear mean_col sd_col ci_lo ci_hi
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
        mean_col.TestAcc, sd_col.TestAcc, ci_lo.TestAcc, ci_hi.TestAcc, ...
        mean_col.HC_Acc,  sd_col.HC_Acc,  ci_lo.HC_Acc,  ci_hi.HC_Acc, ...
        mean_col.FEP_Acc, sd_col.FEP_Acc, ci_lo.FEP_Acc, ci_hi.FEP_Acc, ...
        mean_col.Dprime,  sd_col.Dprime,  ci_lo.Dprime,  ci_hi.Dprime, ...
        'VariableNames', { ...
            'EdgeSet', 'Model', ...
            'TestAcc_Mean', 'TestAcc_SD', 'TestAcc_CI_Lo', 'TestAcc_CI_Hi', ...
            'HC_Acc_Mean',  'HC_Acc_SD',  'HC_Acc_CI_Lo',  'HC_Acc_CI_Hi', ...
            'FEP_Acc_Mean', 'FEP_Acc_SD', 'FEP_Acc_CI_Lo', 'FEP_Acc_CI_Hi', ...
            'Dprime_Mean',  'Dprime_SD',  'Dprime_CI_Lo',  'Dprime_CI_Hi'});

    per_model_path = fullfile(output_dir, ['FEP_FinalSummary_FEPvsHC_PerModel_MeanAcrossRepeats' ncomp_tag '.csv']);
    writetable(per_model_summary, per_model_path);
    fprintf('Per-model summary (mean/SD/CI across %d repeats) saved to: %s\n', ...
        n_repeats_observed, per_model_path);

    % --- Keep this component count's results for the comparison table ---
    dprime_by_ncomp(:, ic) = per_model_summary.Dprime_Mean;
    dprime_sd_ncomp(:, ic) = per_model_summary.Dprime_SD;
    acc_by_ncomp(:, ic)    = per_model_summary.TestAcc_Mean;

    fprintf('===== Finished %d component(s) in %.1f h =====\n', max_ncomp, toc(t_ncomp)/3600);
end   % component loop

%%

%% =========================================================================
% PERMUTATION TEST — FEP vs HC (between-subject), 1, 2 AND 3 PLS COMPONENTS
% -------------------------------------------------------------------------
%
% Design:
%   - ONE fixed 10-fold partition (cv_seed). One row = one subject, so
%     rows are partitioned directly, stratified by GROUP (FEP/HC).
%   - Fold-wise residualization + z-scoring do not use labels, so they are
%     computed ONCE per edge set and reused for every permutation.
%   - Permutation 0 = true labels; 1..n_perm = FEP/HC labels permuted
%      across subjects, the SAME permutations for every edge set AND every component count.
%   - Each PLS fit uses 3 components; the 1- and 2-component models are
%     rebuilt from its first components (PLS components are sequential, so
%     this is identical to fitting 1 or 2 components directly). All three
%     component counts therefore cost about the same as one.
%
% p-values: one-sided, (#null >= observed + 1) / (n_perm + 1).
% FDR: BH across the 326 edge sets, SEPARATELY for each component count.
%
%
% OUTPUT (in <output_dir>/PermutationTest):
%   FEP_PermutationTest_FEPvsHC_Summary_1comp.csv  
%   FEP_PermutationTest_FEPvsHC_Summary_2comp.csv
%   FEP_PermutationTest_FEPvsHC_Summary_3comp.csv
%   FEP_PermutationTest_FEPvsHC_ByNComp_1-2-3comp.csv         (long, NComp column)
%   FEP_PermutationTest_FEPvsHC_NCompComparison_1-2-3comp.csv (wide)
%   FEP_PermutationTest_FEPvsHC_1-2-3comp.mat                 (incl. full null distributions)
% =========================================================================
n_perm       = 5000;
K            = 10;
cv_seed      = 2026;     % fixed partition
perm_seed    = 42;       % fixed label permutations
ncomp_list   = [1 2 3];
n_nc         = numel(ncomp_list);
ncomp_str    = strjoin(arrayfun(@num2str, ncomp_list, 'UniformOutput', false), '-');   % '1-2-3'
model_labels = cellstr(model_labels);
perm_output_dir = '/Users/haldo67/Desktop/FEP_LSD_Classifers/FEP_PermutationOutput';   
if ~exist(perm_output_dir, 'dir'), mkdir(perm_output_dir); end
checkpoint_path = fullfile(perm_output_dir, ['FEP_PermutationTest_FEPvsHC_checkpoint_' ncomp_str 'comp.mat']);

num_edge_sets = numel(all_net_edges);
num_subs      = numel(label_vector);
fprintf('%d subjects (between-subject): FEP=%d | HC=%d\n', ...
    num_subs, sum(label_vector == 1), sum(label_vector == 0));

%% === Fixed GROUP-stratified partition (one row = one subject) ===
rng(cv_seed);
cv_partition  = cvpartition(label_vector, 'KFold', K, 'Stratify', true);
train_idx_all = cell(K, 1);
test_idx_all  = cell(K, 1);
for k = 1:K
    train_idx_all{k} = training(cv_partition, k);
    test_idx_all{k}  = test(cv_partition, k);
end

%% === Fixed between-subject label permutations ===
rng(perm_seed);
perm_idx = zeros(num_subs, n_perm);
for ip = 1:n_perm
    perm_idx(:, ip) = randperm(num_subs)';                       % free permutation across subjects
end
Y_all = [label_vector(:), label_vector(perm_idx)];               % [subjects x (1+n_perm)], col 1 = true labels

%% === Containers: [edge set x permutation x ncomp] ===
if exist(checkpoint_path, 'file')
    C = load(checkpoint_path);
    % Refuse to resume from a checkpoint made with different subjects / settings
    if ~isfield(C, 'ckpt_subjects') || ~isequal(string(C.ckpt_subjects(:)), string(subject_id_vector(:))) ...
            || ~isequal(C.ckpt_ncomp, ncomp_list) || C.ckpt_nperm ~= n_perm
        error(['Checkpoint %s was made with different subjects or settings.\n' ...
               'Delete or rename it to start fresh.'], checkpoint_path);
    end
    obs_acc = C.obs_acc;  obs_dprime = C.obs_dprime;
    null_acc = C.null_acc;  null_dprime = C.null_dprime;  done = C.done;
    clear C
    fprintf('Resuming from checkpoint: %d of %d edge sets already done.\n', sum(done), num_edge_sets);
else
    obs_acc     = nan(num_edge_sets, n_nc);
    obs_dprime  = nan(num_edge_sets, n_nc);
    null_acc    = nan(num_edge_sets, n_perm, n_nc, 'single');
    null_dprime = nan(num_edge_sets, n_perm, n_nc, 'single');
    done        = false(num_edge_sets, 1);
end
ckpt_subjects = subject_id_vector;  ckpt_ncomp = ncomp_list;  ckpt_nperm = n_perm;

z_rate = @(k, n) norminv((k + 0.5) ./ (n + 1));   % log-linear correction
n1 = sum(label_vector == 1);   % FEP (unchanged by permutation)
n0 = sum(label_vector == 0);   % HC
max_nc = max(ncomp_list);

%% === Main loop: edge sets ===
t_start = tic;
for inet = 1:num_edge_sets
    if done(inet), continue; end
    t_set = tic;

    % --- Raw features ---
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

    % --- Permutation 0 (true labels) + n_perm nulls, all ncomp at once ---
    acc_tmp = nan(n_perm + 1, n_nc);
    dp_tmp  = nan(n_perm + 1, n_nc);
    parfor p = 1:(n_perm + 1)          % runs as a normal for-loop without the Parallel Toolbox
        y = Y_all(:, p);
        pred = nan(num_subs, n_nc);
        for k = 1:K
            tr = train_idx_all{k};  te = test_idx_all{k};
            ncomp_fit = min([sum(tr) - 1, size(x_train_all{k}, 2), max_nc]);
            [~, YL, ~, ~, ~, ~, ~, stats] = plsregress(x_train_all{k}, y(tr), ncomp_fit);
            for c = 1:n_nc
                nc   = min(ncomp_list(c), ncomp_fit);   % capped for edge sets with < 3 edges
                B    = stats.W(:, 1:nc) * YL(:, 1:nc)';
                BETA = [mean(y(tr)) - mean(x_train_all{k}, 1) * B; B];
                pred(te, c) = double([ones(sum(te), 1), x_test_all{k}] * BETA >= 0.5);
            end
        end
        acc_row = nan(1, n_nc);
        dp_row  = nan(1, n_nc);
        for c = 1:n_nc
            acc_row(c) = mean(pred(:, c) == y);
            dp_row(c)  = z_rate(sum(pred(:, c) == 1 & y == 1), n1) - z_rate(sum(pred(:, c) == 1 & y == 0), n0);
        end
        acc_tmp(p, :) = acc_row;
        dp_tmp(p, :)  = dp_row;
    end

    obs_acc(inet, :)        = acc_tmp(1, :);
    obs_dprime(inet, :)     = dp_tmp(1, :);
    null_acc(inet, :, :)    = reshape(acc_tmp(2:end, :), 1, n_perm, n_nc);
    null_dprime(inet, :, :) = reshape(dp_tmp(2:end, :),  1, n_perm, n_nc);
    done(inet) = true;

    % --- Checkpoint + progress ---
    save(checkpoint_path, 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', 'done', ...
        'ckpt_subjects', 'ckpt_ncomp', 'ckpt_nperm', '-v7.3');
    n_done_now = sum(done);
    elapsed = toc(t_start);
    fprintf('[%03d/%03d] %-35s d'' 1/2/3comp = %5.2f %5.2f %5.2f | %.0fs | ~%.1f h left\n', ...
        inet, num_edge_sets, model_labels{inet}, obs_dprime(inet, :), toc(t_set), ...
        elapsed / max(n_done_now, 1) * (num_edge_sets - n_done_now) / 3600);
    drawnow;
end

%% === p-values + BH-FDR (across edge sets, separately per ncomp) ===
p_acc    = nan(num_edge_sets, n_nc);  q_acc    = nan(num_edge_sets, n_nc);
p_dprime = nan(num_edge_sets, n_nc);  q_dprime = nan(num_edge_sets, n_nc);
for c = 1:n_nc
    p_acc(:, c)    = (sum(null_acc(:, :, c)    >= obs_acc(:, c),    2) + 1) / (n_perm + 1);
    p_dprime(:, c) = (sum(null_dprime(:, :, c) >= obs_dprime(:, c), 2) + 1) / (n_perm + 1);
    [ps, idx] = sort(p_acc(:, c));     q_acc(idx, c)    = min(1, cummin(ps .* num_edge_sets ./ (1:num_edge_sets)', 'reverse'));
    [ps, idx] = sort(p_dprime(:, c));  q_dprime(idx, c) = min(1, cummin(ps .* num_edge_sets ./ (1:num_edge_sets)', 'reverse'));
end

%% === One summary CSV per component count + long table ===
perm_results = table();
for c = 1:n_nc
    T = table(repmat(ncomp_list(c), num_edge_sets, 1), (1:num_edge_sets)', model_labels(:), ...
        100*obs_acc(:, c), p_acc(:, c), q_acc(:, c), ...
        obs_dprime(:, c), p_dprime(:, c), q_dprime(:, c), ...
        100*mean(null_acc(:, :, c), 2), mean(null_dprime(:, :, c), 2), ...
        'VariableNames', {'NComp', 'EdgeSet', 'Model', ...
        'Acc_FixedPartition', 'Acc_perm_p', 'Acc_perm_q', ...
        'Dprime_FixedPartition', 'Dprime_perm_p', 'Dprime_perm_q', ...
        'NullAcc_Mean', 'NullDprime_Mean'});
    writetable(T, fullfile(perm_output_dir, sprintf('FEP_PermutationTest_FEPvsHC_Summary_%dcomp.csv', ncomp_list(c))));
    perm_results = [perm_results; T];
end
writetable(perm_results, fullfile(perm_output_dir, ['FEP_PermutationTest_FEPvsHC_ByNComp_' ncomp_str 'comp.csv']));

%% === Wide table: one row per model, ncomp side by side ===
perm_cmp = table((1:num_edge_sets)', model_labels(:), 'VariableNames', {'EdgeSet', 'Model'});
for c = 1:n_nc
    t = sprintf('_%dcomp', ncomp_list(c));
    perm_cmp.(['Dprime' t])        = obs_dprime(:, c);
    perm_cmp.(['Dprime_perm_p' t]) = p_dprime(:, c);
    perm_cmp.(['Dprime_perm_q' t]) = q_dprime(:, c);
end
writetable(perm_cmp, fullfile(perm_output_dir, ['FEP_PermutationTest_FEPvsHC_NCompComparison_' ncomp_str 'comp.csv']));

save(fullfile(perm_output_dir, ['FEP_PermutationTest_FEPvsHC_' ncomp_str 'comp.mat']), ...
    'perm_results', 'perm_cmp', 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', ...
    'p_acc', 'p_dprime', 'q_acc', 'q_dprime', 'ncomp_list', 'perm_idx', 'n_perm', 'K', ...
    'cv_seed', 'perm_seed', 'train_idx_all', 'test_idx_all', 'model_labels', ...
    'subject_id_vector', '-v7.3');

fprintf('\n=== Permutation test by number of PLS components ===\n');
for c = 1:n_nc
    fprintf('%d comp: %3d of %d edge sets significant (d'', FDR q < .05) | mean d'' = %.2f\n', ...
        ncomp_list(c), sum(q_dprime(:, c) < .05), num_edge_sets, mean(obs_dprime(:, c)));
end
fprintf('Saved to: %s\n', perm_output_dir);