%% LSD-Placebo Classifer 
% 1) loads data
% 2) PLS-DA classifer
% 3) Permuatation tesing 

%% Loadin Dating 
% Three Basel studies 
%   5-HT, LAM, LPM
% Define directories for LSD datasets (including placebo scans)
data_root = '/Volumes/toshiba/LSD_November11/LSD_BASELS/LSD_Basel/combinedsubcorticalOnly_NO_GSR';
HT5_dir = fullfile(data_root, '5-HT');
LAM_dir = fullfile(data_root, 'LAM');
LPM_dir = fullfile(data_root, 'LPM');

% Covariates
covariates_dir = '/Volumes/toshiba/LSD_November11/LSD_BASELS/Covariates';

% Output directory
output_dir = '/Volumes/toshiba/LSD_November11/LSD_BASELS/OutputBaselOnly/schaeferSubcortical_Nov17';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Setup Paths
% Placebo covariate files
HT5_Placebo_cov_file = fullfile(covariates_dir, 'Basel_5-HT_placebo_only_passed.csv');
LAM_Placebo_cov_file = fullfile(covariates_dir, 'Basel_LAM_placebo_only_passed.csv');
LPM_Placebo_cov_file = fullfile(covariates_dir, 'LPM_covariates_Placebo_only.csv');

% LSD covariate files
HT5_LSD_cov_file = fullfile(covariates_dir, 'Basel_5-HT_covariates_LSD_only.csv');
LAM_LSD_cov_file = fullfile(covariates_dir, 'Basel_LAM_covariates_LSD_only.csv');
LPM_LSD_cov_file = fullfile(covariates_dir, 'LPM_covariates_LSD_only.csv');


%% =========================== LSD LOADING =========================== %%
% 5-HT LSD 
HT5_LSD_cov_table = readtable(HT5_LSD_cov_file);
HT5_Placebo_cov_table_check = readtable(HT5_Placebo_cov_file); % to identify placebo-only subs

ht5_lsd_files = dir(fullfile(HT5_dir, 'sub-*_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_ht5_lsd = regexp({ht5_lsd_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_ht5_lsd = unique(cellfun(@(x) x{1}, all_ids_ht5_lsd, 'UniformOutput', false));

% Identify subjects present in Placebo covariates but absent from LSD file list
placebo_ids_all = HT5_Placebo_cov_table_check.Subject;
placebo_ids_all = placebo_ids_all(~cellfun(@isempty, placebo_ids_all));
placebo_ids_numeric = regexp(placebo_ids_all, 'sub-(\d+)', 'tokens', 'once');
placebo_ids_numeric = cellfun(@(x) x{1}, placebo_ids_numeric, 'UniformOutput', false);
ht5_placebo_only_ids = setdiff(placebo_ids_numeric, all_ids_ht5_lsd);
if ~isempty(ht5_placebo_only_ids)
    fprintf('Excluding %d placebo-only 5-HT subject(s) not present in LSD: %s\n', ...
        length(ht5_placebo_only_ids), strjoin(ht5_placebo_only_ids, ', '));
end

LSD_UsedHT5 = sort(setdiff(all_ids_ht5_lsd, ht5_placebo_only_ids));
num_HT5_LSD = length(LSD_UsedHT5);

fc_matrices_HT5_LSD = cell(num_HT5_LSD, 1);
subject_ids_HT5_LSD = strings(num_HT5_LSD, 1);
mean_fd_HT5_LSD = nan(num_HT5_LSD, 1);
age_HT5_LSD = nan(num_HT5_LSD, 1);
male_HT5_LSD = nan(num_HT5_LSD, 1);

for i = 1:num_HT5_LSD
    sub_id = LSD_UsedHT5{i};
    pattern = sprintf('sub-%s_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(HT5_dir, pattern));
    if isempty(file)
        warning('File not found for 5-HT subject %s', sub_id);
        continue;
    end
    file_path = fullfile(HT5_dir, file(1).name);
    fc_matrices_HT5_LSD{i} = readmatrix(file_path);
    subject_ids_HT5_LSD(i) = string(['5-HT_sub-' sub_id]);

    cov_idx = find(strcmp(HT5_LSD_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    age_HT5_LSD(i) = HT5_LSD_cov_table.Age(cov_idx);
    male_HT5_LSD(i) = HT5_LSD_cov_table.Male(cov_idx);
    mean_fd_HT5_LSD(i) = HT5_LSD_cov_table.Mean_framewise_displacement(cov_idx);
end

save(fullfile(output_dir, '5HT_LSD_Selections.mat'), ...
    'LSD_UsedHT5', 'fc_matrices_HT5_LSD', 'subject_ids_HT5_LSD', ...
    'mean_fd_HT5_LSD', 'age_HT5_LSD', 'male_HT5_LSD');

fprintf('Loaded %d 5-HT LSD matrices and matched covariates.\n', num_HT5_LSD);


% LAM LSD 
LAM_LSD_cov_table = readtable(LAM_LSD_cov_file);
LAM_Placebo_cov_table_check = readtable(LAM_Placebo_cov_file); % to identify placebo-only subs

lam_lsd_files = dir(fullfile(LAM_dir, 'sub-*_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_lam_lsd = regexp({lam_lsd_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_lam_lsd = unique(cellfun(@(x) x{1}, all_ids_lam_lsd, 'UniformOutput', false));

% Identify subjects present in Placebo covariates but absent from LSD file list
placebo_ids_all = LAM_Placebo_cov_table_check.Subject;
placebo_ids_all = placebo_ids_all(~cellfun(@isempty, placebo_ids_all));
placebo_ids_numeric = regexp(placebo_ids_all, 'sub-(\d+)', 'tokens', 'once');
placebo_ids_numeric = cellfun(@(x) x{1}, placebo_ids_numeric, 'UniformOutput', false);
lam_placebo_only_ids = setdiff(placebo_ids_numeric, all_ids_lam_lsd);
if ~isempty(lam_placebo_only_ids)
    fprintf('Excluding %d placebo-only LAM subject(s) not present in LSD: %s\n', ...
        length(lam_placebo_only_ids), strjoin(lam_placebo_only_ids, ', '));
end

LSD_UsedLAM = sort(setdiff(all_ids_lam_lsd, lam_placebo_only_ids));
num_LAM_LSD = length(LSD_UsedLAM);

fc_matrices_LAM_LSD = cell(num_LAM_LSD, 1);
subject_ids_LAM_LSD = strings(num_LAM_LSD, 1);
mean_fd_LAM_LSD = nan(num_LAM_LSD, 1);
age_LAM_LSD = nan(num_LAM_LSD, 1);
male_LAM_LSD = nan(num_LAM_LSD, 1);

for i = 1:num_LAM_LSD
    sub_id = LSD_UsedLAM{i};
    pattern = sprintf('sub-%s_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(LAM_dir, pattern));
    if isempty(file)
        warning('File not found for LAM subject %s', sub_id);
        continue;
    end
    file_path = fullfile(LAM_dir, file(1).name);
    fc_matrices_LAM_LSD{i} = readmatrix(file_path);
    subject_ids_LAM_LSD(i) = string(['LAM_sub-' sub_id]);

    cov_idx = find(strcmp(LAM_LSD_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    age_LAM_LSD(i) = LAM_LSD_cov_table.Age(cov_idx);
    male_LAM_LSD(i) = LAM_LSD_cov_table.Male(cov_idx);
    mean_fd_LAM_LSD(i) = LAM_LSD_cov_table.Mean_framewise_displacement(cov_idx);
end

save(fullfile(output_dir, 'LAM_LSD_Selections.mat'), ...
    'LSD_UsedLAM', 'fc_matrices_LAM_LSD', 'subject_ids_LAM_LSD', ...
    'mean_fd_LAM_LSD', 'age_LAM_LSD', 'male_LAM_LSD');

fprintf('Loaded %d LAM LSD matrices and matched covariates.\n', num_LAM_LSD);


% LPM LSD Select 
LPM_LSD_cov_table = readtable(LPM_LSD_cov_file);

lpm_lsd_files = dir(fullfile(LPM_dir, 'sub-*_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_lpm_lsd = regexp({lpm_lsd_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_lpm_lsd = unique(cellfun(@(x) x{1}, all_ids_lpm_lsd, 'UniformOutput', false));

LSD_UsedLPM = sort(all_ids_lpm_lsd);
num_LPM_LSD = length(LSD_UsedLPM);

fc_matrices_LPM_LSD = cell(num_LPM_LSD, 1);
subject_ids_LPM_LSD = strings(num_LPM_LSD, 1);
mean_fd_LPM_LSD = nan(num_LPM_LSD, 1);
age_LPM_LSD = nan(num_LPM_LSD, 1);
male_LPM_LSD = nan(num_LPM_LSD, 1);

for i = 1:num_LPM_LSD
    sub_id = LSD_UsedLPM{i};
    pattern = sprintf('sub-%s_ses-LSD*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(LPM_dir, pattern));
    if isempty(file)
        warning('File not found for LPM subject %s', sub_id);
        continue;
    end
    file_path = fullfile(LPM_dir, file(1).name);
    fc_matrices_LPM_LSD{i} = readmatrix(file_path);
    subject_ids_LPM_LSD(i) = string(['LPM_sub-' sub_id]);

    cov_idx = find(strcmp(LPM_LSD_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    age_LPM_LSD(i) = LPM_LSD_cov_table.Age(cov_idx);
    male_LPM_LSD(i) = LPM_LSD_cov_table.Male(cov_idx);
    mean_fd_LPM_LSD(i) = LPM_LSD_cov_table.Mean_framewise_displacement(cov_idx);
end

save(fullfile(output_dir, 'LPM_LSD_Selections.mat'), ...
    'LSD_UsedLPM', 'fc_matrices_LPM_LSD', 'subject_ids_LPM_LSD', ...
    'mean_fd_LPM_LSD', 'age_LPM_LSD', 'male_LPM_LSD');

fprintf('Loaded %d LPM LSD matrices and matched covariates.\n', num_LPM_LSD);



%% =========================== PLACEBO LOADING =========================== %%
% Every LSD subject needs a matched placebo scan, and vice versa within study.

% 5-HT Placebo Select 
HT5_Placebo_cov_table = readtable(HT5_Placebo_cov_file);

ht5_placebo_files = dir(fullfile(HT5_dir, 'sub-*_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_ht5_plcb = regexp({ht5_placebo_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_ht5_plcb = unique(cellfun(@(x) x{1}, all_ids_ht5_plcb, 'UniformOutput', false));

Placebo_UsedHT5 = sort(all_ids_ht5_plcb);
num_HT5_Placebo = length(Placebo_UsedHT5);

fc_matrices_HT5_Placebo = cell(num_HT5_Placebo, 1);
subject_ids_HT5_Placebo = strings(num_HT5_Placebo, 1);
mean_fd_HT5_Placebo = nan(num_HT5_Placebo, 1);
age_HT5_Placebo = nan(num_HT5_Placebo, 1);
male_HT5_Placebo = nan(num_HT5_Placebo, 1);

for i = 1:num_HT5_Placebo
    sub_id = Placebo_UsedHT5{i};
    pattern = sprintf('sub-%s_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(HT5_dir, pattern));
    if isempty(file)
        warning('File not found for 5-HT subject %s', sub_id);
        continue;
    end
    file_path = fullfile(HT5_dir, file(1).name);
    fc_matrices_HT5_Placebo{i} = readmatrix(file_path);
    subject_ids_HT5_Placebo(i) = string(['5-HT_sub-' sub_id]);

    cov_idx = find(strcmp(HT5_Placebo_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    mean_fd_HT5_Placebo(i) = HT5_Placebo_cov_table.Mean_framewise_displacement(cov_idx);
    age_HT5_Placebo(i) = HT5_Placebo_cov_table.Age(cov_idx);
    male_HT5_Placebo(i) = HT5_Placebo_cov_table.Male(cov_idx);
end

fprintf('5-HT Placebo: Missing FD values: %d out of %d\n', sum(isnan(mean_fd_HT5_Placebo)), num_HT5_Placebo);

save(fullfile(output_dir, '5HT_Placebo_Selections.mat'), ...
    'Placebo_UsedHT5', 'fc_matrices_HT5_Placebo', 'subject_ids_HT5_Placebo', ...
    'mean_fd_HT5_Placebo', 'age_HT5_Placebo', 'male_HT5_Placebo');

fprintf('Loaded %d 5-HT placebo matrices and matched covariates.\n', num_HT5_Placebo);


% LAM Placebo Select 
LAM_Placebo_cov_table = readtable(LAM_Placebo_cov_file);

lam_placebo_files = dir(fullfile(LAM_dir, 'sub-*_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_lam_plcb = regexp({lam_placebo_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_lam_plcb = unique(cellfun(@(x) x{1}, all_ids_lam_plcb, 'UniformOutput', false));

Placebo_UsedLAM = sort(all_ids_lam_plcb);
num_LAM_Placebo = length(Placebo_UsedLAM);

fc_matrices_LAM_Placebo = cell(num_LAM_Placebo, 1);
subject_ids_LAM_Placebo = strings(num_LAM_Placebo, 1);
mean_fd_LAM_Placebo = nan(num_LAM_Placebo, 1);
age_LAM_Placebo = nan(num_LAM_Placebo, 1);
male_LAM_Placebo = nan(num_LAM_Placebo, 1);

for i = 1:num_LAM_Placebo
    sub_id = Placebo_UsedLAM{i};
    pattern = sprintf('sub-%s_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(LAM_dir, pattern));
    if isempty(file)
        warning('File not found for LAM subject %s', sub_id);
        continue;
    end
    file_path = fullfile(LAM_dir, file(1).name);
    fc_matrices_LAM_Placebo{i} = readmatrix(file_path);
    subject_ids_LAM_Placebo(i) = string(['LAM_sub-' sub_id]);

    cov_idx = find(strcmp(LAM_Placebo_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    mean_fd_LAM_Placebo(i) = LAM_Placebo_cov_table.Mean_framewise_displacement(cov_idx);
    age_LAM_Placebo(i) = LAM_Placebo_cov_table.Age(cov_idx);
    male_LAM_Placebo(i) = LAM_Placebo_cov_table.Male(cov_idx);
end

fprintf('LAM Placebo: Missing FD values: %d out of %d\n', sum(isnan(mean_fd_LAM_Placebo)), num_LAM_Placebo);

save(fullfile(output_dir, 'LAM_Placebo_Selections.mat'), ...
    'Placebo_UsedLAM', 'fc_matrices_LAM_Placebo', 'subject_ids_LAM_Placebo', ...
    'mean_fd_LAM_Placebo', 'age_LAM_Placebo', 'male_LAM_Placebo');

fprintf('Loaded %d LAM placebo matrices and matched covariates.\n', num_LAM_Placebo);


% LPM Placebo Select 
LPM_Placebo_cov_table = readtable(LPM_Placebo_cov_file);

lpm_placebo_files = dir(fullfile(LPM_dir, 'sub-*_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv'));
all_ids_lpm_plcb = regexp({lpm_placebo_files.name}, 'sub-(\d+)_', 'tokens', 'once');
all_ids_lpm_plcb = unique(cellfun(@(x) x{1}, all_ids_lpm_plcb, 'UniformOutput', false));

Placebo_UsedLPM = sort(all_ids_lpm_plcb);
num_LPM_Placebo = length(Placebo_UsedLPM);

fc_matrices_LPM_Placebo = cell(num_LPM_Placebo, 1);
subject_ids_LPM_Placebo = strings(num_LPM_Placebo, 1);
mean_fd_LPM_Placebo = nan(num_LPM_Placebo, 1);
age_LPM_Placebo = nan(num_LPM_Placebo, 1);
male_LPM_Placebo = nan(num_LPM_Placebo, 1);

for i = 1:num_LPM_Placebo
    sub_id = Placebo_UsedLPM{i};
    pattern = sprintf('sub-%s_ses-plcb*combinedsubcorticalOnly_TS_denoised_NO_GSR*.csv', sub_id);
    file = dir(fullfile(LPM_dir, pattern));
    if isempty(file)
        warning('File not found for LPM subject %s', sub_id);
        continue;
    end
    file_path = fullfile(LPM_dir, file(1).name);
    fc_matrices_LPM_Placebo{i} = readmatrix(file_path);
    subject_ids_LPM_Placebo(i) = string(['LPM_sub-' sub_id]);

    cov_idx = find(strcmp(LPM_Placebo_cov_table.Subject, ['sub-' sub_id]));
    if isempty(cov_idx)
        warning('Covariates not found for %s', sub_id);
        continue;
    end
    mean_fd_LPM_Placebo(i) = LPM_Placebo_cov_table.Mean_framewise_displacement(cov_idx);
    age_LPM_Placebo(i) = LPM_Placebo_cov_table.Age(cov_idx);
    male_LPM_Placebo(i) = LPM_Placebo_cov_table.Male(cov_idx);
end

fprintf('LPM Placebo: Missing FD values: %d out of %d\n', sum(isnan(mean_fd_LPM_Placebo)), num_LPM_Placebo);

save(fullfile(output_dir, 'LPM_Placebo_Selections.mat'), ...
    'Placebo_UsedLPM', 'fc_matrices_LPM_Placebo', 'subject_ids_LPM_Placebo', ...
    'mean_fd_LPM_Placebo', 'age_LPM_Placebo', 'male_LPM_Placebo');

fprintf('Loaded %d LPM placebo matrices and matched covariates.\n', num_LPM_Placebo);


%% === Input for LSD vs Placebo Classification  ===
% Combine FC matrices (LSD first, then Placebo; 5-HT, LAM, LPM — matches loading order above)
fc_matrices_all = [fc_matrices_HT5_LSD;     fc_matrices_LAM_LSD;     fc_matrices_LPM_LSD; ...
                   fc_matrices_HT5_Placebo; fc_matrices_LAM_Placebo; fc_matrices_LPM_Placebo];
% Combine subject IDs 
subject_ids = [subject_ids_HT5_LSD;     subject_ids_LAM_LSD;     subject_ids_LPM_LSD; ...
               subject_ids_HT5_Placebo; subject_ids_LAM_Placebo; subject_ids_LPM_Placebo];

% Combine ALL covariates: Age, Male, Mean FD
age_all = [age_HT5_LSD;     age_LAM_LSD;     age_LPM_LSD; ...
           age_HT5_Placebo; age_LAM_Placebo; age_LPM_Placebo];
male_all = [male_HT5_LSD;     male_LAM_LSD;     male_LPM_LSD; ...
            male_HT5_Placebo; male_LAM_Placebo; male_LPM_Placebo];
mean_fd_all = [mean_fd_HT5_LSD;     mean_fd_LAM_LSD;     mean_fd_LPM_LSD; ...
               mean_fd_HT5_Placebo; mean_fd_LAM_Placebo; mean_fd_LPM_Placebo];

% Subject ID vector for within-subject pairing in CV folds
subject_id_vector = subject_ids;

% === Define Binary Class Labels for LSD and Placebo ===
num_LSD     = num_HT5_LSD + num_LAM_LSD + num_LPM_LSD;
num_Placebo = num_HT5_Placebo + num_LAM_Placebo + num_LPM_Placebo;

num_subs = num_LSD + num_Placebo;

% Sanity check: label count must match loaded FC matrices
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

% === Additional sanity checks and final save ===
if num_subs ~= length(subject_ids) || num_subs ~= length(age_all) ...
        || num_subs ~= length(male_all) ...
        || num_subs ~= length(mean_fd_all) || num_subs ~= length(label_vector)
    error('Mismatch in number of subjects across matrices, covariates, labels, or IDs');
end

if any(isnan(age_all)) || any(isnan(male_all)) || any(isnan(mean_fd_all))
    error('NaNs present in covariates — fix before residualizing (B_cov would be all NaN)');
end

% Build final FC array: ROI x ROI x subjects
sfc_mat = cat(3, fc_matrices_all{:});

% Build covariates matrix: columns = [Age, Male, Mean_FD]
covariates = [age_all, male_all, mean_fd_all];

% Keep only subjects with BOTH an LSD and a Placebo scan (fully paired)
[unique_ids, ~, ic] = unique(subject_id_vector);
keep_mask = ismember(subject_id_vector, unique_ids(accumarray(ic, 1) == 2));
fprintf('Dropping %d unpaired scan(s): %s\n', sum(~keep_mask), strjoin(subject_id_vector(~keep_mask), ', '));
sfc_mat           = sfc_mat(:,:,keep_mask);
covariates        = covariates(keep_mask, :);
label_vector      = label_vector(keep_mask);
subject_id_vector = subject_id_vector(keep_mask);
subject_ids       = subject_ids(keep_mask);     % filtered too, so it never goes stale
num_subs          = size(sfc_mat, 3);

% Save outputs
save(fullfile(output_dir, 'BASEL_DATASETS_sfc_mat_LSD_vs_Placebo.mat'), 'sfc_mat', '-v7.3');
save(fullfile(output_dir, 'BASEL_DATASETS_covariates_LSD_vs_Placebo.mat'), 'covariates');
save(fullfile(output_dir, 'BASEL_DATASETS_label_vector_LSD_vs_Placebo.mat'), 'label_vector');
save(fullfile(output_dir, 'BASEL_DATASETS_subject_id_vector_LSD_vs_Placebo.mat'), 'subject_id_vector');
writematrix(subject_ids, fullfile(output_dir, 'BASEL_DATASETS_subject_ids_LSD_vs_Placebo.txt'));

fprintf('Saved sFC matrix, full covariates (Age/Male/FD), labels, and subject IDs for LSD vs Placebo classification.\n');

%% === Define Paths ===
basel_output_dir = output_dir;   % same folder as the loading block
edge_dir       = '/Volumes/toshiba/LSD_November11/Within_Between_allEdge/Schaefer100_Tian38';
edge_file      = fullfile(edge_dir, 'edge_sets_and_labels.mat');
num_rois    = size(sfc_mat, 1);
num_subs    = size(sfc_mat, 3);
num_LSD     = sum(label_vector == 1);   % recomputed AFTER the pairing filter
num_Placebo = sum(label_vector == 0);

%% === Load Edge Sets + Model Labels (built and saved TOGETHER) ===
if ~exist(edge_file, 'file'), error('Missing file: %s', edge_file);
end
E = load(edge_file, 'edge_sets', 'model_labels', 'model_type', ...
                    'model_i', 'model_j', 'group_labels', 'roi_group');

all_net_edges = E.edge_sets;          % 326 x 1: full, 25 within, 300 between
model_labels  = E.model_labels;       % 326 x 1, same order as all_net_edges
model_type    = E.model_type;         % 'full' / 'within' / 'between'
group_labels  = E.group_labels;       % the 25 group names
num_nets      = numel(group_labels);

%% === Sanity Checks ===
% Scans: matrices, labels, covariates and IDs all line up
assert(num_subs == numel(label_vector) && num_subs == size(covariates, 1) ...
    && num_subs == numel(subject_id_vector), ...
    'Scan count mismatch: sfc_mat %d, labels %d, covariates %d, IDs %d', ...
    num_subs, numel(label_vector), size(covariates, 1), numel(subject_id_vector));
assert(all(ismember(label_vector, [0 1])), 'label_vector must contain only 0 (Placebo) and 1 (LSD)');
assert(isequal(label_vector(:), [ones(num_LSD,1); zeros(num_Placebo,1)]), ...
    'label_vector is not ordered LSD first, then Placebo (expected by the loading block)');
assert(size(covariates, 2) == 3, 'covariates should be [Age, Male, Mean_FD]');
assert(~any(isnan(covariates(:))), 'NaNs present in covariates');

% Within-subject pairing: every subject has exactly one LSD and one Placebo scan
[unique_ids, ~, ic] = unique(subject_id_vector);
assert(all(accumarray(ic, 1) == 2), 'Some subjects do not have exactly 2 scans');
assert(all(accumarray(ic, label_vector) == 1), 'Some subjects do not have one LSD + one Placebo scan');

% Every subject ID carries its study label (5-HT_ / LAM_ / LPM_)
assert(all(startsWith(subject_id_vector, ["5-HT_", "LAM_", "LPM_"])), ...
    'Subject IDs must start with 5-HT_, LAM_ or LPM_');

fprintf('Loaded sfc_mat [%d x %d x %d]: LSD = %d | Placebo = %d | %d subjects\n', ...
    num_rois, num_rois, num_subs, num_LSD, num_Placebo, numel(unique_ids));
fprintf('Loaded %d edge sets (%d full, %d within, %d between) for %d ROIs\n', ...
    numel(all_net_edges), sum(strcmp(model_type,'full')), sum(strcmp(model_type,'within')), ...
    sum(strcmp(model_type,'between')), num_rois);
fprintf('%d model labels match %d edge sets (verified edge-by-edge).\n', ...
    numel(model_labels), numel(all_net_edges));

%% === Run Cross-Validation Binary PLS Classifier (LSD vs Placebo, Within-Subject) ===
% ===============================================================
% LSD vs Placebo: Binary PLS Classifier with FOLD-WISE Covariate Control
% WITHIN-SUBJECT design — each subject contributes an LSD scan AND a
% Placebo scan. Both scans of a subject are ALWAYS kept in the same fold.
% Run once each for 1, 2 and 3 PLS components (ncomp_list).
% ===============================================================

% Configuration
num_repeats = 100;
K = 10;
ncomp_list = [1 2 3];    % comparing results by number of componets 
verbose_folds = false;   % true = print every fold (326 sets x 100 repeats x 10 folds = 326,000 lines per ncomp)

output_dir   = basel_output_dir;          
model_labels = cellstr(model_labels);     

% Number of edge/network models
num_edge_sets = length(all_net_edges);
assert(numel(model_labels) == num_edge_sets, ...
    'model_labels (%d) and all_net_edges (%d) differ in length', numel(model_labels), num_edge_sets);
num_subs = length(label_vector); % total rows (subjects x 2 sessions)
assert(num_subs == size(sfc_mat, 3) && num_subs == size(covariates, 1) && num_subs == numel(subject_id_vector), ...
    'label_vector (%d), sfc_mat (%d), covariates (%d) and subject_id_vector (%d) differ', ...
    num_subs, size(sfc_mat, 3), size(covariates, 1), numel(subject_id_vector));

% Unique subjects, fixed order (study-prefixed IDs, e.g. 5-HT_sub-001)
unique_subjects = unique(subject_id_vector, 'stable');
num_unique_subs = length(unique_subjects);
fprintf('%d scans across %d subjects (within-subject): LSD=%d | Placebo=%d\n', ...
    num_subs, num_unique_subs, sum(label_vector == 1), sum(label_vector == 0));

% -------------------------------------------------------
% STUDY for each subject, from the ID prefix (5-HT_ / LAM_ / LPM_)
% -------------------------------------------------------
site_labels_unique = strings(num_unique_subs, 1);
site_labels_unique(startsWith(unique_subjects, "5-HT_")) = "5-HT";
site_labels_unique(startsWith(unique_subjects, "LAM_"))  = "LAM";
site_labels_unique(startsWith(unique_subjects, "LPM_"))  = "LPM";
if any(site_labels_unique == "")
    error('Unrecognized study prefix for subject(s): %s', strjoin(unique_subjects(site_labels_unique == ""), ', '));
end
fprintf('Study breakdown: 5-HT=%d | LAM=%d | LPM=%d (total %d)\n', ...
    sum(site_labels_unique == "5-HT"), sum(site_labels_unique == "LAM"), ...
    sum(site_labels_unique == "LPM"), num_unique_subs);

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

    % held-out predictions across repeats (both LSD and Placebo rows)
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
        % PAIRED, STUDY-STRATIFIED SUBJECT-LEVEL K-FOLD PARTITION
        %
        % Built ONCE per repeat and shared by all edge sets. Partition is on
        % unique SUBJECTS, stratified by STUDY (5-HT / LAM / LPM); both scans
        % of a subject are then expanded into the same fold below.
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
            % loop over CV folds
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

                if verbose_folds, fprintf(' [%dcomp] Edge %d (%s) | Fold %d — TrainAcc: %.3f | TestAcc: %.3f | Placebo: %.3f, LSD: %.3f\n', ...
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

            % Sanity check: every scan should have one held-out prediction
            if any(isnan(heldout_class_this_repeat)) || any(isnan(heldout_score_this_repeat))
                warning('Missing held-out predictions in Repeat %d, Edge Set %d (%dcomp)', r, inet, max_ncomp);
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
    % with the FEP/HC pipeline, and means the saved model was never fit on
    % 100% of the data.
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

    model_path = fullfile(output_dir, ['AveragedPLSModels_LSDvsPlacebo' ncomp_tag '.mat']);
    save(model_path, 'avg_BETA_final', 'avg_B_cov_all', 'avg_x_mean_all', 'avg_x_std_all', ...
        'model_labels', 'max_ncomp');
    fprintf('Saved averaged models to: %s\n', model_path);


    % === Save all subject-level held-out predictions (LSD vs Placebo) ===
    heldout_save_path = fullfile(output_dir, ['BASEL_HeldOutPredictions_AllEdgeSets_LSDvsPlacebo' ncomp_tag '.mat']);
    save(heldout_save_path, 'heldout_class_all', 'heldout_score_all', ...
        'label_vector', 'subject_id_vector', 'model_labels', 'max_ncomp', '-v7.3');
    fprintf('Saved subject-level held-out predictions to: %s\n', heldout_save_path);


    % === Subject-level summary CSV: average across repeats ===
    subject_id_cell = cellstr(subject_id_vector);

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

    summary_csv_path = fullfile(output_dir, ['BASEL_HeldOutPredictions_SubjectSummary_LSDvsPlacebo' ncomp_tag '.csv']);
    writetable(T_summary, summary_csv_path);
    fprintf('Saved subject-summary CSV: %s\n', summary_csv_path);


    % =========================================================================
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
            misses = sum(pred == 0 & label_vector == 1);   % LSD called placebo
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
        per_iteration_summary{:, {'TestAcc','Placebo_Acc','LSD_Acc'}} * 100;  

    per_iteration_path = fullfile(output_dir, ['BASEL_FinalSummary_LSDvsPlacebo_PerIteration_WithDprime' ncomp_tag '.csv']);
    writetable(per_iteration_summary, per_iteration_path);
    fprintf('Per-iteration summary (%d rows) saved to: %s\n', n_iter_rows, per_iteration_path);


    % =========================================================================
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

    per_model_path = fullfile(output_dir, ['BASEL_FinalSummary_LSDvsPlacebo_PerModel_MeanAcrossRepeats' ncomp_tag '.csv']);
    writetable(per_model_summary, per_model_path);
    fprintf('Per-model summary (mean/SD/CI across %d repeats) saved to: %s\n', ...
        n_repeats_observed, per_model_path);

    % --- Keep this component count's results for the comparison table ---
    dprime_by_ncomp(:, ic) = per_model_summary.Dprime_Mean;
    dprime_sd_ncomp(:, ic) = per_model_summary.Dprime_SD;
    acc_by_ncomp(:, ic)    = per_model_summary.TestAcc_Mean;

    fprintf('===== Finished %d component(s) in %.1f h =====\n', max_ncomp, toc(t_ncomp)/3600);
end   % component loop


%% === Compare component counts  ===
comp_table = table((1:num_edge_sets)', model_labels(:), 'VariableNames', {'EdgeSet', 'Model'});
for ic = 1:numel(ncomp_list)
    t = sprintf('_%dcomp', ncomp_list(ic));
    comp_table.(['Dprime_Mean' t])  = dprime_by_ncomp(:, ic);
    comp_table.(['Dprime_SD' t])    = dprime_sd_ncomp(:, ic);
    comp_table.(['TestAcc_Mean' t]) = acc_by_ncomp(:, ic);
    fprintf('%d component(s): full-connectome d'' = %.2f | mean d'' across edge sets = %.2f\n', ...
        ncomp_list(ic), dprime_by_ncomp(1, ic), mean(dprime_by_ncomp(:, ic)));
end
comp_path = fullfile(output_dir, sprintf('BASEL_CompareNcomp_LSDvsPlacebo_%scomp.csv', ...
    strjoin(arrayfun(@num2str, ncomp_list, 'UniformOutput', false), '-')));
writetable(comp_table, comp_path);
fprintf('Saved component comparison to: %s\n', comp_path);

%% =========================================================================
% PERMUTATION TEST — LSD vs Placebo (5-HT / LAM / LPM), 1, 2 AND 3 PLS COMPONENTS
% -------------------------------------------------------------------------
% Run AFTER the paths / edge-set block (needs sfc_mat, covariates,
% label_vector, subject_id_vector, all_net_edges, model_labels,
% basel_output_dir).
%
% Design:
%   - ONE fixed 10-fold partition (cv_seed), subject-level, stratified by
%     STUDY (5-HT / LAM / LPM). Both scans of a subject stay in one fold.
%   - Fold-wise residualization + z-scoring do not use labels, so they are
%     computed ONCE per edge set and reused for every permutation.
%   - Permutation 0 = true labels; 1..n_perm = within-subject swaps
%     (each subject's LSD/Placebo labels swapped with p = .5), the SAME
%     swaps for every edge set AND every component count.
%   - Each PLS fit uses 3 components; the 1- and 2-component models are
%     rebuilt from its first components (PLS components are sequential, so
%     this is identical to fitting 1 or 2 components directly). All three
%     component counts therefore cost about the same as one.
%
% p-values: one-sided, (#null >= observed + 1) / (n_perm + 1).
% FDR: BH across the 326 edge sets, SEPARATELY for each component count.
%
% Checkpoints after every edge set; a restart resumes where it stopped.
%
% OUTPUT (in <basel_output_dir>/PermutationTest):
%   BASEL_PermutationTest_LSDvsPlacebo_Summary_1comp.csv   (one per ncomp)
%   BASEL_PermutationTest_LSDvsPlacebo_Summary_2comp.csv
%   BASEL_PermutationTest_LSDvsPlacebo_Summary_3comp.csv
%   BASEL_PermutationTest_LSDvsPlacebo_ByNComp_1-2-3comp.csv     (long, NComp column)
%   BASEL_PermutationTest_LSDvsPlacebo_NCompComparison_1-2-3comp.csv (wide)
%   BASEL_PermutationTest_LSDvsPlacebo_1-2-3comp.mat         (incl. full null distributions)
% =========================================================================
n_perm       = 5000;
K            = 10;
cv_seed      = 2026;     % fixed partition
perm_seed    = 42;       % fixed label swaps
ncomp_list   = [1 2 3];
n_nc         = numel(ncomp_list);
ncomp_str    = strjoin(arrayfun(@num2str, ncomp_list, 'UniformOutput', false), '-');   % '1-2-3'
model_labels = cellstr(model_labels);
output_dir   = fullfile(basel_output_dir, 'PermutationTest');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
checkpoint_path = fullfile(output_dir, ['BASEL_PermutationTest_LSDvsPlacebo_checkpoint_' ncomp_str 'comp.mat']);

num_edge_sets = numel(all_net_edges);
num_subs      = numel(label_vector);

%% === Fixed subject-level, STUDY-stratified partition ===
unique_subjects    = unique(subject_id_vector, 'stable');
num_unique_subs    = numel(unique_subjects);
site_labels_unique = strings(num_unique_subs, 1);
site_labels_unique(startsWith(unique_subjects, "5-HT_")) = "5-HT";
site_labels_unique(startsWith(unique_subjects, "LAM_"))  = "LAM";
site_labels_unique(startsWith(unique_subjects, "LPM_"))  = "LPM";
if any(site_labels_unique == "")
    error('Unrecognized study prefix for subject(s): %s', strjoin(unique_subjects(site_labels_unique == ""), ', '));
end
fprintf('Study breakdown: 5-HT=%d | LAM=%d | LPM=%d (total %d)\n', ...
    sum(site_labels_unique == "5-HT"), sum(site_labels_unique == "LAM"), ...
    sum(site_labels_unique == "LPM"), num_unique_subs);

rng(cv_seed);
cv_partition  = cvpartition(site_labels_unique, 'KFold', K, 'Stratify', true);
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
Y_all = [label_vector(:), double(xor(label_vector(:), flips(subj_row, :)))];   % col 1 = true labels

%% === Containers: [edge set x permutation x ncomp] (or resume) ===
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
n1 = sum(label_vector == 1);   % LSD scans (unchanged by within-subject swaps)
n0 = sum(label_vector == 0);   % Placebo scans
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
    writetable(T, fullfile(output_dir, sprintf('BASEL_PermutationTest_LSDvsPlacebo_Summary_%dcomp.csv', ncomp_list(c))));
    perm_results = [perm_results; T]; %#ok<AGROW>
end
writetable(perm_results, fullfile(output_dir, ['BASEL_PermutationTest_LSDvsPlacebo_ByNComp_' ncomp_str 'comp.csv']));

%% === Wide table: one row per model, ncomp side by side ===
perm_cmp = table((1:num_edge_sets)', model_labels(:), 'VariableNames', {'EdgeSet', 'Model'});
for c = 1:n_nc
    t = sprintf('_%dcomp', ncomp_list(c));
    perm_cmp.(['Dprime' t])        = obs_dprime(:, c);
    perm_cmp.(['Dprime_perm_p' t]) = p_dprime(:, c);
    perm_cmp.(['Dprime_perm_q' t]) = q_dprime(:, c);
end
writetable(perm_cmp, fullfile(output_dir, ['BASEL_PermutationTest_LSDvsPlacebo_NCompComparison_' ncomp_str 'comp.csv']));

save(fullfile(output_dir, ['BASEL_PermutationTest_LSDvsPlacebo_' ncomp_str 'comp.mat']), ...
    'perm_results', 'perm_cmp', 'obs_acc', 'obs_dprime', 'null_acc', 'null_dprime', ...
    'p_acc', 'p_dprime', 'q_acc', 'q_dprime', 'ncomp_list', 'flips', 'n_perm', 'K', ...
    'cv_seed', 'perm_seed', 'train_idx_all', 'test_idx_all', 'model_labels', ...
    'subject_id_vector', '-v7.3');

fprintf('\n=== Permutation test by number of PLS components ===\n');
for c = 1:n_nc
    fprintf('%d comp: %3d of %d edge sets significant (d'', FDR q < .05) | mean d'' = %.2f\n', ...
        ncomp_list(c), sum(q_dprime(:, c) < .05), num_edge_sets, mean(obs_dprime(:, c)));
end
fprintf('Saved to: %s\n', output_dir);