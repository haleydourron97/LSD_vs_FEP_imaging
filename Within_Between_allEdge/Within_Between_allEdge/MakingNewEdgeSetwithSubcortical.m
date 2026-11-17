%% Build edge sets for Schaefer-100 (17 networks) + Tian 38 subcortical

%% === Paths ===
base_dir     = '/Volumes/toshiba/LSD_November11/Within_Between_allEdge';
atlas_dir    = '/Volumes/toshiba/LSD_November11/Within_Between_allEdge/AtlasLabels';
out_dir      = fullfile(base_dir, 'Schaefer100_Tian38');
schaefer_txt = fullfile(atlas_dir, 'Schaefer2018_100Parcels_17Networks_order_info.txt');
tian_txt     = fullfile(atlas_dir, 'Tian_subcortical38_labels.txt');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

n_cort = 100;
n_sub  = 38;
n_roi  = n_cort + n_sub;   % 138

%% === Groups, in figure order ===
group_labels = {'VIS A','VIS B','SMN A','SMN B','DAN A','DAN B','SAL A','SAL B', ...
                'LIM B','LIM A','FPN A','FPN B','FPN C','DN A','DN B','DN C','TPar', ...
                'HIP','AMY','NAc','GP','PUT','CAU','vTHA','dTHA'};
n_groups = numel(group_labels);   % 25

% Schaefer network name (as written in the parcel names) for groups 1-17, same order
cort_net_names = ["VisCent","VisPeri","SomMotA","SomMotB","DorsAttnA","DorsAttnB", ...
                  "SalVentAttnA","SalVentAttnB","LimbicB","LimbicA", ...
                  "ContA","ContB","ContC","DefaultA","DefaultB","DefaultC","TempPar"];

% Tian label patterns for groups 18-25, same order (anchored so each label matches one group)
sub_patterns = {'^[ap]HIP-', '^[lm]AMY-', '^NAc-', '^[ap]GP-', '^[ap]PUT-', '^[ap]CAU-', ...
                '^THA-V', '^THA-D'};

roi_group = zeros(n_roi, 1);
roi_name  = strings(n_roi, 1);

%% === Cortical ROIs: parse Schaefer names ===
lut = strtrim(readlines(schaefer_txt));
lut = lut(lut ~= "");
assert(numel(lut) == 2*n_cort, 'Expected %d lines (name + index per parcel), found %d', 2*n_cort, numel(lut));

for k = 1:n_cort
    name    = lut(2*k - 1);
    numline = sscanf(lut(2*k), '%d');
    % The line after each name starts with the parcel's atlas index: must equal its position
    assert(numline(1) == k, 'Parcel "%s" has atlas index %d but sits at position %d', name, numline(1), k);

    tok = regexp(name, '^17Networks_(LH|RH)_([A-Za-z]+)_', 'tokens', 'once');
    assert(~isempty(tok), 'Could not parse Schaefer name "%s"', name);
    g = find(cort_net_names == tok{2});
    assert(numel(g) == 1, 'Unknown Schaefer network "%s" in "%s"', tok{2}, name);

    roi_group(k) = g;
    roi_name(k)  = name;
end

%% === Subcortical ROIs: parse Tian names ===
tian_names = strtrim(readlines(tian_txt));
tian_names = tian_names(tian_names ~= "");
assert(numel(tian_names) == n_sub, 'Expected 38 Tian labels, found %d', numel(tian_names));

for k = 1:n_sub
    hit = find(cellfun(@(p) ~isempty(regexp(tian_names(k), p, 'once')), sub_patterns));
    assert(numel(hit) == 1, 'Tian label "%s" matched %d groups', tian_names(k), numel(hit));
    roi_group(n_cort + k) = 17 + hit;
    roi_name(n_cort + k)  = tian_names(k);
end

assert(all(roi_group > 0), 'Some ROIs were not assigned to a group');

%% === Check ROI counts per group against the official atlas ===
% Cortical counts from the CBIG Schaefer-100 17-network file; Tian counts from the S2/S3 labels
expected_counts = [7 6 6 8 6 7 9 5 2 3 6 5 5 7 10 4 4, ...   % VIS A ... TPar
                   4 4 4 4 4 4 8 6];                            % HIP AMY NAc GP PUT CAU vTHA dTHA
actual_counts = accumarray(roi_group, 1, [n_groups 1])';
if ~isequal(actual_counts, expected_counts)
    bad = find(actual_counts ~= expected_counts);
    error('ROI counts differ from the official atlas for: %s', ...
        strjoin(arrayfun(@(g) sprintf('%s (%d vs %d expected)', group_labels{g}, ...
        actual_counts(g), expected_counts(g)), bad, 'UniformOutput', false), ', '));
end

% Spot-check the four parcels the old reordering CSV got wrong
assert(roi_group(89) == find(cort_net_names == "DefaultA"), 'ROI 89 should be DefaultA');
assert(roi_group(93) == find(cort_net_names == "DefaultB"), 'ROI 93 should be DefaultB');
assert(roi_group(96) == find(cort_net_names == "DefaultC"), 'ROI 96 should be DefaultC');
assert(roi_group(98) == find(cort_net_names == "TempPar"),  'ROI 98 should be TempPar');

%% === Full connectome edge list ===
[row_idx, col_idx] = find(triu(ones(n_roi), 1));
all_edges = [row_idx, col_idx];            % 9453 x 2
g1 = roi_group(all_edges(:,1));
g2 = roi_group(all_edges(:,2));

%% === Within-group edges ===
within_net_idx = cell(n_groups, 1);
for g = 1:n_groups
    within_net_idx{g} = all_edges(g1 == g & g2 == g, :);
end

%% === Between-group edges (symmetric cell array) ===
btw_net_idx = cell(n_groups, n_groups);
for i = 1:n_groups
    for j = (i+1):n_groups
        is_btw = (g1 == i & g2 == j) | (g1 == j & g2 == i);
        btw_net_idx{i,j} = all_edges(is_btw, :);
        btw_net_idx{j,i} = btw_net_idx{i,j};
    end
end

%% === Edge sets + model labels, built TOGETHER ===
% Order: FullConnectome, Within-<g> (1..25), Between-<i>-<j> (i<j, row-major)
edge_sets    = {all_edges};
model_labels = {'FullConnectome'};
model_type   = {'full'};
model_i      = 0;
model_j      = 0;

for g = 1:n_groups
    edge_sets{end+1,1}    = within_net_idx{g};
    model_labels{end+1,1} = sprintf('Within-%s', group_labels{g});
    model_type{end+1,1}   = 'within';
    model_i(end+1,1) = g;
    model_j(end+1,1) = g;
end

for i = 1:n_groups
    for j = (i+1):n_groups
        edge_sets{end+1,1}    = btw_net_idx{i,j};
        model_labels{end+1,1} = sprintf('Between-%s-%s', group_labels{i}, group_labels{j});
        model_type{end+1,1}   = 'between';
        model_i(end+1,1) = i;
        model_j(end+1,1) = j;
    end
end

%% === Sanity checks ===
n_models = numel(edge_sets);
assert(n_models == 1 + n_groups + n_groups*(n_groups-1)/2, 'Unexpected number of edge sets');
assert(numel(model_labels) == n_models);
assert(all(~cellfun(@isempty, edge_sets)), 'An edge set is empty');

n_within = sum(cellfun(@(x) size(x,1), within_net_idx));
n_btw    = sum(cellfun(@(x) size(x,1), edge_sets(strcmp(model_type, 'between'))));
assert(n_within + n_btw == size(all_edges,1), 'Within + between edges do not cover the connectome');

for m = 2:n_models
    e = edge_sets{m};
    pair = sort([roi_group(e(:,1)), roi_group(e(:,2))], 2);
    assert(all(pair(:,1) == model_i(m) & pair(:,2) == model_j(m)), ...
        'Edge set %d (%s) contains edges from the wrong groups', m, model_labels{m});
end

fprintf('%d ROIs, %d groups, %d edges, %d edge sets (1 full, %d within, %d between)\n', ...
    n_roi, n_groups, size(all_edges,1), n_models, n_groups, n_groups*(n_groups-1)/2);

%% === Print group membership for a visual check ===
for g = 1:n_groups
    idx = find(roi_group == g);
    fprintf('%-6s (%2d ROIs, %3d within-edges): %s\n', group_labels{g}, numel(idx), ...
        size(within_net_idx{g}, 1), mat2str(idx'));
end

%% === Save ===
save(fullfile(out_dir, 'all_edges.mat'),      'all_edges');
save(fullfile(out_dir, 'within_net_idx.mat'), 'within_net_idx');
save(fullfile(out_dir, 'btw_net_idx.mat'),    'btw_net_idx');
save(fullfile(out_dir, 'edge_sets_and_labels.mat'), ...
    'edge_sets', 'model_labels', 'model_type', 'model_i', 'model_j', ...
    'group_labels', 'roi_group', 'roi_name');

writetable(table((1:n_roi)', roi_name, roi_group, group_labels(roi_group)', ...
    'VariableNames', {'ROI','AtlasLabel','GroupIndex','GroupLabel'}), ...
    fullfile(out_dir, 'roi_group_assignment.csv'));

writetable(table((1:n_models)', model_labels, model_type, model_i, model_j, ...
    cellfun(@(x) size(x,1), edge_sets), ...
    'VariableNames', {'EdgeSetIndex','Model','Type','GroupI','GroupJ','NumEdges'}), ...
    fullfile(out_dir, 'edge_set_labels.csv'));

fprintf('Saved to %s\n', out_dir);