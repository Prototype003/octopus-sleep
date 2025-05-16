%% Load

source_file = 'hctsa_space/example_hctsa.mat';

hctsa = load(source_file);

%%

% Keep track of features which are removed
feature_list = (1:size(hctsa.TS_DataMat, 2));

% Condition labels for each epoch
tmp = cellfun(@(x) strsplit(x, ','), hctsa.TimeSeries.Name, 'UniformOutput', false);
epoch_labels = cellfun(@(x) x{1}, tmp, 'UniformOutput', false);

% Colours for each condition
cLegend = dictionary(); % [0 1 2] = QS, W, AS
cLegend('QS') = 1; % blue
cLegend('W') = 2; % red
cLegend('AS') = 3; % purple

cmap = brewermap(3, 'set1'); % set1 = red, blue, green
cmap = cmap([2 1 3], :); % reorder colours

%% Scale hctsa values

hctsa_mat = BF_NormalizeMatrix(hctsa.TS_DataMat, 'mixedSigmoid');

%% Cluster features (columns)

fCorr = corr(hctsa_mat, 'Type', 'Spearman', 'rows', 'pairwise');

invalid_features = all(isnan(fCorr), 1);

fCorr(invalid_features, :) = [];
fCorr(:, invalid_features) = [];
hctsa_mat(:, invalid_features) = [];
feature_list(:, invalid_features) = [];

fDistances = 1 - fCorr;

fTree = linkage(squareform(fDistances), 'average'); % note - distances must be pdist vector (treats matrix as data instead of distances)

% Sorted features
f = figure('visible', 'on'); % we want the order, not the actual plot
[h, T, fOrder] = dendrogram(fTree, 0);
close(f);

%%

figure;

imagesc(fCorr(fOrder, fOrder));

%% Show hctsa matrix with clustered features

figure;
set(gcf, 'Color', 'w');

imagesc(hctsa_mat(:, fOrder));
c = colorbar;

set(gca, 'YTick', (1:5:20), 'YTickLabel', {'QS', 'AS', 'QS', 'W'});
ylabel('epoch');

xlabel('feature');

title(c, 'norm. value');

%% Identify "top" features to focus on (for clustering)

topN = 40;

% Conduct signrank/ranksum test for each feature
% Compare wake and quiet sleep
compare_states = [1 2];
compare_vals = cell(size(compare_states));
for s = 1 : numel(compare_states)
	rows = cLegend(epoch_labels) == compare_states(s);

	compare_vals{s} = hctsa_mat(rows, :);
end

% Test
f_zStats = nan(1, size(hctsa_mat, 2));
for f = 1 : size(hctsa_mat, 2)

	[p, h, stats] = ranksum(compare_vals{1}(:, f), compare_vals{2}(:, f));
	f_zStats(f) = stats.ranksum;

	%[p, h, ci, stats] = ttest2(compare_vals{1}(:, f), compare_vals{2}(:, f));
	%f_zStats(f) = stats.tstat;

end

% Take N features with the greatest test stats
[sorted, sort_order] = sort(f_zStats, 'descend');
topFeatures = feature_list(sort_order(1:topN));

% Get names of the top features
topFeatures_names = hctsa.Operations.Name(topFeatures);

%% Cluster only the top features

fCorr = fCorr(sort_order(1:topN), sort_order(1:topN));

fDistances = 1 - fCorr;

fTree = linkage(squareform(fDistances), 'average'); % note - distances must be pdist vector (treats matrix as data instead of distances)

%%

% Sorted features
f = figure('color', 'w', 'visible', 'on'); % we want the order, not the actual plot
[h, T, fOrder] = dendrogram(fTree, 0);

view([270, 270]);

ylabel('dist');
set(gca, 'XTick', (1:topN), 'XTickLabel', topFeatures_names(fOrder), 'TickLabelInterpreter', 'none');
set(gca,'xaxisLocation','top')

box on;

%% Plot

figure('color', 'w');

imagesc(fCorr(fOrder, fOrder));
c = colorbar;
title(c, 'r');

title('feature correlation');

xlabel('feature');
ylabel('feature');

axis square

%% Cluster epochs (rows)

eCorr = corr(hctsa_mat', 'Type', 'Spearman');

eDistances = 1 - eCorr;

eTree = linkage(squareform(eDistances), 'average'); % note - distances must be pdist vector (treats matrix as data instead of distances)

% Sorted features
f = figure('color', 'w', 'visible', 'on'); % we want the order, not the actual plot
subplot(5, 6, [1 25]);
[h, T, eOrder] = dendrogram(eTree, 0);
view([270, 90]);
%close(f);

%% Dendrogram labels

ylabel('dist');
set(gca, 'XTick', (1:numel(epoch_labels)), 'XTickLabel', epoch_labels(eOrder));

box on;

%% Plot epoch correlations

subplot(5, 6, [2 30]);

imagesc(eCorr(eOrder, eOrder));
title('epoch correlation');

c = colorbar;
title(c, 'r');

set(gca, 'XTick', (1:numel(epoch_labels)), 'XTickLabel', epoch_labels(eOrder));
set(gca, 'YTick', (1:numel(epoch_labels)), 'YTickLabel', epoch_labels(eOrder));

axis square

%% MDS across epochs

nDims = 2;

Y = cmdscale(eDistances, nDims);

%% Show MDS

figure('color', 'w');
%scatter(Y(:, 1), Y(:, 2), [], cmap(cLegend(epoch_labels), :), 'filled');
classes = cLegend(epoch_labels);

hold on;
for c = 1 : numel(cLegend.keys)
	
	plot_points = classes == c;

	h(c) = scatter(Y(plot_points, 1), Y(plot_points, 2), 100, cmap(cLegend(epoch_labels(plot_points)), :), 'filled');

end
legend(h, cLegend.keys);

title('cMDS');

xlabel('dim1');
ylabel('dim2');