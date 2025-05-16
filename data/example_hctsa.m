%%

source_file = 'hctsa_space/example_hctsa.mat';

hctsa = load(source_file);

%% Scale hctsa values

hctsa_mat = BF_NormalizeMatrix(hctsa.TS_DataMat, 'mixedSigmoid');

%% Cluster columns

fCorr = corr(hctsa_mat, 'Type', 'Spearman', 'rows', 'pairwise');

invalid_features = all(isnan(fCorr), 1);

fCorr(invalid_features, :) = [];
fCorr(:, invalid_features) = [];
hctsa_mat(:, invalid_features) = [];

distances = 1 - fCorr;

tree = linkage(squareform(distances), 'average'); % note - distances must be pdist vector (treats matrix as data instead of distances)

% Sorted features
f = figure('visible', 'on'); % we want the order, not the actual plot
[h, T, order] = dendrogram(tree, 0);
close(f);

%%

figure;

imagesc(fCorr(order, order));

%%

figure;
set(gcf, 'Color', 'w');

imagesc(hctsa_mat(:, order));
c = colorbar;

set(gca, 'YTick', (1:5:20), 'YTickLabel', {'QS', 'AS', 'QS', 'W'});
ylabel('epoch');

xlabel('feature');

title(c, 'norm. value');