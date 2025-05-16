
%lfp: local field potential from the sFL, recorded at 1000 Hz, in uV
%skinBrightness: mean mantle skin brightness (z scored)
%behaveState: vector defining behavioral state (manually categorized)
            % 0: quiet sleep
            % 1: wake
            % 2: active sleep
            
% data example has an active sleep bout surrounded by quiet sleep, 
% ending with the animal waking up. 
% Note: similarity of wake-active sleep activity
%       presence of  oscillatory burst events during quiet sleep 
%       wake-like activity during brief color blashes during quiet sleep

lfp=h5read('exampleEphys.h5','/lfp');
skinBrightness=h5read('exampleEphys.h5','/skinBrightness');
behaveState=h5read('exampleEphys.h5','/behaveState');

%% Plot the data

figure;
set(gcf, 'Color', 'w');
hold on
plot(zscore(lfp))
plot(skinBrightness,'r','linewidth',3)
plot(behaveState,'k','linewidth',3)

%% Split the time-series based on behaveState

states = [0 1 2];
state_labels = {'QS', 'W', 'AS'}; % quit sleep, wake, active sleep

% Identify state with least amount of data
nSamples_perState = nan(size(states));
for s = 1 : length(states)
    nSamples_perState(s) = sum(behaveState==(states(s)));
end

% Find transitions between states
behaveState_change = find(diff(behaveState) ~= 0);
% Mark the end of the recording as a change
%	Otherwise, the last contiguous state gets dropped
behaveState_change = [behaveState_change size(lfp, 2)];
behaveState_start = [1 behaveState_change(1:end-1)+1];

% Split the time-series at each transition
field_names = {'start', 'behaveState', 'lfp'};
field_types = {'double', 'double', 'cell'};
data_long = table(...
	'Size', [numel(behaveState_change) numel(field_names)],...
	'VariableNames', field_names,...
	'VariableTypes', field_types);
for segment = 1 : numel(behaveState_change)
	
	start = behaveState_start(segment);
	
	data_long.start(segment) = start;
	data_long.behaveState(segment) = behaveState(start);
	data_long.lfp(segment) = {lfp(start:behaveState_change(segment))};
	
end

%% Split each segment into smaller time-series

tlength = 10000; % 10 seconds (@1000hz)

% Figure out how many total epochs there will be
nEpochs_total = 0;
nEpochs_perSegment = nan(size(behaveState_start));
for s = 1 : size(data_long, 1)
	
	tmp = data_long(s, :);
	
	nEpochs = floor(length(tmp.lfp{1}) / tlength);
	nEpochs_total = nEpochs_total + nEpochs;
	nEpochs_perSegment(s) = nEpochs;
end

% Preallocate table
field_names = {'start', 'behaveState', 'epoch', 'lfp'};
field_types = {'double', 'double', 'double', 'cell'};
data = table(...
	'Size', [nEpochs_total numel(field_names)],...
	'VariableNames', field_names,...
	'VariableTypes', field_types);

row_counter = 1;
for segment = 1 : size(data_long, 1)
	
	tmp = data_long(segment, :);
	
	% Trim segment so that it can be equally divided
	nEpochs = floor(length(tmp.lfp{1}) / tlength);
	trimmed_length = tlength * nEpochs;
	tmp.lfp{1} = tmp.lfp{1}(1:trimmed_length);

	% Separate segment into smaller epochs
	%	(1 x N) -> (tlength x nEpochs)
	%	test: reshape((1:10), [5 2]) to confirm that time order is
	%		preserved
	tmp.lfp{1} = reshape(tmp.lfp{1}, [tlength nEpochs]);
	
	% Turn epochs into rows of table
	segment_table = table(...
		'Size', [nEpochs numel(field_names)],...
		'VariableNames', field_names,...
		'VariableTypes', field_types);
	
	for e = 1 : size(tmp.lfp{1}, 2)
		
		segment_table.start(e) = tmp.start;
		segment_table.behaveState(e) = tmp.behaveState;
		segment_table.epoch(e) = e;
		segment_table.lfp{e} = tmp.lfp{1}(:, e);
		
	end
	
	% Add segment table to overall table
	data(row_counter : row_counter + size(segment_table, 1) - 1, :) = segment_table;
	
	row_counter = row_counter + size(segment_table, 1);
	
end

%% Figure out which segment is the shortest
% Take the same number of epochs from each segment to analyse

sample_nEpochs = min(nEpochs_perSegment);

% Preallocate table
field_names = {'start', 'behaveState', 'epoch', 'lfp'};
field_types = {'double', 'double', 'double', 'cell'};
data_sampled = table(...
	'Size', [sample_nEpochs*numel(behaveState_start) numel(field_names)],...
	'VariableNames', field_names,...
	'VariableTypes', field_types);

% Take the middle N epochs of each segment
row_counter = 1;
for s = 1 : numel(behaveState_start)
	
	segment_indices = find(data.start == behaveState_start(s));
	midpoint = ceil(length(segment_indices)/2);
	
	data_sampled(row_counter:row_counter+sample_nEpochs-1, :) = data(segment_indices(midpoint-2:midpoint+2), :);
	
	row_counter = row_counter + sample_nEpochs;
	
	% Track where in the overall LFP the sampling is from
	line(data.start([segment_indices(midpoint) segment_indices(midpoint)]) + tlength*data.epoch(segment_indices(midpoint)), [-25 20],...
		'Color', 'k', 'LineWidth', 1);
end

%%

axis tight

sampleRate = 1000;

tstep = sampleRate*60*5;
taxis = (1:tstep:length(lfp));
taxis_label = round(taxis / (sampleRate*60));

set(gca, 'XTick', taxis, 'XTickLabel', taxis_label);

xlabel('t (mins)');

legend('LFP (z)', 'skinBrightness', 'behaveState', 'Location', 'south');

%% Setup structures for hctsa

% timeSeriesData: Nx1 cell, or NxM matrix
% labels: Nx1 cell of unique strings
% keywords: Nx1 cell of comma-delimited keywords as a string

timeSeriesData = data_sampled.lfp;

% Generate labels
labels = cell(size(data_sampled, 1), 1);
keywords = cell(size(labels));
for e = 1 : size(data_sampled, 1)
	
	label = [...
		state_labels{data_sampled.behaveState(e)+1} ','...
		'behaveState' num2str(data_sampled.behaveState(e)) ','...
		'epoch' num2str(data_sampled.epoch(e)) ','...
		'start' num2str(data_sampled.start(e))...
		];
	
	labels{e} = label;
	keywords{e} = label;
	
end

%% Save variables to .mat

init_save = 1;

init_dir = 'hctsa_space/';
init_file = 'example_hctsa_init.mat';
if init_save == 1
	
	if ~isfolder(init_dir)
		mkdir(init_dir);
	end
	
	save([init_dir init_file], 'timeSeriesData', 'labels', 'keywords');
	
end

%% Run hctsa initialisation

init = 1;

hctsa_file = 'example_hctsa.mat';

if init == 1
	TS_Init([init_dir init_file], [], [], [false false false], [init_dir hctsa_file]);
end

%% Compute feature values

compute = 1;

if compute == 1
	TS_Compute(true, [], [], [], [init_dir hctsa_file]);
end
