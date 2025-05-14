
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

% Split the time-series at each transition
start = 1;
data = table;
data.start = {};
data.behaveState = {};
data.lfp = {};
for segment = 1 : numel(behaveState_change)
	row = {...
		start,...
		behaveState(start),...
		lfp(start:behaveState_change(segment))...
		};
	data = [data; row];
	start = behaveState_change(segment)+1;
end

%% Split each segment into smaller time-series

tlength = 10000; % 10 seconds (@1000hz)

data_long = data;

data = table;
data.behaveState = {};
data.lfp = {};

for segment = 1 : size(data_long, 1)
	
	tmp = data_long(segment, :);
	
	% Trim segment so that it can be equally divided
	nEpochs = floor(length(tmp.lfp) / tlength);
	trimmed_length = tlength * nEpochs;
	tmp.lfp = tmp.lfp(1:trimmed_length);

	% Separate segment into smaller epochs
	%	(1 x N) -> (tlength x nEpochs)
	%	test: reshape((1:10), [5 2]) to confirm that time order is
	%		preserved
	tmp.lfp = reshape(tmp.lfp, [tlength nEpochs]);
	
	% Turn epochs into rows of table
	segment_table = table;
	
end