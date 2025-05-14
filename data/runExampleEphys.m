%% create sample data
% This code was used to create the example data exampleEphys.h5

cd('/home/sam/bucket/octopus/electrophysiology/2022_03_07')
load('rec3_g0_selectChannel.mat')
load('rec3_g0_behaviorAnalysis.mat')

currT=1.5*10^7:1.9*10^7;
lfp=allLFP(2,currT);
skinBrightness=zscore(mi(currT));
behaveState=wakeVec(currT);

h5create('exampleEphys.h5','/lfp',size(lfp))
h5write('exampleEphys.h5','/lfp',lfp)
h5create('exampleEphys.h5','/skinBrightness',size(skinBrightness))
h5write('exampleEphys.h5','/skinBrightness',skinBrightness)
h5create('exampleEphys.h5','/behaveState',size(behaveState))
h5write('exampleEphys.h5','/behaveState',behaveState)

%% look at sample data
% This code is an example for looking at the example data exampleEphys.h5

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

pathToExampleData='/home/sam/pophale2023/example_data';
cd(pathToExampleData)
lfp=h5read('exampleEphys.h5','/lfp');
skinBrightness=h5read('exampleEphys.h5','/skinBrightness');
behaveState=h5read('exampleEphys.h5','/behaveState');

figure
hold on
plot(zscore(lfp))
plot(skinBrightness,'r','linewidth',3)
plot(behaveState,'k','linewidth',3)