
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

figure
hold on
plot(zscore(lfp))
plot(skinBrightness,'r','linewidth',3)
plot(behaveState,'k','linewidth',3)