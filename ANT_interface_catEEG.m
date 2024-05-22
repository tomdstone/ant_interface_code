function [EEG] = ANT_interface_catEEG(EEG, EEG_new)
%
% ANT INTERFACE CODES - CATEEG
%
% % - function to concatenate a segment EEG structure to the master EEG
% structure. This is used in reading .cnt files as segments of 30min for
% long recordings. This function also corrects sample latencies of event
% logs to make them accurate in the resultant master EEG structure. 
%
% Last edit: Alex He 05/21/2024

%%
% Confirm that the sampling rates are identical and the same original files
assert(EEG.srate == EEG_new.srate, 'Different sampling rate, please check!')
assert(strcmp(EEG.comments, EEG_new.comments), 'Different original file!')
assert(size(EEG.data, 1) == size(EEG_new.data, 1), 'Different numbers of channels!')

% % Manually updating each field in the EEG structure to concatenate the
% % second segment encoded in EEG_new
% for i = 1:length(EEG_new.event)
%     EEG_new.event(i).latency = EEG_new.event(i).latency + EEG.pnts;
%     EEG.event(length(EEG.event)+1) = EEG_new.event(i);
% end
% 
% EEG.pnts = EEG.pnts + EEG_new.pnts;
% EEG.data = [EEG.data, EEG_new.data];
% if EEG_new.xmin == 0
%     EEG.times = [EEG.times, EEG_new.times + EEG.times(end) + 1/EEG.srate*1000];
% else
%     EEG.times = [EEG.times, EEG_new.times];
% end
% EEG.xmax = EEG_new.xmax;

% Use EEGLAB function to merge the structures
EEG = pop_mergeset(EEG, EEG_new);

% Check if impedance measures should be incorporated
if ~isempty(EEG_new.endimp)
    EEG.endimp = EEG_new.endimp;
end
if ~isempty(EEG_new.initimp)
    EEG.initimp = EEG_new.initimp;
end

end
