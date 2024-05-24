function [ EEG ] = ANT_interface_setmontage(EEG, montage)
%
% ANT INTERFACE CODES - SETMONTAGE
%
% - a new function to automatically detect which cap montage was used for
% the recording and to **load** the correct template channel locations. An
% empty reference channel is also added to the data to preserve rank in
% subsequent processing such as re-referencing and for forward modeling.
%
% This function now replaces sections calling ANT_interface_eego2asa(),
% loading the Duke template channel location (chanlocs), and updating the
% EEG structure .ref and .refscheme fields in ANT_interface_readcnt().
%
% A key change is that the recording channel order is kept intact. Instead,
% corresponding coordinates are grabbed from the template chanlocs. This is
% so that the exported data can be handled similarly either as .cnt file
% and read with ANT_interface_readcnt(), or as BrainVision files and read
% directly into MNE-Python.
%
% Last edit: Alex He 05/22/2024
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - EEG:          an EEG structure with EEG.data in the order of
%                           eego mylab EEG system used to collect the data.
%           - montage:      a string specifying which montage to use,
%                           default to 'auto' that automatically detects
%                           the montage used.
%
% Output:
%           - EEG:          an EEG structure with EEG.data in the same
%                           order as the input structure, with an empty
%                           reference channel added to the end, and with
%                           correct channel location info filled in the
%                           .chanlocs field.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Make sure the correct number of channels is read from the .cnt file
assert(EEG.nbchan == length(EEG.chanlocs), 'Inconsistent number of channels read from .cnt file.')
assert(EEG.nbchan == 128, 'Montage setting is only supported for the 128-channel caps for now.')

if nargin < 2
    montage = 'auto';
end

%% Automatically detect the montage used
if strcmp(montage, 'auto')
    % Extract the current channel labels
    labels = {EEG.chanlocs.labels};
    
    % Auto-detect among waveguard montages
    if strcmp(labels{1}, 'Lm') && ...
            strcmp(labels{85}, 'VEOGL') && ...
            ~any(cellfun(@(x) strcmp(x, 'Z3'), labels)) && ...
            ~any(cellfun(@(x) strcmp(x, 'LL14'), labels))
        montage = 'gelDuke-Z3';
        
    elseif strcmp(labels{1}, 'Z2') && ...
            strcmp(labels{66}, 'LL14') && ...
            ~any(cellfun(@(x) strcmp(x, 'Z7'), labels)) && ...
            ~any(cellfun(@(x) strcmp(x, 'VEOGL'), labels))
        montage = 'salineNet-Z7';
        
    else
        error('Could not detect a compatible montage. Please manually check!')
        
    end
end

%% Add the reference channel back to the data
% Move the extra channels (such as bipolar) if there is any
if EEG.nbchan >= 129
    EEG.data(130:EEG.nbchan+1, :) = EEG.data(129:end, :);
    EEG.chanlocs(130:EEG.nbchan+1) = EEG.chanlocs(129:end);
    if ~isempty(EEG.initimp)
        EEG.initimp(130:EEG.nbchan+1) = EEG.initimp(129:end);
    end
    if ~isempty(EEG.endimp)
        EEG.endimp(130:EEG.nbchan+1) = EEG.endimp(129:end);
    end
end

% Insert a reference channel as zero recording
EEG.data(129, :) = zeros(1, EEG.pnts);

% Update the reference field
EEG.ref = 'see refscheme';

% Add the reference channel to chanlocs
switch montage
    case 'gelDuke-Z3'
        EEG.chanlocs(129).labels = 'Z3';
        EEG.refscheme = 'Z3';
    case 'salineNet-Z7'
        EEG.chanlocs(129).labels = 'Z7';
        EEG.refscheme = 'Z7';
end

if ~isempty(EEG.initimp)
    EEG.initimp(129) = NaN;
end

if ~isempty(EEG.endimp)
    EEG.endimp(129) = NaN;
end

% Update the number of channels
EEG.nbchan = size(EEG.data, 1);

%% Fill in the channel location info from template
% this is not the digitized channel location for each individual, just the
% channel coordinates from montage templates

% Load an appropriate montage template
switch montage
    case 'gelDuke-Z3'
        chanlocs = load('ANT_montage_templates.mat', 'chanlocs_dukeZ3');
        chanlocs = chanlocs.chanlocs_dukeZ3;
    case 'salineNet-Z7'
        chanlocs = load('ANT_montage_templates.mat', 'chanlocs_netZ7');
        chanlocs = chanlocs.chanlocs_netZ7;
end
labels = {chanlocs.labels};

% Create an empty chanlocs using the same fields as the template chanlocs
fns = fieldnames(chanlocs)';
fns{2, 1} = {};
tmp_chanlocs = struct(fns{:});

for ii = 1:EEG.nbchan
    current_label = EEG.chanlocs(ii).labels;
    template_idx = find(cellfun(@(x) strcmp(x, current_label), labels));
    if isempty(template_idx) % channel label not found in template
        tmp_chanlocs(ii).labels = current_label;
    else
        tmp_chanlocs(ii) = chanlocs(template_idx);
    end
end

% Update EEG.chanlocs with location info populated from an template
EEG.chanlocs = tmp_chanlocs;

end
