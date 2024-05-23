 function [ EEG ] = ANT_interface_readcnt(filename, filepath, dsrate, verbose)
%
% ANT INTERFACE CODES - READCNT
%
% - this function is used to read a .cnt file exported from the ANT eego
% mylab system (2019). It expects a .cnt file, a .evt file, and a .seg file
% (if the recording was broken into multiple segments) with the same naming
% before the file extension.
%
% Last edit: Alex He 05/22/2024
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - filename:     file name of the .cnt file, the .evt and .seg
%                           files must have the same naming before the
%                           extension.
%
%           - filepath:     full path to the folder containing the .cnt
%                           file (and .evt, .seg files).
%
%           - dsrate:       whether to downsample and desired sampling
%                           rate. This is useful when the recording is
%                           very long such that we can't keep a large data
%                           structure in the original sampling rate in
%                           cache memory.
%                           default: [false, 0]
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
% Output:
%           - EEG:          an EEGLAB structure containing all information
%                           of the recording in .cnt file.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% set default values of flag arguments
if nargin < 3
    dsrate = [false, 0];
    verbose = true;
elseif nargin < 4
    verbose = true;
end

%% Addpath to required toolbox folders
% Note: this section expects these two folders to exist in the same
% directory as the ANT_interface_readcnt.m function. If not, you would have
% to manually addpath to the folder containing these two toolbox folders,
% or you can modify the directories to be full paths to your desired
% locations.

% Let's first attempt to use the butter function SleepEEG_addpath(). You
% should be able to call it if you are using ANT_interface_readcnt() from
% any of the SleepEEG_ functions.

try
    SleepEEG_addpath(matlabroot);
catch
    % if using SleepEEG_addpath() fails, we will assume the current directory
    % has the ANT_interface_readcnt.m or at least the folder containing it has
    % been added to path when calling this function. We will try to addpath to
    % ANT importer and EEGLAB directly.
    
    ANTinterface_path = which('ANT_interface_readcnt');
    temp = strsplit(ANTinterface_path, 'ANT_interface_readcnt.m');
    
    % Add path to ANT importer
    addpath(fullfile(temp{1}, 'ANTeepimport1.13'))
    
    % Add path to EEGLAB
    addpath(fullfile(temp{1}, 'eeglab14_1_2b'))
    
    % Add path to EEGLAB firfilt plugin for pop_resample()
    addpath(fullfile(temp{1}, 'eeglab14_1_2b/plugins/firfilt1.6.2'))
end

%% Start EEGLab
eeglab nogui; close;

%% Construct data filename
datafn = fullfile(filepath, filename);

%% Report EEG record length
if verbose
    disp(' ')
    if dsrate(1)
        disp('Loading and Downsampling...')
    else
        disp('Loading data file...')
    end
    disp(' ')
    
    if ~isfile(datafn) % if .cnt file not on disk report error
        error('Cannot locate the .cnt file. Please check the path and filename :\n%s.',datafn)
    else
        disp('Loading .cnt file from:')
        disp(datafn)
    end
end

cnt_info = eepv4_read_info(datafn); % this function is case-insensitive!
tot_sample_point = cnt_info.sample_count;
record_time = tot_sample_point/cnt_info.sample_rate/60;

if verbose
    disp(['Total recording duration = ' num2str(record_time) ' minutes = ', num2str(record_time/60) ' hours.'])
    disp(['Total sample points = ' num2str(cnt_info.sample_count) '.'])
    disp(' ')
end

%% Load data in pieces
% construct 30-min segment sample points for looping
sample_point = 0:cnt_info.sample_rate*60*30:cnt_info.sample_count;
sample_point = [sample_point cnt_info.sample_count]; % add in the last sample point
if verbose
    if length(sample_point) > 2
        disp(['A total of ' num2str(length(sample_point)-1) ' segments of 30min'])
    else
        disp(['A total of ' num2str(length(sample_point)-1) ' segment of ' num2str(record_time) 'min'])
    end
end

% load each segment sequentially and save each segment into a cell store
EEG_store = cell(1, length(sample_point)-1);

if verbose; tic; end
for i = 1:length(sample_point)-1
    % load in one segment of data
    if i == 1
        sample1 = 1;
    else
        sample1 = sample_point(i)+1;
    end
    sample2 = sample_point(i+1);
    if verbose
        disp(' ')
        disp(['Loading the ' num2str(i) 'th segment...'])
        disp(['sample1: ', num2str(sample1)])
        disp(['sample2: ', num2str(sample2)])
    end
    
    % use pop_loadeep_v4.m function to load in the segment
    EEG = pop_loadeep_v4(datafn, 'sample1', sample1, 'sample2', sample2);
    
    % test whether export montage and filtering was applied
    if i == 1 % we only need to check the first segment to raise the error
        [pxx,f] = pwelch(EEG.data(1,:),[],[],[],EEG.srate);
        % if filtering was applied, the high frequency content will be less
        % than -100dB! This likely implies that post-hoc montage and
        % filtering was applied when exporting the raw data from eego
        % software into .cnt and .evt files. Please re-export!
        pxx_db_high = pow2db(pxx(f>=(EEG.srate/2 * 0.8))); % check in the 400-500Hz range if srate is 1000Hz
        mean_db_high = mean(pxx_db_high(isfinite(pxx_db_high)));
        assert(mean_db_high >= -100, 'The high frequency power in LM is abnormally low. Do not apply montages and filtering when exporting raw data in eego!')
    end
    
    if dsrate(1) % if downsampling flag is true
        % For using pop_resample function: trials number is updated to be 1
        if EEG.trials == 0
            EEG.trials = 1;
        end
        
        % downsample to desired sampling rate in Hz as specified by dsrate(2)
        % pop_resample function
        % Inputs:
        %   INEEG      - input dataset
        %   freq       - frequency to resample (Hz)
        %
        % Optional inputs:
        %   fc         - anti-aliasing filter cutoff (pi rad / sample)
        %                {default 0.9}
        %   df         - anti-aliasing filter transition band width (pi rad /
        %                sample) {default 0.2}        
        if verbose; disp(' '); disp('Downsampling...'); end
        EEG = pop_resample(EEG, dsrate(2), 0.9, 0.1);
        
        % When EEG.event is resampled, there could be non-integer latency
        % values. We write a for loop to convert them to next larger
        % integer using the ceil() function.
        for j = 1:length(EEG.event)
            EEG.event(j).latency = ceil(EEG.event(j).latency);
        end
        
    else
        if verbose
            disp(' ')
            disp('Downsampling is turned off. Original sampling rate will be retained.')
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end
    
    % save the downsampled data to a temporary cell
    EEG_store{i} = EEG;
    
    % clearvars to free up memory
    clearvars EEG
end
if verbose
    disp(' ')
    if dsrate(1)
        disp('Total time taken in Loading and Downsampling...')
    else
        disp('Total time taken in Loading...')
    end
    disp(' ')
    toc
end

%% Cascade downsampled data
if verbose
    disp(' ')
    if dsrate(1)
        disp('Concatenating downsampled data...')
    else
        disp('Concatenating data in original sampling rate...')
    end
    disp(' ')
end

if verbose; tic; end
% Pick up the first EEG structure
EEG = EEG_store{1};
if verbose; disp('Concatenating the 1st segment...'); end
% Loop through the rest EEG structures
for i = 2:length(EEG_store) % starts from 2nd index
    EEG = ANT_interface_catEEG(EEG, EEG_store{i});
    if verbose; disp(['Concatenating the ' num2str(i) 'th segment...']); end
end
if verbose
    disp('Total time taken in Concatenating Data...')
    disp(' ')
    toc
    disp(' ')
end

% Clear EEG_store to free up memory
clearvars EEG_store

%% Set montage
% This function does two things:
%   - adds an empty reference channel back
%   - fills in the channel location info from a template
EEG = ANT_interface_setmontage(EEG, 'auto');

end
