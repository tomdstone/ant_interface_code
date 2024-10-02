function [ EEG ] = ANT_interface_loadset(filename, filepath, verbose, todouble)
%
% ANT INTERFACE CODES - LOADSET
%
% - used to load an EEGLAB format .set file containing the EEG structure
% with the data and other recording information.
%
% Last edit: Alex He 05/04/2024
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - filename:     file name of the .set file.
%
%           - filepath:     full path to the folder containing the .set
%                           file
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
%           - todouble      whether convert EEG.data to double from single.
%                           default: false
%
% Output:
%           - EEG:          an EEGLAB structure containing all information
%                           of the recording in .cnt file. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin < 3
    verbose = true;
    todouble = false;
elseif nargin < 4
    todouble = false;
end

% addpath to the appropriate folders
try 
    SleepEEG_addpath(matlabroot);
    
catch
    % if using SleepEEG_addpath() fails, we will assume the current directory
    % has the ANT_interface_loadset.m or at least the folder containing it has
    % been added to path when calling this function. We will try to addpath to
    % EEGLAB directly.

    ANTinterface_path = which('ANT_interface_loadset');
    temp = strsplit(ANTinterface_path, 'ANT_interface_loadset.m');
    
    % Add path to EEGLAB
    addpath(fullfile(temp{1}, 'eeglab14_1_2b'))
end

%% Load data
% Start EEGLab
eeglab nogui;

if verbose; tic; end
% Call pop_loadset.m function from EEGLAB to load .set file
EEG = pop_loadset(filename, filepath);
if verbose
    disp(' ')
    disp('Total time taken in Loading the dataset...')
    disp(' ')
    toc
end

%% Change EEG data from single to double
if todouble && ~isa(EEG.data, 'double')
    EEG.data = double(EEG.data);
end

end
