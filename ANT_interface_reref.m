function [ EEG, data ] = ANT_interface_reref(data, ref_method, G, EEG, verbose)
%
% ANT INTERFACE CODES - REREF
%
% - function to re-reference EEG recording.
%            - AR (common average)
%            - Z3 (recording reference at Z3)
%            - [CH] (arbitrary channel number as reference)
%            - left mastoid
%            - linked mastoid
%            - REST (reference electrode standardization
%                    technique)
%            - contral mastoid
%            - LP (Laplacian reference based on duke layout)
%
% Last edit: Alex He 05/22/2024
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - data:         a double matrix containing recording data.
%
%           - ref_method:   re-reference method.
%                              - AR (common average) [default]
%                              - [CH] (arbitrary channel number as unipolar reference)
%                              - LP (Laplacian reference)
%                              - REST (reference electrode standardization technique)
%                              - rREST (regularized REST)
%
%           - G:            lead field matrix for REST-referencing. If not
%                           provided, the lead field matrix from example
%                           data will be used for REST referencing.
%                           default: []
%
%           - EEG:          data structure containing the EEG data and
%                           reference scheme specification.
%                           default: empty structure
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
% Output:
%           - EEG:          data structure containing the updated
%                           reference scheme specification and data.
%
%           - data:         a double matrix containing re-referenced data.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin < 2 % for simple use on raw data provided as first input
    ref_method = 'AR';
    G = [];
    EEG = struct; EEG.refscheme = 'NaN';
    verbose = true;
else
    if ~exist('G', 'var')
        G = [];
    end
    if ~exist('EEG', 'var')
        EEG = struct; EEG.refscheme = 'NaN';
    end
    if ~exist('verbose', 'var')
        verbose = true;
    end
end
if isempty(data)
    justEEG = true;
else
    justEEG = false;
end

%% Load the reference matrix
ANTinterface_path = which('ANT_interface_reref');
temp = strsplit(ANTinterface_path, 'ANT_interface_reref.m');
ANTinterface_path = temp{1};
load(fullfile(ANTinterface_path, 'gelDuke-128_noEOG_refmatrix.mat'), 'ref_matrix')

%% set up data for re-referencing if not provided as first input
if isempty(data)
    % the referencing matrices in gelDuke-128_noEOG_refmatrix.mat assume a
    % particular ordering of channels. Let's pull the data into the
    % expected ordering.
    
    assert(isfield(EEG, 'chanlocs'), 'Cannot proceed with referencing as EEG is missing chanlocs.')
    
    load('ANT_montage_templates.mat', 'chanlocs_dukeZ3')
    template_labels = {chanlocs_dukeZ3.labels};
    
    % drop the EOG channel
    template_labels = template_labels(1:128);
    
    % build the channel indices
    channelindex = nan(size(template_labels));
    labels = {EEG.chanlocs.labels};
    
    for ii = 1:length(channelindex)
        idx = find(cellfun(@(x) strcmp(x, template_labels{ii}), labels));
        if ~isempty(idx)
            channelindex(ii) = idx;
        end
    end
    
    assert(~any(isnan(channelindex)), 'Some channels required for referencing cannot be found.')
    
    % Pull the data into the expected ordering of the referencing matrices
    data = EEG.data(channelindex, :);
end

%% Referencing
if isa(ref_method, 'double')
    assert(length(ref_method) == 1, 'More than one channel provided for unipolar re-referencing.')
    
    if strcmp(EEG.refscheme, EEG.chanlocs(ref_method).labels)
        disp(['EEG data already under the specified referencing scheme: ', EEG.refscheme])
        disp('ANT_interface_reref() existed without referencing.')
    else
        EEG.refscheme = EEG.chanlocs(ref_method).labels;
        
        if verbose
            disp(' ')
            disp(['Re-referencing EEG data: unipolar ', num2str(ref_method)])
            disp(' ')
        end
        
        R2 = eye(size(data,1));
        R2(:, ref_method) = R2(:, ref_method) - 1;
        data = R2 * data;
    end
else
    assert(isa(ref_method,'char'), 'Check ref_method argument, not double or string character')
    if strcmp(EEG.refscheme, ref_method)
        disp(['EEG data already under the specified referencing scheme: ', EEG.refscheme])
        disp('ANT_interface_reref() existed without referencing.')
    else
        EEG.refscheme = ref_method;
        
        switch ref_method
            case 'AR'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: common average')
                    disp(' ')
                end
                data = squeeze(ref_matrix(1,:,:)) * data;
                
            case 'Z3'
                if verbose
                    warning('Are you sure you want to reference to Z3? Raw data should be recorded under Z3 reference already. Check EEG.ref_scheme!')
                    disp(' ')
                    disp('Re-referencing EEG data: recording reference at Z3')
                    disp(' ')
                end
                data = squeeze(ref_matrix(2,:,:)) * data;
                
            case 'left mastoid'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: left mastoid')
                    disp(' ')
                end
                data = squeeze(ref_matrix(3,:,:)) * data;
                
            case 'linked mastoid'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: linked mastoid')
                    disp(' ')
                end
                data = squeeze(ref_matrix(4,:,:)) * data;

            case 'REST'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: REST')
                    disp(' ')
                end
                if ~isempty(G)
                    assert(size(G,1) == size(data,1), 'Number of channels in G doesn"t match that of data.')
                    fREST = pinv(G)'*pinv(G)*ones(size(G,1), 1) ./ (ones(size(G,1), 1)'*pinv(G)'*pinv(G)*ones(size(G,1), 1));
                    R5 = eye(size(G,1)) - ones(size(G,1), 1)*fREST';
                    data = R5 * data;
                else
                    disp('Lead field matrix G not provided, will use the LFM from example data for REST referencing.')
                    data = squeeze(ref_matrix(5,:,:)) * data;
                    EEG.refscheme = [ref_method, '_exampleG'];
                end
            
            case 'contral mastoid'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: linked mastoid')
                    disp(' ')
                end
                data = squeeze(ref_matrix(6,:,:)) * data;
                
            case 'LP'
                if verbose
                    disp(' ')
                    disp('Re-referencing EEG data: manual laplacian')
                    disp(' ')
                end
                data = squeeze(ref_matrix(7,:,:)) * data;
                
        end
    end
end

%% update the EEG structure with re-referenced data
if justEEG
    EEG.data(channelindex,:) = data;
elseif size(EEG.data, 1) == size(data, 1)
    EEG.data = data;
else
    disp('Data inputed have different numbers of channel from EEG.data. EEG struct will not be updated.')
    assert(nargout > 1, 'Only 1 output specified, re-referenced data cannot be accessed outside this function!')
end

end
