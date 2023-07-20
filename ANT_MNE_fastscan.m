function [ fastscan ] = ANT_MNE_fastscan(pcfn, markfn, fsfn, filepath, criteria, verbose)
%
% ANT MNE CODES - FASTSCAN
%
% - function used to prepare the digitization of electrode locations from
% Polhemus FastScanII scannner into data arrays and electrode labels to be
% read by mne.channels.read_dig_montage to create the montage for
% assembling a _raw.fif file for input in mne.gui.coregistration when
% generating the -trans.fif file during forward modeling in MNE python.
%
% Last edit: Alex He 12/07/2022
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - pcfn:         filename of the .mat file exported from
%                           Polhemus FastScanII software containing the
%                           point clouds of the head surface as "Cloud of 
%                           Points" (not triangulation face indices).
%
%           - markfn:       filename of the .txt file exported from
%                           Polhemus FastScanII software (in a separate
%                           step from exporting the point cloud) containing
%                           the coordinates of the markers manually placed
%                           by RAs according to the instruction manual. It
%                           must follow a certain order when placing the
%                           markers in the Polhemus FastScanII software in
%                           order to be read correctly. 
%
%           - fsfn:         filename of the .mat file created by this
%                           function containing the various electrode and
%                           fiducial landmark coordinates both in native
%                           FastScan acquisition coordinate and in the
%                           newly transformed Duke template coordinate
%                           system, which has the origin (0,0,0) in the
%                           center of the head. This filename typically
%                           should have the same naming as pcfn and markfn,
%                           but with a _dig suffix appended behind to
%                           indicate it contains information about
%                           digitization of electrodes. However, you can
%                           give it arbitrary filename in this function.
%                           
%                           This same fsfn filename will also be used to
%                           name the .csv file exported from this function
%                           containing the digitization electrode
%                           coordinates and fiducial landmark coordinates
%                           read by functions in ANT_MNE_python_util.py
%                           to create the digitized montage object in MNE
%                           python.
%
%           - filepath:     full path to the folder containing the FastScan
%                           exported files (.mat and .txt) as well as the
%                           location to save the output quality-check
%                           figures produced from this function and the
%                           final fastscan structure (name specified by
%                           fsfn) to be saved to a .mat file.
%
%           - criteria:     a 1x2 vector specifying the tolerable ranges of
%                           angle deviation and distance deviation in
%                           transforming from FastScan native coordinate
%                           system to the Duke coordinate system. 
%                           defaul: [2, 0.2]
%
%           - verbose:      whether print messages and plotting during
%                           processing subject's electrode location files
%                           default: true
%
% Output:
%           - fastscan:     a structure containing all information of the
%                           digitization of electrodes acquired by Polhemus
%                           FastScanII scanner. The fields of this
%                           structure are as following:
%
%                           fastscan.head 
%                               - point cloud of head surface
%                           fastscan.electrode 
%                               - coordinates of electrodes in native
%                               FastScan coordinate system
%                           fastscan.landmark 
%                               - coordinates of right preauricular point,
%                               nasion, and left preauricular point in
%                               native FastScan coordinate system
%                           fastscan.electrode_dukexyz 
%                               - coordinates of electrodes in the new Duke
%                               Waveguard configuration coordinate system
%                           fastscan.landmark_dukexyz 
%                               - coordinates of right preauricular
%                               point, nasion, and left preauricular point
%                               in the new Duke Waveguard configuration
%                               coordinate system
%                           fastscan.elc_labels 
%                               - names of the channels in both
%                               fastscan.electrode and
%                               fastscan.electrode_dukexyz
%                           fastscan.landmark_labels 
%                               - names of the fiducial landmarks in both
%                               fastscan.landmark and
%                               fastscan.landmark_dukexyz
%                           fastscan.chanlocs_fs
%                               - a structure containing the electrode
%                               coordinates in the new Duke Waveguard
%                               configuration coordinate system and the
%                               same electrode order as EEG.data to be
%                               integrated with an EEGLAB structure (EEG)
%                               for the field chanlocs
%                           fastscan.chanlocs_duke 
%                               - a structure containing the electrode
%                               coordinates of the Duke Waveguard template
%                               contained in the field chanlocs in an
%                               EEGLAB structure (EEG) obtained from
%                               calling ANT_interface_readcnt.m function 
%                           fastscan.chanlocs_duke_reord
%                               - a structure containing the electrode
%                               coordinates of the Duke Waveguard template
%                               contained in the field chanlocs in an
%                               EEGLAB structure (EEG) obtained from
%                               calling ANT_interface_readcnt.m function,
%                               re-ordered to the same ordering as the
%                               order of electrodes marked in Polhemus
%                               FastScanII software (the order of 
%                               fastscan.electrode,
%                               fastscan.electrode_dukexyz, and
%                               fastscan.elc_labels)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin < 5
    criteria = [2,0.2];
    verbose = true;
elseif nargin < 6
    verbose = true;
end

angle_tol = criteria(1);
distance_tol = criteria(2);

% addpath to the appropriate folders
try 
    SleepEEG_addpath(matlabroot);
    
    ANTinterface_path = which('ANT_MNE_fastscan');
    temp = strsplit(ANTinterface_path, 'ANT_MNE_fastscan.m');
    ANTinterface_path = temp{1};
    
catch
    % if using SleepEEG_addpath() fails, we will assume the current directory
    % has the ANT_MNE_fastscan.m or at least the folder containing it has
    % been added to path when calling this function. We will try to addpath to
    % EEGLAB directly.
    
    ANTinterface_path = which('ANT_MNE_fastscan');
    temp = strsplit(ANTinterface_path, 'ANT_MNE_fastscan.m');
    ANTinterface_path = temp{1};
    
    % Add path to EEGLAB
    addpath(fullfile(ANTinterface_path, 'eeglab14_1_2b'))
end

% Add paths to dependent EEGLAB functions
eeglabpath = which('eeglab');
temp = strsplit(eeglabpath, 'eeglab.m');
addpath(fullfile(temp{1}, 'functions', 'sigprocfunc'))
addpath(fullfile(temp{1}, 'functions', 'guifunc'))
addpath(fullfile(temp{1}, 'functions', 'adminfunc'))

% Set warning of graphical display errors in command window to off
warning('off', 'MATLAB:callback:error')

% Use default naming of fsfn
if isempty(fsfn)
    temp = strsplit(pcfn, '.mat');
    fsfn = [temp{1},'_dig.mat'];
end

% Extract a general fileID for saved figures
temp = strsplit(fsfn, '.mat');
fileID = temp{1};

%% Visualize the data point clouds and confirm anatomical landmarks
fsfn = [fileID '.mat'];
fsfn_full = fullfile(filepath, fsfn);

if exist(fsfn_full, 'file') == 2 % if already on disk, load it
    if verbose
        disp(' ')
        disp([fsfn, ' is already on disk. Loading from:'])
        disp(' ')
        disp(fsfn_full)
    end
    load(fsfn_full, 'fastscan')
else
    if verbose
        disp(' ')
        disp([fsfn, ' is not created yet.'])
        disp(' ')
        disp('We will now create it from the .mat and .txt files from Polhemus FastScanII...')
    end
    
    % Load the exported files from Polhemus FastScanII
    pcloud = load(fullfile(filepath, pcfn));
    
    fileID_f = fopen(fullfile(filepath, markfn),'r');
    startRow = 4;
    formatSpec = '%10f%10f%f%[^\n\r]';
    dataArray = textscan(fileID_f, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID_f);
    markers = [dataArray{1:end-1}];
    clearvars filename startRow formatSpec fileID_f dataArray ans;    
%     markers = readtable(fullfile(filepath, markfn));
%     markers = table2array(markers);
    
    % Show the point cloud of head surfaces
    figure
    set(gcf, 'units', 'pixels', 'Position', [0 0 1400 1000]);
    pcshow(pcloud.Points');
    hold on
    pcshow([0,0,0], 'g', 'MarkerSize', 6000)
    grid off; axis off
    
    % Display the landmarks
    pcshow(markers, 'w', 'MarkerSize', 5000);
    
    % ------------------
    % Manual checkpoint!
    % ------------------
    % Check if the first three correspond to anatomical landmarks
    landmarks = markers(1:3, :);
    ldindex = 1:3;
    pcshow(landmarks, 'r', 'MarkerSize', 5000);
    title('Head Surface, Electrodes, and Anatomical Landmarks', 'FontSize', 16)
    endsearch = input('Are the anatomical landmarks correctly colored red? (y/n): ', 's');
    if strcmpi(endsearch, 'n')
        % try the last three as anatomical landmarks
        hold on
        pcshow(landmarks, 'k', 'MarkerSize', 5000);
        landmarks = markers(end-2:end, :);
        ldindex = size(markers,1)-2:size(markers,1);
        hold on
        pcshow(landmarks, 'r', 'MarkerSize', 5500);
        endsearch = input('What about now, are the anatomical landmarks correctly colored red? (y/n): ', 's');
        if strcmpi(endsearch, 'n')
            % If still not, it means the electrode ordering is messed up!
            save(fullfile(filepath, [fileID '_landmark_order_debug_workspace']))
            error('Manual check failed. Workspace is saved. Please debug...')
            %             hold on
            %             pcshow(landmarks, 'k', 'MarkerSize', 5000);
            %             % Accept inputs for three landmark XYZ values
            %             XYZvalue = input('Please enter the XYZ value of first landmark separated by comma (X,Y,Z): ', 's');
            %             loc(1,:) = cellfun(@str2double, strsplit(XYZvalue, [" ",',']));
            %             XYZvalue = input('Please enter the XYZ value of second landmark separated by comma (X,Y,Z): ', 's');
            %             loc(2,:) = cellfun(@str2double, strsplit(XYZvalue, [" ",',']));
            %             XYZvalue = input('Please enter the XYZ value of third landmark separated by comma (X,Y,Z): ', 's');
            %             loc(3,:) = cellfun(@str2double, strsplit(XYZvalue, [" ",',']));
            %
            %             % Find these landmark points
            %             ldindex = [];
            %             landmarks = zeros(size(loc));
            %             for i = 1:size(markers,1)
            %                 for j = 1:size(loc,1)
            %                     if sum(abs(markers(i,:) - loc(j,:))) < 1
            %                         landmarks(j,:) = landmarks(i,:);
            %                         deleteindex = [deleteindex, i]; %#ok<AGROW>
            %                     end
            %                 end
            %             end
            %             hold on
            %             pcshow(landmarks, 'r', 'MarkerSize', 5500);
        end
    end
    close all
    
    %% Create a structure to store these XYZ values
    % remove the landmarks from markers, the remaining ones should be electrodes
    markers(ldindex, :) = [];
    fastscan = struct;
    fastscan.head = pcloud.Points;
    fastscan.electrode = markers;
    fastscan.landmark = landmarks;
    
    %% Visualize just the electrodes and landmarks and save the plot
    %
    %     figure
    %     pcshow(fastscan.electrode, 'w', 'MarkerSize', 2000);
    %     hold on
    %     pcshow(fastscan.landmark, 'r', 'MarkerSize', 4000);
    %     grid off; axis off
    %     title('Electrodes and Landmarks', 'FontSize', 16)
    %
    %     savefig(fullfile(dataDir, subID, 'fastscan', [fileID '_visual_elc.fig']))
    %     close all
    %
    %% Prepare a montage file for mne.gui.coregistration
    
    % Need to convert the XYZ to a new coordinate system with different origin
    % and direction vectors based on the Duke template in EEGLAB
    
    % We define the new Duke coordinate system the following way:
    % - We find the centroid on the left side among LD3, LD4, and LC3, call
    % it Lmid; we find the centroid on the right side among RD3, RD4, and RC3
    % call it Rmid. We connect Lmid and Rmid, take the midpoint of this line as
    % the origin.
    % - Y direction is pointing in Lmid. Then we connect Z direction
    % as pointing towards the vertex electrode Z5. And X is normal to the plane
    % and pointing towards the front, such that Z1-Z4 should have positive X
    % values.
    
    % FASTSCAN digitization in the old coordinate system (native tracker based)
    oldxyz = fastscan.electrode;
    
    % Electrodes in the Duke template coordinates (target coordinate system)
    load(fullfile(ANTinterface_path, 'duke_129_template.mat'), 'chanlocs')
    duke = zeros(length(chanlocs), 3);
    for i = 1:length(chanlocs)
        duke(i,1) = chanlocs(i).X;
        duke(i,2) = chanlocs(i).Y;
        duke(i,3) = chanlocs(i).Z;
    end
    
    % Let's first visualize the involved electrodes in the Duke template. We
    % won't try to convert the Duke template coordinate system to this new
    % coordinate system as well since we are trying to emulate it. But it gives
    % us a sense of how good the approximation is in the coordiante system that
    % we are trying to mimic from.
    if verbose
        figure
        set(gcf, 'units', 'pixels', 'Position', [200 200 1600 600]);
        subplot(1,2,1)
        pcshow(duke, 'k', 'MarkerSize', 2000);
        xlabel('X'); ylabel('Y'), zlabel('Z')
        hold on
        pcshow(duke(124,:), 'b', 'MarkerSize', 2000); % LD3
        pcshow(duke(125,:), 'b', 'MarkerSize', 2000); % LD4
        pcshow(duke(7,:), 'b', 'MarkerSize', 2000); % LC3
        pcshow(duke(50,:), 'b', 'MarkerSize', 2000); % RD3
        pcshow(duke(51,:), 'b', 'MarkerSize', 2000); % RD4
        pcshow(duke(40,:), 'b', 'MarkerSize', 2000); % RC3
        pcshow(duke(86,:), 'b', 'MarkerSize', 2000); % Z5
        % Origin
        pcshow([0,0,0], 'm', 'MarkerSize', 2000);
        % Direction vectors
        plot3([0,100], [0,0], [0,0], 'r', 'LineWidth', 3)
        plot3([0,0], [0,100], [0,0], 'g', 'LineWidth', 3)
        plot3([0,0], [0,0], [0,100], 'b', 'LineWidth', 3)
        title('Duke Template', 'FontSize', 16)
        set(gca, 'Color', 'w')
        
        % Let's visualize these key electrodes in the FASTSCAN digitalization
        subplot(1,2,2)
        pcshow(oldxyz, 'k', 'MarkerSize', 2000);
        xlabel('X'); ylabel('Y'), zlabel('Z')
        hold on
        pcshow(fastscan.landmark, 'r', 'MarkerSize', 4000);
        pcshow(oldxyz(119,:), 'b', 'MarkerSize', 2000); % LD3
        pcshow(oldxyz(120,:), 'b', 'MarkerSize', 2000); % LD4
        pcshow(oldxyz(112,:), 'b', 'MarkerSize', 2000); % LC3
        pcshow(oldxyz(8,:), 'b', 'MarkerSize', 2000); % RD3
        pcshow(oldxyz(9,:), 'b', 'MarkerSize', 2000); % RD4
        pcshow(oldxyz(15,:), 'b', 'MarkerSize', 2000); % RC3
        pcshow(oldxyz(62,:), 'c', 'MarkerSize', 2000); % Z5
        % Origin
        pcshow([0,0,0], 'm', 'MarkerSize', 2000);
        % Direction vectors
        plot3([0,100], [0,0], [0,0], 'r', 'LineWidth', 3)
        plot3([0,0], [0,100], [0,0], 'g', 'LineWidth', 3)
        plot3([0,0], [0,0], [0,100], 'b', 'LineWidth', 3)
        title('FastScan Digitization', 'FontSize', 16)
        set(gca, 'Color', 'w')
    end
    
    % Define centroids of LD3, LD4, LC3, and RD3, RD4, RC3
    Lmid = [mean([oldxyz(119,1), oldxyz(120,1), oldxyz(112,1)]),...
        mean([oldxyz(119,2), oldxyz(120,2), oldxyz(112,2)]),...
        mean([oldxyz(119,3), oldxyz(120,3), oldxyz(112,3)])];
    Rmid = [mean([oldxyz(8,1), oldxyz(9,1), oldxyz(15,1)]),...
        mean([oldxyz(8,2), oldxyz(9,2), oldxyz(15,2)]),...
        mean([oldxyz(8,3), oldxyz(9,3), oldxyz(15,3)])];
    
    % Find the new origin point
    neworigin = (Lmid + Rmid)./2;
    
    % Display these centroid points
    if verbose
        hold on
        pcshow([Lmid;Rmid], 'c', 'MarkerSize', 2000);
        % Plot the directions of the new coordinate system
        plot3([Lmid(1), neworigin(1)], [Lmid(2), neworigin(2)], [Lmid(3), neworigin(3)], 'g', 'LineWidth', 3) % Y direction
        pcshow(neworigin, 'm', 'MarkerSize', 3000);
        plot3([neworigin(1), oldxyz(62,1)], [neworigin(2), oldxyz(62,2)], [neworigin(3), oldxyz(62,3)], 'b', 'LineWidth', 3) % Z direction
        newnosep = cross(Lmid-neworigin, oldxyz(62,:)-neworigin)./100 + neworigin;
        plot3([newnosep(1), neworigin(1)], [newnosep(2), neworigin(2)], [newnosep(3), neworigin(3)], 'r', 'LineWidth', 3) % X direction
        set(gca, 'Color', 'w')
    end
    
    % Compute the new unit vectors XYZ
    newuniX = (newnosep-neworigin)./norm(newnosep-neworigin);
    newuniY = (Lmid-neworigin)./norm(Lmid-neworigin);
    newuniZ = (oldxyz(62,:)-neworigin)./norm(oldxyz(62,:)-neworigin);
    
    % Confirm that the new unit vectors are roughly orthogonal
    xy_off = abs(90 - atan2d(norm(cross(newuniX,newuniY)),dot(newuniX,newuniY)));
    yz_off = abs(90 - atan2d(norm(cross(newuniY,newuniZ)),dot(newuniY,newuniZ)));
    zx_off = abs(90 - atan2d(norm(cross(newuniZ,newuniX)),dot(newuniZ,newuniX)));
    if ~all([xy_off, yz_off, zx_off] < angle_tol) % off angle should be less than [angle_tol]degrees (default = 2deg)
        disp(['xy angle off: ', num2str(xy_off), 'deg'])
        disp(['yz angle off: ', num2str(yz_off), 'deg'])
        disp(['zx angle off: ', num2str(zx_off), 'deg'])
        error(['new coordinate direction vectors are not orthogonal (enough). off angle > ', num2str(angle_tol), 'deg. This is usually due to erroneous definition of Lmid or Rmid.'])
    end
    
    % Now transform electrode coordinates into the new Duke coordinate system
    newxyz = (oldxyz-neworigin) * [newuniX', newuniY', newuniZ'];
    
    % transform the landmarks as well
    newlandmark = (fastscan.landmark-neworigin) * [newuniX', newuniY', newuniZ'] ;
    
    % Visualize the electrode locations in the new Duke coordinate system
    if verbose
        figure
        set(gcf, 'units', 'pixels', 'Position', [200 200 1600 600]);
        subplot(1,2,1)
        pcshow(oldxyz, 'k', 'MarkerSize', 2000);
        xlabel('X'); ylabel('Y'), zlabel('Z')
        hold on
        pcshow(fastscan.landmark, 'r', 'MarkerSize', 4000);
        pcshow(oldxyz(119,:), 'b', 'MarkerSize', 2000); % LD3
        pcshow(oldxyz(120,:), 'b', 'MarkerSize', 2000); % LD4
        pcshow(oldxyz(112,:), 'b', 'MarkerSize', 2000); % LC3
        pcshow(oldxyz(8,:), 'b', 'MarkerSize', 2000); % RD3
        pcshow(oldxyz(9,:), 'b', 'MarkerSize', 2000); % RD4
        pcshow(oldxyz(15,:), 'b', 'MarkerSize', 2000); % RC3
        pcshow(oldxyz(62,:), 'c', 'MarkerSize', 2000); % Z5
        % Origin
        pcshow([0,0,0], 'm', 'MarkerSize', 2000);
        % Direction vectors
        plot3([0,100], [0,0], [0,0], 'r', 'LineWidth', 3)
        plot3([0,0], [0,100], [0,0], 'g', 'LineWidth', 3)
        plot3([0,0], [0,0], [0,100], 'b', 'LineWidth', 3)
        title('FastScan Native Coordinate System', 'FontSize', 16)
        set(gca, 'Color', 'w')
        
        subplot(1,2,2)
        pcshow(newxyz, 'k', 'MarkerSize', 2000);
        xlabel('X'); ylabel('Y'), zlabel('Z')
        hold on
        pcshow(newlandmark, 'r', 'MarkerSize', 4000);
        pcshow(newxyz(119,:), 'b', 'MarkerSize', 2000); % LD3
        pcshow(newxyz(120,:), 'b', 'MarkerSize', 2000); % LD4
        pcshow(newxyz(112,:), 'b', 'MarkerSize', 2000); % LC3
        pcshow(newxyz(8,:), 'b', 'MarkerSize', 2000); % RD3
        pcshow(newxyz(9,:), 'b', 'MarkerSize', 2000); % RD4
        pcshow(newxyz(15,:), 'b', 'MarkerSize', 2000); % RC3
        pcshow(newxyz(62,:), 'c', 'MarkerSize', 2000); % Z5
        % Origin
        pcshow([0,0,0], 'm', 'MarkerSize', 2000);
        % Direction vectors
        plot3([0,100], [0,0], [0,0], 'r', 'LineWidth', 3)
        plot3([0,0], [0,100], [0,0], 'g', 'LineWidth', 3)
        plot3([0,0], [0,0], [0,100], 'b', 'LineWidth', 3)
        title('New Duke Coordinate System', 'FontSize', 16)
        set(gca, 'Color', 'w')
    end
    % Make sure that distance between preauricular points are not signicantly
    % altered during this coordiante system transformation
    preauc_distance = norm(newlandmark(1,:)-newlandmark(3,:)) - norm(fastscan.landmark(1,:)-fastscan.landmark(3,:));
    if ~(preauc_distance < distance_tol) % distance change should be less than [distance_tol]mm (default = 0.2mm)
        disp(['Preauricular point distance changed by: ', num2str(preauc_distance), 'mm'])
        error(['Distance between PA points significantly changed during coordiante transformation. Distance off > ', num2str(distance_tol), 'mm.'])
    end
    
    %% Now, let's do some quality control inspections
    % The next two steps have to be done manually...
    figure
    set(gcf, 'units', 'pixels', 'Position', [200 200 800 600]);
    scatter3(newxyz(:,1), newxyz(:,2), newxyz(:,3), 400, 'k', 'filled');
    axis equal
    rotate3d on
    xlabel('X'); ylabel('Y'), zlabel('Z')
    hold on
    scatter3(newlandmark(1,1), newlandmark(1,2), newlandmark(1,3), 600, 'r', 'filled');
    scatter3(newlandmark(2,1), newlandmark(2,2), newlandmark(2,3), 600, 'g', 'filled');
    scatter3(newlandmark(3,1), newlandmark(3,2), newlandmark(3,3), 600, 'b', 'filled');
    % Direction vectors
    plot3([0,100], [0,0], [0,0], 'r', 'LineWidth', 3)
    plot3([0,0], [0,100], [0,0], 'g', 'LineWidth', 3)
    plot3([0,0], [0,0], [0,100], 'b', 'LineWidth', 3)
    legend('Electrodes', 'Right PA', 'Nasion', 'Left PA', 'X', 'Y', 'Z')
    title('Landmarks in Duke Coordinates')
    set(gca, 'FontSize', 20)
    view(270, 90)
    
    % ------------------
    % Manual checkpoint!
    % ------------------
    checkok =  input('Do the anatomical landmarks look ok? (y/n): ', 's');
    if strcmpi(checkok, 'n')
        save(fullfile(filepath, [fileID '_landmark_debug_workspace']))
        error('Manual check failed. Workspace is saved. Please debug...')
    else
        % Save the plot as .png for reference
        figure(3)
        view(270, 90) % reset to top-down view point
        saveas(gcf, fullfile(filepath, [fileID '_landmark_check.png']))
        close all
    end
    
    %% Create a chanlocs structure for the subject
    % loads in the labels of fastscan electrodes done by manual labelling
    % of the RAs
    load(fullfile(ANTinterface_path, 'duke_129_template.mat'), 'fs_label_order') 
    for i = 1:length(fs_label_order)
        fs_label_order(i).X = newxyz(i,1);
        fs_label_order(i).Y = newxyz(i,2);
        fs_label_order(i).Z = newxyz(i,3);
    end
    
    % Convert to EEGLAB 2D polar coordinates for topoplot
    fschanlocs = convertlocs(fs_label_order, 'cart2all');
    
    % Load reordered Duke template chanlocs such that it has the same order
    % as the desired FastScan labeling of electrodes (as in fs_template)
    load(fullfile(ANTinterface_path, 'duke_129_template.mat'), 'chanlocs_reord')
    
    % Make EEGLAB 2D topoplots of the electrode locations in the FastScan
    % order for
    % 1) Duke template electrodes
    % 2) Digitization that was acquired for each subject
    % We want to manually compare the numbers are approximately identify at
    % different electrode positions to make sure no labeling error was made
    % when using stylus pen to mark the electrode positions in FastScanII
    % software
    
    figure;
    set(gcf, 'Position', [200 200 1600 800])
    
    % Duke template
    subplot(1,2,1)
    topoplot([],chanlocs_reord,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});
    L = findobj(gcf, 'type', 'Text');
    for ind = 1:length(chanlocs_reord)
        set(L(length(chanlocs_reord)+1-ind), 'FontSize', 14)
    end
    title([fileID ' Template Positions'], 'FontSize', 30, 'Interpreter', 'none')
    
    % FastScan Digitization
    subplot(1,2,2)
    topoplot([],fschanlocs,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});
    L = findobj(gcf, 'type', 'Text');
    for ind = 1:length(fschanlocs)
        set(L(length(fschanlocs)+1-ind), 'FontSize', 14)
        set(L(length(fschanlocs)+1-ind), 'Color', [0,0,1])
    end
    title([fileID ' Digitization Positions'], 'FontSize', 30, 'Color', [0,0,1], 'Interpreter', 'none')
    
    % ------------------
    % Manual checkpoint!
    % ------------------
    checkok =  input('Does the digitization have identical numbering order as the template? (y/n): ', 's');
    if strcmpi(checkok, 'n')
        save(fullfile(filepath, [fileID '_elcorder_debug_workspace']))
        error('Manual check failed. Workspace is saved. Please debug...')
    else
        % Save the plot as .png for reference
        saveas(gcf, fullfile(filepath, [fileID '_elcorder_check.png']))
        close all
    end
    
    %% Now that both landmark orders and electrode orders are vetted, store them
    % Digitization in Duke coordinates
    fastscan.electrode_dukexyz = newxyz;
    fastscan.landmark_dukexyz = newlandmark;
    
    % Digitization labels
    elclabels = cell(length(fschanlocs), 1);
    for i = 1:length(fschanlocs)
        elclabels{i} = fschanlocs(i).labels;
    end
    fastscan.elc_labels = elclabels;
    fastscan.landmark_labels = {'rpa'; 'nasion'; 'lpa'}; % Right Preauricular Point, Nasion, Left Preauricular Point
    
    % Re-order into the Duke template electrode order that is also the
    % order of collected EEG.data
    chanlocs_fs = fschanlocs;
    for i = 1:length(fschanlocs)
        chanlocs_fs(fschanlocs(i).ognum) = fschanlocs(i);
    end
    fastscan.chanlocs_fs = chanlocs_fs;
    
    % Duke template values
    fastscan.chanlocs_duke = chanlocs;
    fastscan.chanlocs_duke_reord = chanlocs_reord;
    
    % Save the structure so we don't have to do these manual checkings again
    if verbose
        disp(' ')
        disp('FastScan digitization manual check completed, fastscan file saved to: ')
        disp(' ')
        disp(fsfn_full)
    end
    save(fsfn_full, 'fastscan')
    
end

%% Store in formats readable by Python
% We have confirmed the correctness of electrode labelling by
% transforming the digitization locations into EEGLAB 2D topoplot of
% the Duke coordinate system of the template. We are confident in the
% electrode and landmark labels.

% Rather than using the transformed locations, we can save the raw
% digitzation electrode + landmark locations in its native coordinate
% system based on the tracker. mne.channels.read_dig_montage will
% construct its own transformation based on the head coordinate system
% in MNE-Python, which is defined by the 3 anatomical landmarks.

% Now we just have to save these raw XYZ values and electrode +
% landmark labels to formats readable into creating Python arrays. We
% can then pass these arrays to mne.channels.read_dig_montage to create
% the montage object for mne.gui.coregistration, which is trying to
% create the _trans.fif file for forward modeling solution. We will save
% as .csv file here. 

% First let's store the labels and XYZ values into a table. 
fs_dig = table;
fs_dig.labels = char([fastscan.landmark_labels; fastscan.elc_labels]);
fs_dig.X = [fastscan.landmark(:,1); fastscan.electrode(:,1)];
fs_dig.Y = [fastscan.landmark(:,2); fastscan.electrode(:,2)];
fs_dig.Z = [fastscan.landmark(:,3); fastscan.electrode(:,3)];

% Write to .csv file

csvfn = fullfile(filepath, [fileID '.csv']);
writetable(fs_dig, csvfn, 'Delimiter', ',')

if verbose
    disp(' ')
    disp('Labels and XYZ values for sleepeeg_create_montage is saved to .csv file: ')
    disp(' ')
    disp(csvfn)
end

%%
% turn the warning setting back on
warning('on', 'MATLAB:callback:error')

end
