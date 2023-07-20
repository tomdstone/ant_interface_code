function [ EEG ] = ANT_interface_eego2asa(EEG, capconfig)
%
% ANT INTERFACE CODES - EEGO2ASA
%
% - function used to reorder some of the channels to convert the data
% collected on eego mylab EEG system to the same numbering order as the
% asalab and on the Duke Waveguard cap (old version).

% % Update 10/03/2019: when using the new ANT Duke Waveguard Caps with Z3
% reference, we need a different set of tranformation of data, as we now
% have an extra dropdown eye-electrode. We use this same function to
% configure EEG struct to the correct channel order and give back the
% reference channel as a zero-line recording channel.
%
% Last edit: Alex He 04/19/2023
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - EEG:          an EEG structure with EEG.data in the order of
%                           eego mylab EEG system.
%
% Output:
%           - EEG:          an EEG structure with EEG.data in the order of
%                           asalab / Duke Waveguard cap configuration.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Change the order of data and impedance measures

disp('eego mylab data loaded. ANT_interface_eego2asa() is called to make channel order consistent with Duke configuration.')

switch capconfig
    case 'old-asa' % old asa cap collected on eego mylab system with adaptors
        asa_order = [1:73,...   % up to Channel 73 the same
            95,...     % R1(#74) is Channel 95 on eego
            74:93,...  % Channels 75-94 are shifted upwards by 1
            96,...     % L4(#95) is Channel 96 on eego
            94,...     % L5(#96) is Channel 94 on eego
            97:128];   % Channels 97-128 are the same
        
        % Change the data channel order
        EEG.data = EEG.data(asa_order, :);
        
        % Change the impedance channel order
        if ~isempty(EEG.initimp)
            EEG.initimp = EEG.initimp(asa_order);
        end
        
        if ~isempty(EEG.endimp)
            EEG.endimp = EEG.endimp(asa_order);
        end
        
    case 'new-Z3' % new ANT Duke cap collected on eego mylab system, with dropdown EOG sensor
        % the only change now is channel 85 is recorded to be EOG, and the
        % original channel 85 (Z3) electrode is not explicitly recorded
        % since it is the reference electrode. We will add it back in for
        % the ease of fastscan / forward modeling.
        
        tempdata = EEG.data;
        
        % store away extra channels when there are any
        if EEG.nbchan >= 129
            tempdata(130:EEG.nbchan+1,:) = tempdata(129:end,:);
        end
        
        tempdata(129,:) = tempdata(85,:); % now number 129 is the EOG
        tempdata(85,:) = zeros(1, size(tempdata,2)); % reference channel as zero recording
        EEG.data = tempdata;
                
        % Add one more channel in the channel location
        if EEG.nbchan >= 129
            EEG.chanlocs(130:EEG.nbchan+1) = EEG.chanlocs(129:end);
        end
        EEG.chanlocs(129) = EEG.chanlocs(85);
        EEG.chanlocs(85).labels = 'Z3';
        
        % Change the impedance channel order
        if ~isempty(EEG.initimp)
            tempimp = EEG.initimp;
            if EEG.nbchan >= 129
                tempimp(130:EEG.nbchan+1) = tempimp(129:end);
            end
            tempimp(129) = tempimp(85);
            tempimp(85) = NaN;
            EEG.initimp = tempimp;
        end
        
        if ~isempty(EEG.endimp)
            tempimp = EEG.endimp;
            if EEG.nbchan >= 129
                tempimp(130:EEG.nbchan+1) = tempimp(129:end);
            end
            tempimp(129) = tempimp(85);
            tempimp(85) = NaN;
            EEG.endimp = tempimp;
        end
        
        % Update number of channels
        EEG.nbchan = size(EEG.data,1);
        
    case 'asa-EEG' % old asa cap collected on old asalab EEG system
        % no action needed
end

end
