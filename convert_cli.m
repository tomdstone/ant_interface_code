function convert_cli(old_file_name, new_file_name)
%     disp(old_file_name)
   
    old_path_parts = strsplit(old_file_name, '/');
    old_path = strjoin(old_path_parts(1:end-1), '/');
    old_file = old_path_parts{end};
%     disp(old_file);

    new_path_parts = strsplit(new_file_name, '/');
    new_path = strjoin(new_path_parts(1:end-1), '/');
    new_file = new_path_parts{end};
 
    eeg = ANT_interface_readcnt(old_file, old_path);

    [~] = ANT_interface_saveset(eeg, new_file, new_path);
end