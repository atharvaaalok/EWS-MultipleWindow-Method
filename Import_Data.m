function Data = Import_Data(DataFolder_path, System_name)
    
    % Call the import function for the respective system
    switch System_name
        case 'PowerSystem'
            Data = Import_PowerSystem_Data(DataFolder_path);
        case 'RijkeTube'
            Data = Import_RijkeTube_Data(DataFolder_path);
        otherwise
            warning('Unexpected system name. Check name again.');
            exit();
    end


end