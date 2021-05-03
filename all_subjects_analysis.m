clear

%% loop through all subject foldes
dir_name = 'E:\Min and Steph\subjects_data\2018';
D = dir(dir_name);
for k = 3:length(D)
    % Get the current subdirectory name
    currD = D(k).name;
    currD_name = [dir_name, '\', currD];
    disp(currD_name);
    
    % change to finall data folder
    % and loop through 4 conditions
    sub_D_name = [currD_name, '\data_for_LD'];
    cd(sub_D_name);
    sub_D = dir(sub_D_name);
    
    % begin loop through each condition and location
    for l = 3 : length(sub_D)
        condition_D = sub_D(l).name;
        condition_D_name = [sub_D_name, '\', condition_D];
        cd(condition_D_name);
        
        % condition code and location
        condition_location = split(condition_D, '_');
        condition = condition_location{1};
        location = condition_location{2};
        disp(condition)
        disp(location)
        
        % load velocity, diameter, labchart
        velocity_name = 'velocity.mat';
        diameter_name = 'diameter.csv';
        labchart_name = 'clean.mat';
        destination_folder_name = 'E:\Min and Steph\subjects_data\final_data_2018';
        
        % begin aligment scripts
        align_labchart_velocity(velocity_name, diameter_name, labchart_name, destination_folder_name);
        
        % begin average scripts
        
        
        
        
        
        
        
        
        % change back to origianl directory
        cd(sub_D_name);
    end
    
    
    
    % change back to orignal directory
    cd(dir_name);
end
