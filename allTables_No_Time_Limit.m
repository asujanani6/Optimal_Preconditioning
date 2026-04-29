% allTablesNew.m No Time Limit Working
% Run all table scripts from the codesPCGandNCG root folder.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));

scripts = {
    fullfile(rootDir, 'Tables_1-5_Paper_One_Sided', 'New_Tables_1_2_Draft2', 'table_1_Suite.m')
    fullfile(rootDir, 'Tables_1-5_Paper_One_Sided', 'New_Tables_1_2_Draft2', 'table_2_Rand.m')
    fullfile(rootDir, 'Tables_1-5_Paper_One_Sided', 'table_3_Suite.m')
    fullfile(rootDir, 'Tables_1-5_Paper_One_Sided', 'table_4_Rand.m')
    fullfile(rootDir, 'Tables_1-5_Paper_One_Sided', 'table_5_PCG.m')
    fullfile(rootDir, 'Tables_6-8_Paper_Two_Sided', 'table_6_Suite.m')
    fullfile(rootDir, 'Tables_6-8_Paper_Two_Sided', 'table_7_Rand.m')
    fullfile(rootDir, 'Tables_6-8_Paper_Two_Sided', 'table_8_Suite.m')
    fullfile(rootDir, 'Tables_9-10_Paper_PCG', 'table_9_PCG.m')
    fullfile(rootDir, 'Tables_9-10_Paper_PCG', 'table_PCG_Linux.m')
};

fprintf('Running all table scripts from:\n%s\n\n', rootDir);

for k = 1:numel(scripts)
    thisScript = scripts{k};

    if ~isfile(thisScript)
        error('File not found:\n%s', thisScript);
    end

    fprintf('\n============================================================\n');
    fprintf('Running %d/%d:\n%s\n', k, numel(scripts), thisScript);
    fprintf('============================================================\n');

    scriptDir = fileparts(thisScript);
    oldDir = pwd;

    try
        cd(scriptDir);
        run(thisScript);
        cd(oldDir);
    catch ME
        cd(oldDir);
        fprintf('\nERROR while running:\n%s\n', thisScript);
        rethrow(ME);
    end
end

fprintf('\nAll table scripts finished successfully.\n');
