% allTablesNew.m
% Cross-platform runner with 10-second timeout per script.
%
% Mac/Linux:
%   - shows live output
%   - saves log files
%   - kills only the child MATLAB process on timeout
%
% Windows:
%   - saves output to log files
%   - kills only the child MATLAB process on timeout
%   - does not require Parallel Computing Toolbox

clear; clc;

rootDir = fileparts(mfilename('fullpath'));

%%%Set Time Limit for Each call to Table File
maxTime = 10; % seconds

if ispc
    matlabExec = fullfile(matlabroot, 'bin', 'matlab.exe');
else
    matlabExec = fullfile(matlabroot, 'bin', 'matlab');
end

if ~isfile(matlabExec)
    error('MATLAB executable not found:\n%s', matlabExec);
end

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
    fullfile(rootDir, 'Tables_9-10_Paper_PCG', 'table_10_PCG_Linux.m')
};

fprintf('Root folder:\n%s\n', rootDir);
fprintf('Timeout per script: %d seconds\n', maxTime);
fprintf('MATLAB executable:\n%s\n\n', matlabExec);

for k = 1:numel(scripts)
    thisScript = scripts{k};

    if ~isfile(thisScript)
        warning('File not found. Skipping:\n%s', thisScript);
        continue;
    end

    scriptDir = fileparts(thisScript);
    [~, scriptName, ~] = fileparts(thisScript);

    logFile = fullfile(rootDir, sprintf('log_%02d_%s.txt', k, scriptName));

    fprintf('\n============================================================\n');
    fprintf('Running %d/%d with %d-second timeout:\n%s\n', ...
        k, numel(scripts), maxTime, thisScript);
    fprintf('Log file:\n%s\n', logFile);
    fprintf('============================================================\n\n');

    matlabCmd = sprintf("cd('%s'); run('%s');", ...
        escapeForMatlab(scriptDir), ...
        escapeForMatlab(thisScript));

    if ispc
        % Windows: safe PID-based timeout using PowerShell.
        % Output is saved to logFile, but not streamed live.
        psCmd = sprintf([ ...
            '$p = Start-Process -FilePath ''%s'' ' ...
            '-ArgumentList ''-batch "%s"'' ' ...
            '-RedirectStandardOutput ''%s'' ' ...
            '-RedirectStandardError ''%s.err'' ' ...
            '-PassThru; ' ...
            'if (-not $p.WaitForExit(%d000)) { ' ...
            'Add-Content -Path ''%s'' -Value "`nTIMEOUT: killed after %d seconds."; ' ...
            'Stop-Process -Id $p.Id -Force; ' ...
            '} else { ' ...
            'Add-Content -Path ''%s'' -Value "`nFinished normally."; ' ...
            '}'], ...
            escapeForPowerShell(matlabExec), ...
            escapeForPowerShell(matlabCmd), ...
            escapeForPowerShell(logFile), ...
            escapeForPowerShell(logFile), ...
            maxTime, ...
            escapeForPowerShell(logFile), maxTime, ...
            escapeForPowerShell(logFile));

        shellCmd = sprintf('powershell -NoProfile -Command "%s"', psCmd);

    else
        % macOS/Linux: live output via tee, plus safe PID-based timeout.
        pidFile = fullfile(tempdir, sprintf('matlab_child_pid_%d.txt', k));

        shellCmd = sprintf([ ...
            '( "%s" -batch "%s" 2>&1 & echo $! > "%s" ) | tee "%s" & ' ...
            'tee_pid=$!; ' ...
            'sleep %d; ' ...
            'matlab_pid=$(cat "%s" 2>/dev/null); ' ...
            'if [ -n "$matlab_pid" ] && kill -0 $matlab_pid 2>/dev/null; then ' ...
            'echo "\nTIMEOUT: killed after %d seconds." | tee -a "%s"; ' ...
            'kill -9 $matlab_pid 2>/dev/null; ' ...
            'fi; ' ...
            'wait $tee_pid 2>/dev/null; ' ...
            'rm -f "%s"'], ...
            escapeForShell(matlabExec), ...
            escapeForShell(matlabCmd), ...
            escapeForShell(pidFile), ...
            escapeForShell(logFile), ...
            maxTime, ...
            escapeForShell(pidFile), ...
            maxTime, ...
            escapeForShell(logFile), ...
            escapeForShell(pidFile));
    end

    status = system(shellCmd);

    if status ~= 0
        warning('System command returned nonzero status for:\n%s', thisScript);
    end
end

fprintf('\nDone running all scripts.\n');

function s = escapeForMatlab(s)
    s = strrep(s, '''', '''''');
end

function s = escapeForShell(s)
    s = strrep(s, '"', '\"');
end

function s = escapeForPowerShell(s)
    s = strrep(s, '''', '''''');
end

