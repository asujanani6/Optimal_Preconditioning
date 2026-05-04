function all_tables()

originalPath = path;

scripts = {
    'Tables_1-5_Paper_One_Sided/table_1_Suite.m'
    'Tables_1-5_Paper_One_Sided/table_2_Rand.m'
    'Tables_1-5_Paper_One_Sided/table_3_Suite.m'
    'Tables_1-5_Paper_One_Sided/table_4_Rand.m'
    'Tables_1-5_Paper_One_Sided/table_5_PCG.m'
    'Tables_6-8_Paper_Two_Sided/table_6_Suite.m'
    'Tables_6-8_Paper_Two_Sided/table_7_Rand.m'
    'Tables_6-8_Paper_Two_Sided/table_8_Suite.m'
    'Tables_9-10_Paper_PCG/table_9_PCG.m'
    'Tables_9-10_Paper_PCG/table_10_PCG_Linux.m'
};

for k = 1:numel(scripts)
    fprintf('\nRunning %s\n', scripts{k});

    try
        evalin('base', sprintf('run(''%s'')', scripts{k}));
    catch ME
        path(originalPath);
        rethrow(ME);
    end

    path(originalPath);
end

end
