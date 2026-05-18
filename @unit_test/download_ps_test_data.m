function download_ps_test_data()
%DOWNLOAD_PS_TEST_DATA  Ensure ps/ test .mat files are available locally.
%   Downloads ps.zip from Zenodo and unzips into data/ps/ if not already present.

ps_dir = fullfile(ustb_path(), 'data', 'ps');
if isfolder(ps_dir) && ~isempty(dir(fullfile(ps_dir, '*.mat')))
    return
end

url = 'https://zenodo.org/records/20269473/files';
data_dir = fullfile(ustb_path(), 'data');
zipfile = fullfile(data_dir, 'ps.zip');

tools.download('ps.zip', url, data_dir);
unzip(zipfile, data_dir);
if isfile(zipfile)
    delete(zipfile);
end
end
