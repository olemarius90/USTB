function correction_download_uff_with_fallbacks(filename, dst_dir)
%CORRECTION_DOWNLOAD_UFF_WITH_FALLBACKS  Download UFF: Zenodo and/or ustb.no mirrors.
%
%   Several Vrålstad (T-USON) speckle chamber .uff files are hosted on ustb.no only
%   (not in Zenodo record 19550715). For those names we try mirrors first to avoid 404.

outfile = fullfile(dst_dir, filename);
if isfile(outfile)
    return
end

mirror_bases = {
    'https://www.ustb.no/datasets'
    'http://www.ustb.no/datasets'
    'https://ustb.no/datasets'
    'http://ustb.no/datasets'
    'https://www.ultrasoundtoolbox.com/datasets'
    };
zenodo_base = tools.zenodo_dataset_files_base();

% Hosts prior to Zenodo for datasets known to be missing from 19550715
pfx = 'speckle_sim_FI_P4_';
if numel(filename) >= numel(pfx) && strcmp(filename(1:numel(pfx)), pfx)
    bases = [mirror_bases; {zenodo_base}];
else
    bases = [{zenodo_base}; mirror_bases];
end
lastME = [];

for k = 1:numel(bases)
    base = strtrim(char(bases{k}));
    if isempty(base)
        continue
    end
    try
        tools.download(filename, base, dst_dir);
        if isfile(outfile)
            return
        end
    catch ME
        lastME = ME;
    end
end

if ~isfile(outfile)
    if ~isempty(lastME)
        rethrow(lastME)
    end
    error('correction_download_uff_with_fallbacks:missing', ...
        'Could not download ''%s'' from Zenodo or ustb mirrors.', filename);
end

end
