function correction_download_uff_with_fallbacks(filename, dst_dir)
%CORRECTION_DOWNLOAD_UFF_WITH_FALLBACKS  Download UFF: Zenodo (CI), then ustb.no mirrors.

outfile = fullfile(dst_dir, filename);
if isfile(outfile)
    return
end

bases = {
    tools.zenodo_dataset_files_base()
    'https://www.ustb.no/datasets'
    'http://www.ustb.no/datasets'
    'https://ustb.no/datasets'
    'http://ustb.no/datasets'
    'https://www.ultrasoundtoolbox.com/datasets'
    };
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
