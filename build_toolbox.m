% build_toolbox.m
% Builds the USTB toolbox (.mltbx) for release.
% Version is injected from the TOOLBOX_VERSION environment variable
% when run via GitHub Actions, or defaults to 'dev' for local builds.

version = getenv('TOOLBOX_VERSION');
if isempty(version)
    version = 'dev';
    fprintf('No TOOLBOX_VERSION set, using "%s"\n', version);
else
    fprintf('Building USTB v%s\n', version);
end

output_file = sprintf('USTB_v%s.mltbx', version);

% Load options from the existing .prj file (preserves all your
% file lists, icons, dependencies etc.), then override version + output
opts = matlab.addons.toolbox.ToolboxOptions('USTB.prj');
opts.ToolboxVersion = version;
opts.OutputFile     = output_file;

matlab.addons.toolbox.packageToolbox(opts);

fprintf('Toolbox built: %s\n', output_file);
