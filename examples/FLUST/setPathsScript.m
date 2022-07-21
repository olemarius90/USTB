lastwarn(''); % reset warning

%% set paths here

% PSF generators (Field II/MUST/other)
addpath('PathToField_II')
addpath('PathToMUST');

% FLUST folders
addpath( genpath( pwd) ); % please run FLUST from its root folder
addpath('..\..'); % ustb main folder

%% check if addpath returned warnings
[warnmsg, msgid] = lastwarn;
if strcmp( msgid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
    disp('At least one addpath statement in setPathsScript returned with a warning.')
    edit setPathsScript.m
    error('Check that paths to PSF simulator (Field/MUST/other) are set properly in setPathsScript.m');
end