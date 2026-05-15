function correction_publish_beam_snap(b_obj, headless_tf, varargin)
%CORRECTION_PUBLISH_BEAM_SNAP  beamformed_data.plot tuned for MATLAB publish() in -batch.
% Passing axes on throwaway invisible figures omitted snapshots for HTML embed; passing an
% explicit figure handle (matching the Save PNGs section) restores output. Stripping ui
% controls avoids alternateGetframe "Printing of uicontrols is not supported".

if headless_tf
    fh_snap = figure('Visible', 'off');
    if isempty(varargin)
        b_obj.plot(fh_snap);
    else
        b_obj.plot(fh_snap, varargin{:});
    end
else
    if isempty(varargin)
        b_obj.plot([]);
    else
        b_obj.plot([], varargin{:});
    end
end

drawnow;
if nargin >= 2 && headless_tf && exist('fh_snap', 'var') == 1 && isgraphics(fh_snap)
    delete(findall(fh_snap, 'Type', 'uicontrol'));
end
snapnow;
if nargin >= 2 && headless_tf && exist('fh_snap', 'var') == 1 && isgraphics(fh_snap)
    close(fh_snap);
end

end
