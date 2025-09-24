% eegplugin_escape()
%
% ESCAPE: Entropy and Signal Complexity Analysis Plugin for EEG/MEG
%
% Alternative name ASCENT: Aperiodic Signal Complexity Estimation for Neurophysiological Time series
%
% Copyright (C) - EEGLAB, Cedric Cannard, August 2022-2025
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_escape(fig,try_strings,catch_strings)

% Plugin version
vers = '1.2';

% Add paths to subfolders
p = fileparts(which('eegplugin_escape.m'));
addpath(p);
addpath(fullfile(p,'functions'))

cmd = [ try_strings.check_data '[EEG,LASTCOM] = escape_compute(EEG);' ...
        catch_strings.new_and_hist ];

% create menu
toolsmenu = findobj(fig, 'tag', 'tools');
uimenu(toolsmenu, 'label', 'Compute entropy', 'userdata', 'startup:off;epoch:off;study:off', ...
    'callback', cmd, 'position', 15);

end
