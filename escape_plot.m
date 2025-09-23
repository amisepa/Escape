function escape_plot(entropyData, chanlocs, entropyType, scales, varargin)
% ESCAPE_PLOT Visualize entropy maps (uni-/multi-scale and time-resolved)
% • Uniscale (chan×1): scalp topography.
% • Multiscale (chan×scale): heatmap (channels×scales), line at peak channel, topo at peak scale.
% • Time-resolved (chan×scale×time): adds a time-course at the peak (channel,scale).
%
% escape_plot(entropyData, chanlocs, entropyType, scales)
% escape_plot(entropyData, chanlocs, entropyType, scales, time_sec)
%
% INPUTS
% entropyData : [nChan×1] or [nChan×nScales] or [nChan×nScales×nTime]
% chanlocs : EEGLAB channel locations (with .labels and positions)
% entropyType : title string (e.g., 'MSE (m=2, r=0.15)')
% scales : cellstr of scale labels (e.g., {'1','2',...} or '[lo hi] Hz'); numeric also accepted
% time_sec : (optional) [1×nTime] window centers in seconds for time-resolved data
%
% BEHAVIOR
% • Heatmap uses time-average if a 3-D array is supplied.
% • Peak (channel,scale) is taken from the heatmap; printed to console.
% • If time_sec is provided (and sized nTime), the time axis uses seconds; otherwise window index.
%
% OUTPUTS
% Creates a figure (no return value): heatmap, peak-channel curve, peak-scale topography,
% and (if 3-D input) a time-course at the peak.
%
% EXAMPLE
% escape_plot(EEG.escape.MSE.data, EEG.chanlocs, 'MFE', EEG.escape.MSE.scales)
%
% -------------------------------------------------------------------------
% Copyright (C) 2025 EEGLAB Escape plugin — Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% Optional time vector for time-resolved MSE
time_sec = [];
if ~isempty(varargin)
    time_sec = varargin{1};
end

% Detect uniscale / multiscale / time-resolved
isTimeResolved = ndims(entropyData) == 3;      % [chan x scale x time]
if isTimeResolved
    entropyData3D = entropyData;
    % For the main heatmap, use the time-average (keeps your existing layout)
    entropyData   = mean(entropyData3D, 3, 'omitnan');
end
multiscale = size(entropyData,2) > 1;

% Handle exact zeros (can be artifact of upstream steps)
if ~multiscale
    entropyData(entropyData==0) = NaN;
end

if multiscale
    % ===== Multiscale heatmap + per-channel curve + topo (and optional time plot) =====
    figure('Color','w','InvertHardCopy','off');
    % Main heatmap occupies left 2 columns (all rows)
    subplot(3,3,[1 2 4 5 7 8]); hold on;

    nScales = size(entropyData,2);
    nChan   = size(entropyData,1);

    imagesc(1:nScales, 1:nChan, entropyData); axis tight;
    set(gca,'TickDir','out'); box on;

    % Y ticks (channel labels)
    Yticks = {chanlocs.labels};
    if nChan > 30
        newY = round(linspace(1, nChan, min(20, nChan)));  % a bit more readable than every other
    else
        newY = 1:nChan;
    end
    set(gca,'YTick',newY,'YTickLabel',Yticks(newY),'FontWeight','normal');

    % X ticks (scales)
    if iscell(scales)
        Xticks = scales; nX = numel(scales);
    else
        Xticks = arrayfun(@(x){num2str(x)}, 1:nScales); nX = nScales;
    end
    if nX > 30, newX = round(linspace(1, nX, min(20, nX))); else, newX = 1:nX; end
    set(gca,'XTick',newX,'XTickLabel',Xticks(newX),'FontWeight','normal');
    
    % ===== FIX CROPPING: rotate labels + grow bottom margin =====
    ax = gca;
    ax.TickLabelInterpreter = 'none';          % prevents underscores from becoming subscripts
    ax.XTickLabelRotation  = 45;               % 45° usually enough; use 60 if your labels are very long
    ax.PositionConstraint  = 'outerposition';  % allow us to control outer box
    
    % Lift the axes a bit and shrink height so rotated labels fit
    % drawnow;                                    % ensure TightInset is up to date
    outer = ax.OuterPosition;                   % [left bottom width height] in normalized
    ti    = ax.TightInset;                      % [l b r t] padding needed for tick labels etc.
    extra = 0.02;                               % a touch more bottom space (tweak if needed)
    newBottom = max(outer(2), ti(2) + extra);
    newHeight = max(outer(4) - (newBottom - outer(2)) - (ti(4) + 0.01), 0.1);
    ax.OuterPosition = [outer(1), newBottom, outer(3), newHeight];
    
    % make the colorbar, labels, title as you had -----
    colormap('parula');
    c = colorbar; ylabel(c,'Entropy','FontWeight','bold','FontSize',9);
    xlabel('Scales'); ylabel('EEG channels');
    title(entropyType, 'Interpreter','none');


    % % X ticks (scales)
    % if iscell(scales)
    %     Xticks = scales;
    %     nX = numel(scales);
    % else
    %     Xticks = arrayfun(@(x){num2str(x)}, 1:nScales);
    %     nX = nScales;
    % end
    % if nX > 30
    %     newX = round(linspace(1, nX, min(20, nX)));
    % else
    %     newX = 1:nX;
    % end
    % set(gca,'XTick',newX,'XTickLabel',Xticks(newX),'FontWeight','normal');
    % 
    % % Colormap & colorbar
    % colormap('parula');
    % c = colorbar; ylabel(c,'Entropy','FontWeight','bold','FontSize',9);

    % % Labels & title
    % xlabel('Scales'); ylabel('EEG channels');
    % title(entropyType, 'Interpreter','none');

    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','bold');

    % ---- Peak (from time-averaged map) ----
    [peak_value, linear_idx] = max(entropyData(:));
    [peak_channel, peak_scale] = ind2sub(size(entropyData), linear_idx);

    % Robust scale label for printing
    if iscell(scales)
        sclabel = scales{peak_scale};
    else
        sclabel = num2str(peak_scale);
    end
    fprintf('Peak entropy: %.3f at Scale %s (index %d), Channel %s.\n', ...
        peak_value, sclabel, peak_scale, chanlocs(peak_channel).labels);

    % ---- Per-channel curve at peak channel ----
    subplot(3,3,6); hold on; box on;
    plot(1:nScales, entropyData(peak_channel,:), 'LineWidth', 2);
    xlim([1 nScales]); xlabel('Scale'); ylabel('Entropy');
    title(sprintf('Channel %s', chanlocs(peak_channel).labels), 'Interpreter','none');

    % ---- Topography at peak scale ----
    subplot(3,3,3);
    topoplot(entropyData(:,peak_scale), chanlocs, 'emarker',{'.','k',8,1}, 'electrodes','on');
    clim([min(entropyData(:,peak_scale)) max(entropyData(:,peak_scale))]);
    colormap('parula');
    title(sprintf('%s @ scale %s', entropyType, sclabel), 'Interpreter','none');

    % ---- Time-resolved trace at the peak channel & scale ----
    if isTimeResolved
        subplot(3,3,9); hold on; box on;
        tr = squeeze(entropyData3D(peak_channel, peak_scale, :));
        if isempty(time_sec) || numel(time_sec) ~= numel(tr)
            t = 1:numel(tr);
            xlabel('Window #');
        else
            t = time_sec(:);
            xlabel('Time (s)');
        end
        plot(t, tr, 'LineWidth', 1.5);
        ylabel('Entropy');
        title(sprintf('Time course @ %s, scale %s', chanlocs(peak_channel).labels, sclabel), 'Interpreter','none');
    end

    set(gcf,'Name','Multiscale entropy visualization','Color','w','Toolbar','none','Menu','none','NumberTitle','Off');

else
    % ===== Uniscale topography =====
    figure('Color','w','InvertHardCopy','off');
    topoplot(entropyData, chanlocs, 'emarker', {'.','k',15,1}, 'electrodes','labels');
    colormap('parula');
    clim([min(entropyData)*0.95, max(entropyData)*1.05]);
    c = colorbar; c.Label.String = 'Entropy'; c.Label.FontSize = 11; c.Label.FontWeight = 'bold';
    title(entropyType, 'Interpreter','none');
    set(gcf,'Name','Uniscale entropy visualization','Color','w','Toolbar','none','Menu','none','NumberTitle','Off');
    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','bold');
end
