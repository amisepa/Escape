function ascent_plot(entropyData, chanlocs, entropyType, scales, varargin)
% ascent_plot  Visualize uni/multi/time-resolved entropy with robust NaN handling.
%
% entropyData : [chan x scale] or [chan x scale x time]
% chanlocs    : EEGLAB channel locations
% entropyType : label, e.g., 'RCMFEσ' or 'MSE (std)'
% scales      : cellstr of scale labels or numeric 1:S
% varargin{1} : optional time vector (seconds) for time-resolved data

time_sec = [];
if ~isempty(varargin), time_sec = varargin{1}; end

isTimeResolved = ndims(entropyData) == 3;      % [chan x scale x time]
if isTimeResolved
    entropyData3D = entropyData;
    entropyData   = mean(entropyData3D, 3, 'omitnan'); % for the main heatmap
end
multiscale = size(entropyData,2) > 1;

if ~multiscale
    entropyData(entropyData==0) = NaN;  % handle 0-artifacts
end

if multiscale
    % ===== Multiscale heatmap + per-channel curve + topo (+ optional time) =====
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
        newY = round(linspace(1, nChan, min(20, nChan)));
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

    % --- Fix label cropping: rotate + expand bottom margin
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    ax.XTickLabelRotation  = 45;            % or 60 if labels are very long
    ax.PositionConstraint  = 'outerposition';
    drawnow;
    outer = ax.OuterPosition;
    ti    = ax.TightInset;
    extra = 0.02;
    newBottom = max(outer(2), ti(2) + extra);
    newHeight = max(outer(4) - (newBottom - outer(2)) - (ti(4) + 0.01), 0.1);
    ax.OuterPosition = [outer(1), newBottom, outer(3), newHeight];

    % Colormap & colorbar
    colormap('parula');
    c = colorbar; ylabel(c,'Entropy','FontWeight','bold','FontSize',9);

    % Labels & title
    xlabel('Scales'); ylabel('EEG channels');
    title(entropyType, 'Interpreter','none');

    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','bold');

    % ---- Peak (robust to NaNs; fallback to best scale) ----
    finiteMask = isfinite(entropyData);
    if ~any(finiteMask(:))
        warning('All values are NaN; nothing to plot.');
        return
    end
    tmp = entropyData; tmp(~finiteMask) = -Inf;
    [peak_value, linear_idx] = max(tmp(:));
    [peak_channel, peak_scale] = ind2sub(size(tmp), linear_idx);
    if ~any(isfinite(entropyData(:,peak_scale)))
        [~, peak_scale] = max(sum(isfinite(entropyData),1));
    end
    if iscell(scales), sclabel = scales{peak_scale}; else, sclabel = num2str(peak_scale); end
    fprintf('Peak entropy: %.3f at Scale %s (index %d), Channel %s.\n', ...
        peak_value, sclabel, peak_scale, chanlocs(peak_channel).labels);

    % ---- Per-channel curve at peak channel ----
    subplot(3,3,6); hold on; box on;
    row = entropyData(peak_channel,:);
    if all(~isfinite(row)), row = nan(1,nScales); end
    plot(1:nScales, row, 'LineWidth', 2);
    xlim([1 nScales]); xlabel('Scale'); ylabel('Entropy');
    title(sprintf('Channel %s', chanlocs(peak_channel).labels), 'Interpreter','none');

    % ---- Topography at peak scale ----
    subplot(3,3,3);
    vals = entropyData(:,peak_scale);
    finiteVals = vals(isfinite(vals));
    if isempty(finiteVals)
        axis off
        text(0.5, 0.5, sprintf('No finite values @ scale %s', sclabel), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    else
        try
            topoplot(vals, chanlocs, 'emarker',{'.','k',8,1}, 'electrodes','on');
            lo = min(finiteVals); hi = max(finiteVals);
            if isfinite(lo) && isfinite(hi) && lo < hi, clim([lo hi]); end
            colormap('parula');
            title(sprintf('%s @ scale %s', entropyType, sclabel), 'Interpreter','none');
        catch 
            warning('topoplot failed: %s. Falling back to bar chart.');
            bar(vals); xlim([0 numel(vals)+1]); box on;
            title(sprintf('Bar topography @ scale %s', sclabel), 'Interpreter','none');
            ylabel('Entropy');
        end
    end

    % ---- Time-resolved trace at the peak channel & scale ----
    if isTimeResolved
        subplot(3,3,9); hold on; box on;
        tr = squeeze(entropyData3D(peak_channel, peak_scale, :));
        if isempty(tr) || all(~isfinite(tr))
            axis off
            text(0.5,0.5,'No finite time-resolved values at the selected channel/scale',...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
        else
            if isempty(time_sec) || numel(time_sec) ~= numel(tr)
                t = 1:numel(tr); xlabel('Window #');
            else
                t = time_sec(:); xlabel('Time (s)');
            end
            plot(t, tr, 'LineWidth', 1.5);
            ylabel('Entropy');
            title(sprintf('Time course @ %s, scale %s', chanlocs(peak_channel).labels, sclabel), 'Interpreter','none');
        end
    end

    set(gcf,'Name','Multiscale entropy visualization','Color','w','Toolbar','none','Menu','none','NumberTitle','Off');

else
    % ===== Uniscale topography =====
    figure('Color','w','InvertHardCopy','off');
    vals = entropyData(:);
    finiteVals = vals(isfinite(vals));
    if isempty(finiteVals)
        axis off
        text(0.5,0.5,'No finite values to plot','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    else
        topoplot(vals, chanlocs, 'emarker', {'.','k',15,1}, 'electrodes','labels');
        colormap('parula');
        clim([min(finiteVals)*0.95, max(finiteVals)*1.05]);
        c = colorbar; c.Label.String = 'Entropy'; c.Label.FontSize = 11; c.Label.FontWeight = 'bold';
        title(entropyType, 'Interpreter','none');
    end
    set(gcf,'Name','Uniscale entropy visualization','Color','w','Toolbar','none','Menu','none','NumberTitle','Off');
    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','bold');
end
end
