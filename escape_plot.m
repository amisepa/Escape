%   Copyright (C) Cedric Cannard 2025 – Escape EEGLAB Plugin (https://github.com/amisepa/Escape)


function escape_plot(entropyData, chanlocs, entropyType, scales)

% chanlabels = {chanlocs.labels};
% load("colormap_rufin.mat");
% load("colormap_bwr.mat");
% load("colormap_bgy.mat");

% uniscale or multiscale data
if size(entropyData,2)>1
    multiscale = true;
else
    multiscale = false;
end

% % deal with near-0 values
% idx = entropyData < 0.0001;
% if sum(idx) > 0
%     entropyData(idx) = 0.0001;
%     warning(['Channel ' chanlocs(idx).labels ' is probably a bad channel.'])
% end
if ~multiscale
    entropyData(entropyData==0) = NaN;
    % else
    % entropyData(entropyData==0,:) = NaN;
end
% mymap(1,:) = [.9 .9 .9]; % set NaNs to gray

% % old topo plot
% if ~multiscale
%     plot2 = false;
% end

% mass univariate plot scales (x-axis) x channels (y-axis)
if multiscale

    % main plot
    figure('Color','w','InvertHardCopy','off');
    subplot(3, 3, [1 2 4 5 7 8])
    hold on
    nScales = size(entropyData,2);
    nChan = size(entropyData,1);

    % Main plot
    imagesc(1:nScales, 1:nChan, entropyData);
    
    % y-ticks: EEG channels
    Yticks = {chanlocs.labels};
    if nChan > 30
        newticks = 1:2:length(Yticks);
    else
        newticks = 1:length(Yticks);
    end
    newticks = unique(newticks);
    Yticks  = Yticks(newticks);
    set(gca,'YTick',newticks,'YTickLabel', Yticks,'FontWeight','normal');

    % x-ticks: Time scales
    Xticks = scales;
    if nScales > 30
        newticks = 1:2:length(Xticks);
    else
        newticks = 1:length(Xticks);
    end
    newticks = unique(newticks);
    Xticks  = Xticks(newticks);
    set(gca,'XTick',newticks,'XTickLabel', Xticks,'FontWeight','normal');
    

    % Colormap
    colormap("parula"); 
    c = colorbar; c.Label.String = 'Entropy';

    % labels
    xlabel('Scales'); ylabel('EEG channels');
    ylabel(c, 'Entropy','FontWeight','bold','FontSize',9)
    title(sprintf('%s',entropyType))

    % cleanup
    axis tight
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');

    % Clim
    % maxval = max(abs(entropyData(:)));
    % if max(effects(:)) < 0
    %     clim([-maxval 0])
    % elseif min(effects(:)) > 0
    %     clim([0 maxval])
    % else
    % clim([-maxval maxval])
    % clim([-min(entropyData,[],'all') max(entropyData,[],'all')])


    % Peak entropy, channel, and scale
    % peakVal = max(entropyData, [], 'all');
    [peak_value, linear_idx] = max(entropyData(:));  % Find the max value and its linear index
    [peak_channel, peak_scale] = ind2sub(size(entropyData), linear_idx);
    % fprintf('Peak entropy value is %.2f at Scale %g and Channel %s. \n', peak_value, scales(peak_scale), chanlocs(peak_channel).labels);
    fprintf('Peak entropy value is %.2f at Scale %g (%.2f-%2.f Hz) for Channel %s. \n', peak_value, peak_scale, scales{peak_scale}(1), scales{peak_scale}(2), chanlocs(peak_channel).labels);


    % time series of peak cluster
    if ~isempty(chanlocs)
        subplot(3,3,6)
        hold on
        plot(1:length(scales), entropyData(peak_channel,:),'LineWidth',2);
        % color1 = [0, 0.4470, 0.7410];
        axis tight; box on

        % chanLabel = chanlocs(peakChan).labels;
        % title(sprintf('Course plot: %s',chanLabel),'FontSize',11,'fontweight','bold')
        % % plot(xaxis,stats(cluster_maxe,:),'LineWidth',2);  % plot peak effect of all clusters superimposed
        % % chanLabel = {chanlocs(cluster_maxe).labels};
        % % legend(chanLabel)
        % grid on; axis tight;
        % ylabel('t-values','FontSize',11,'fontweight','bold');
        % xlabel('Frequency (Hz)','FontSize',11,'fontweight','bold')
        %
        % % Plot bars of significnace for peak electrode
        % plotSigBar(mask(peakChan,:)~=0,xaxis);

        % Topography of peak scale
        subplot(3,3,9)
        topoplot(entropyData(:,peak_scale), chanlocs, 'emarker', {'.','k',8,1},'electrodes','on');
        clim([min(entropyData(:,peak_scale)) max(entropyData(:,peak_scale))]);
        colormap('parula');  % 'hot' 'parula'
        % c.Label.FontSize = 11;
        % c.Label.FontWeight = 'bold';

        title(entropyType);

        set(gcf,'Name','Multiscale entropy visualization','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
        set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    end

else
    % Topography of uniscale entropy
    figure('Color','w','InvertHardCopy','off');
    topoplot(entropyData, chanlocs, 'emarker', {'.','k',15,1},'electrodes','labels');
    colormap('parula'); % dmap mymap diverging_bgy 'hot' 'bone' 'winter' 'summer' 'viridis'

    % clim and colorbar
    % clim([min(entropyData)*.9 max(entropyData)*1.1]);
    clim([min(entropyData)*.95 max(entropyData)*1.05]);
    % clim([min(entropyData) max(entropyData)]);
    % clim([0 max(entropyData)]);
    c = colorbar;
    c.Label.String = 'Entropy';
    c.Label.FontSize = 11;
    c.Label.FontWeight = 'bold';

    % Title
    title(entropyType);

    set(gcf,'Name','Uniscale entropy visualization','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');

end

% catch
%     plot2 = true; % use back up plot


% if plot2
%     disp('Not enough electrodes or variance across electrodes to plot the scalp topography. Defaulting to secondary plot. ')
%
%     x = [ chanlocs.X ]';
%     y = [ chanlocs.Y ]';
%     z = [ chanlocs.Z ]';
%
%     % Rotate X Y Z coordinates
%     % rotate = 0;       %nosedir = +x
%     rotate = 3*pi/2;    %nosedir = +y
%     % rotate = pi;      %nosedir = -x
%     % rotate = pi/2;
%     allcoords = (y + x.*sqrt(-1)).*exp(sqrt(-1).*rotate);
%     x = imag(allcoords);
%     y = real(allcoords);
%
%     % Project 3D positions on 2D plane if not already done
%     chanpos(:,1) = x;
%     chanpos(:,2) = y;
%     chanpos(:,3) = z;
%
%     if all(chanpos(:,3)==0)
%         coord = chanpos(:,1:2); % positions already projected on a 2D plane
%     else
%         coord = chanpos; % use 3-D data for plotting
%     end
%
%     % 3D figure allowing to open entropy data for each electrode
%     p = figure('color','w');
%     % p.Position = [100 100 540 400];
%     axis equal
%     axis vis3d
%     axis off
%     hold on
%
%     % adj = mean(entropyData(1,:))*5; % to scale marker size
%
%     for iChan = 1:size(entropyData,1)
%
%         if length(entropyData(iChan,:)) == 1 % measures with one value per channel
%             % 3D plot of entropy values at electrode locations
%             p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
%                 'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
%                 'Marker','o','MarkerSize',5);
%
%             % Display channel label + entropy value for each channel
%             text(coord(iChan,1)-15,coord(iChan,2)+10,coord(iChan,3), ...
%                sprintf('%s: %6.1f',chanlabels{iChan}, ...
%                round(entropyData(iChan,:),2)),'FontSize',10,'fontweight','bold');
%
%         else % for multiscales, take area under the curve as sensor size
%             p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
%                 'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
%                 'Marker','o','MarkerSize', 5, 'UserData',iChan, ...
%                 'ButtonDownFcn', @(~,~,~) buttonCallback(entropyData(iChan,:), coord(iChan,:), chanlabels{iChan}));
%
%             % Display channel label above each electrode
%             text(coord(iChan,1)-7,coord(iChan,2)+10,coord(iChan,3), ...
%                 sprintf('%s %6.3f',chanlabels{iChan}), ...
%                 'FontSize',10,'fontweight','bold');
%             title('[Click on sensors to display entropy values]', ...
%                 'Position', [1 120 1], 'fontweight', 'bold')
%
%         end
%     end
%
%     % set(gcf,'Name','Multiscale entropy visualization','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
%     % set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
%
% end



% % subfunction to display entropy values in the plot where use clicks
% function buttonCallback(tmpdata, coor, label)
% 
% % Entropy measures with only one value per channel
% figure('color','w','Position', [500 500 280 210]);
% plot(tmpdata,'linewidth',2,'color','black'); % blue: [0, 0.4470, 0.7410]
% % area(tmpdata,'linewidth',2);
% title(label,'FontSize',14)
% % xticks(2:nScales); xticklabels(join(string(scales(:,2:end)),1)); xtickangle(45)
% % xlim([2 nScales]);
% xlabel('Time scale','FontSize',12,'fontweight','bold');
% ylabel('Entropy','FontSize',12,'fontweight','bold')

