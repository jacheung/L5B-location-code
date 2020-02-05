%% BUILD
hilbert_feature = {'angle','phase','midpoint','amplitude','velocity'};

for g = 1:2
    wStruct= whisk_location_quant_vsupp(U,1:length(U),hilbert_feature{g},'off');
    cd('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\tmp_structs')
        save(['whisk_' hilbert_feature{g} '_window_v2'],'wStruct')
end
%% OR LOAD
hilbert_feature = 'angle';
mdl_versions = {'instant','window'};

for b = 1:numel(mdl_versions)
    fileName = ['whisk_' hilbert_feature '_' mdl_versions{b}];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
        whisk_struct.(mdl_versions{b}) = wStruct;
    end
end

itune = find(cellfun(@(x) x.is_tuned==1,whisk_struct.instant));
wtune = find(cellfun(@(x) x.is_tuned==1,whisk_struct.window));
both_tuned = intersect(itune,wtune);

%% PLOT RANDOM TRIAL
% cellNum = datasample(intersect(itune,wtune),1);
[r] = numSubplots(numel(both_tuned));
figure(5480);clf
for g =1:numel(both_tuned)
    cellNum = both_tuned(g);
    
    iraw_x = cellfun(@median,(whisk_struct.instant{cellNum}.stim_response.raw_stim));
    iraw_y = cellfun(@mean,(whisk_struct.instant{cellNum}.stim_response.raw_response));
    
    ix = whisk_struct.instant{cellNum}.stim_response.bars_stim;
    iy = whisk_struct.instant{cellNum}.stim_response.bars_fit.mean;
    
    raw_x = cellfun(@median,(whisk_struct.window{cellNum}.stim_response.raw_stim));
    raw_y = cellfun(@mean,(whisk_struct.window{cellNum}.stim_response.raw_response));
    
    x = whisk_struct.window{cellNum}.stim_response.bars_stim;
    y = whisk_struct.window{cellNum}.stim_response.bars_fit.mean;
    
    subplot(r(1),r(2),g)
    hold on; plot(ix,iy,'g')
    hold on; plot(x,y,'r')
    
        hold on; scatter(iraw_x,iraw_y,'g.');alpha(.25)
        hold on; scatter(raw_x,raw_y,'r.'); alpha(.25)
    
    if g == numel(both_tuned)
        legend('instantaneous','lagged window')
    end
    %     xlabel('whisker angle')
    %     ylabel('firing rate')
end

%% INTERSECTION OF PEAKS
i_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.instant(both_tuned),'uniformoutput',0)') ;
w_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],whisk_struct.window(both_tuned),'uniformoutput',0)');

figure(480);clf
errorbar(i_pw(:,1),w_pw(:,1),w_pw(:,2),w_pw(:,3),i_pw(:,2),i_pw(:,3),'ko','capsize',0)
if strcmp(hilbert_feature,'angle')
    hold on; plot([-30 80],[-30 80],'--k')
    set(gca,'xlim',[-30 80],'ylim',[-30 80],'xtick',-40:20:80,'ytick',-40:20:80)
elseif strcmp(hilbert_feature,'phase')
    hold on; plot([-pi pi],[-pi pi],'--k')
    set(gca,'xlim',[-4 4],'ylim',[-4 4],'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,...
        'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
end
xlabel('instantaneous preference')
ylabel('lagged window preference')
axis square

