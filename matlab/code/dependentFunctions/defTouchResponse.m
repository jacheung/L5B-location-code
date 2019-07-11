function U = defTouchResponse(U,confidenceThreshold,displayOpt)

%function that builds upon location code data struct. Adds fields of
%responseWindow under meta. Response window is defined by a confidence
%threshold above baseline response (BL response is from -25:0ms pre touch)

if (nargin < 3), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    rc = numSubplots(length(U));
    figure(3000);clf
end

for rec=1:length(U)
    array = U{rec};
    window = [-25:50];
    touchIdx= [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    spks = squeeze(array.R_ntk);
    
    within_range = ~ (logical(sum((touchIdx + window) > numel(spks),2)) | logical(sum((touchIdx + window) < 0, 2)));
    touchIdx = touchIdx(within_range); 
    
    
    blIdx = window(find(window==-25):find(window==0));
    
    touchSpks = spks(touchIdx+window);
    touchSpksShuff = spks(touchIdx+blIdx);
    
    %calculating x% confidence interval
    SEM = nanstd(touchSpksShuff(:))./ sqrt(sum(~isnan(touchSpksShuff(:,1))));
    ts = tinv(confidenceThreshold,sum(~isnan(touchSpksShuff(:,1))));
    CI = SEM.*ts;
    
    touchResponse = smooth(mean(touchSpks));
    excitThreshold = mean(touchSpksShuff(:)) + CI;
    inhibThreshold = mean(touchSpksShuff(:)) - CI;
    
    %Excitatory threshold defined as the period > mean+x%CI
    excitthreshIdx = touchResponse'>excitThreshold;
    excitthreshIdx(1:find(window==0))=0;
    
    %Inhibitory threshold defined as the period > mean-x%CI
    inhibthreshIdx = touchResponse'<inhibThreshold;
    inhibthreshIdx(inhibthreshIdx==1)= -1;
    inhibthreshIdx(1:find(window==0))=0;
    
    if willdisplay
        figure(3000);subplot(rc(1),rc(2),rec)
        hold on; scatter(window,touchResponse*1000,'k')
        hold on; plot(window,ones(length(window),1).* excitThreshold .* 1000,'r-.')
    end
    
    %Defining touch responses as a period between two points that are less
    %than 5ms apart.
    if sum(excitthreshIdx)>0
        tps = window(excitthreshIdx);
        if ~isempty(tps)
            startPoint = tps(1);
            endPoint = tps(find(diff(diff(tps)<5)==-1,1,'first'));
            if isempty(endPoint)
                endPoint = tps(end);
            end
            
            if willdisplay
                hold on; scatter(window(startPoint+find(window==0):endPoint+find(window==0)),touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000,'b','filled')
            end
            
            %used below to eliminate touch responses that are sig but way
            %too small. 
            if mean(touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000) > 2
                U{rec}.meta.responseWindow=[startPoint endPoint];
            end
            
        end
    end
    
end