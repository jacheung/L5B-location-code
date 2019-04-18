function tunedCells = tuningQuantification(cellStruct,popV,selectedCells,variableFields,touchOrderFields,window)

%% ZSCORE CONVERSION and tuning width

[rc] = numSubplots(length(selectedCells));
% thetabounds = linspace(pi*-1,pi,13);
bounds = popV{1}.allTouches.(variableFields{1}).bounds;
plotrow = rc(1);
plotcolumn = rc(2);

for g = 1:numel(touchOrderFields)
    
    fr= nan(length(selectedCells),2);
    pw=fr;
    pwzs = nan(length(selectedCells),50);
    pwtheta = nan(length(selectedCells),2);
    figure(390+g);clf
    numInterpPts = 16;
    
    for d = 1:length(selectedCells)
        
        %ZSCORREEEEE
        currCell = selectedCells(d);
        
        pOn = round(cellStruct{currCell}.meta.poleOnset(1)*1000);
        restW(d) = nanmean(nanmean(cellStruct{currCell}.S_ctk(1,1:pOn,:)));
        
        counts = popV{currCell}.(touchOrderFields{g}).(variableFields{1}).counts;
        ranges = popV{currCell}.(touchOrderFields{g}).(variableFields{1}).range ;
        thetavect = [];
        for k = 1:length(ranges)
            thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
        end
        tvectNorm = normalize_var(thetavect,-1,1);
        
        allspikes = cell2mat(popV{currCell}.(touchOrderFields{g}).(variableFields{1}).raw);
        meanbase = mean(mean(allspikes(:,1:find(window==0)),2));
        stdbase = std(mean(allspikes(:,1:find(window==0)),2));
        %     postresponses = allspikes(:,find(window==0)+5:find(window==0)+35);
        postresponses = allspikes(:,find(window==0)+ cellStruct{currCell}.meta.responseWindow(1) : find(window==0)+cellStruct{currCell}.meta.responseWindow(2));
        
        xresp = mean(postresponses,2);
        zscore = (xresp - meanbase) ./ stdbase;
        
        for k = 1:size(postresponses,1)
            spkT = find(postresponses(k,:)==1);
            if ~isempty(spkT) & numel(spkT)>1
                ISI = diff(spkT);
                CV(k) = std(ISI)./mean(ISI);
            else
                CV(k) = 0;
            end
        end
        
        [sorted,~,~] = binslin(thetavect,zscore,'equalE',numel(bounds),bounds(1),bounds(end));
        [CVsorted,~,~] = binslin(thetavect,CV','equalE',numel(bounds),bounds(1),bounds(end));
        
        
        samps = cellfun(@numel,sorted);
        selBins = find(samps>sum(samps)./100);
        
        for p=1:length(selBins)
            CVpop(p,:) = [bounds(selBins(p)) mean(CVsorted{selBins(p)})];
        end
        
        
        zraw = nan(length(sorted),2000);
        for b = 1:length(sorted)
            if sum(b == selBins)>0
                currvals = sorted{b};
                if ~isempty(sorted{b})
                    zraw(b,1:length(currvals)) = currvals';
                end
            end
        end
        
        
        [p,~,stats]=anova1(zraw',[],'off');
        vals =multcompare(stats,[],'off');
        popzscore{d} = cellfun(@mean,sorted);
        otune(d)=p;
        
        
        x=cellfun(@str2num,stats.gnames);
        
        %CI BINS
        cibins = nan(size(x,1),1);
        
        SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
        for i = 1:length(x)
            rawx = zraw(x(i),:);
            SEM = SEMg(x(i));
            ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
            cibins(i,:) = ts.*SEM;   %confidence intervals
        end
        
        
        y = stats.means;
        %     ytune(x,d) = y;  %FOR PHASE
        
        x=bounds(x)';
        
        
        zstd = nanstd(zraw,[],2);
        ystd = zstd(~isnan(zstd))';
        
        xy = [x y'];
        
        kxy{d} = xy;
        
        
        subplot(plotrow,plotcolumn,d);
        shadedErrorBar(xy(:,1),xy(:,2),cibins,'-k');
        %     scatter(xy(:,1),xy(:,2),[],[.7 .7 .7],'filled')
        set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
            'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
        
        nxy = [normalize_var(x,0,1),y'];
        iy(:,d) = interp1(nxy(:,1),nxy(:,2),linspace(0,1,numInterpPts));
        
        modely(:,d) = interp1(normalize_var(x,0,1),y,linspace(0,1,numInterpPts));
        modelystd(:,d) = interp1(normalize_var(x,0,1),ystd,linspace(0,1,numInterpPts));
        
        axis square
        
        zs = nanmean(zraw,2);
        zs(isnan(zs))=[];
        
        sigs =  vals(vals(:,end)<.01,:);
        
        clear peaksig
        [~,maxidx]=max(zs);
        midxsel = find(sum(sigs(:,[1 2])==maxidx,2)==1);
        
        if sum(midxsel>0)
            peaksig = sigs(midxsel,1:2);
            
            [~,midx] = min(abs(diff(peaksig')));
            pw(d,:) = peaksig(midx,:);
            figure(390+g);
            hold on; scatter(xy(peaksig(midx,:),1),xy(peaksig(midx,:),2),[],'filled','g')
            
            raws=zs(pw(d,1):pw(d,2));
            pwzs(d,1:length(raws)) = raws;
            pwtheta(d,:) = x(pw(d,:));
        end
        
    end
    
    
    ocellidx = find(~isnan(pw(:,1)));
    tunedCells{g} = selectedCells(ocellidx);
end
