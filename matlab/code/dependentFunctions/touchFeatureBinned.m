%this function will be used to build a struct necessary for feature
%tuning quantification

function popV = touchFeatureBinned(U,viewWindow)

popV = cell(length(U),1);
%viewWindow : range around touch to view spikes (default -25:50ms) 

for rec=1:length(U)
    countThresh = 10; %min touch in each bin to be considered

    [tVar] = atTouch_sorter(U{rec},viewWindow);
    
    % Run this line below if you want to look at ret touches.
    %varspikes(varspikes(:,5)<0,:)=[];
    
    %touchOrder = fields(tVar);
    touchOrder = {'allTouches'};
    
    for g = 1:numel(touchOrder)
        %%
        fields = tVar.(touchOrder{g}).varNames;
        bounds = [{[-100:2:100]} {[-9750:125:9750]} {[-99.5:2.5:99.5]} {[-99.5:2.5:99.5]} {linspace(pi*-1,pi,13)} {[-.95:.05:.95]}];
        for d = [1 2 3 4 5 6] %for variables 1:6
            if d == 2
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).variables(:,d),tVar.allTouches.spikeMat,'equalE',numel(bounds{d})+1,-10000,10000);
            elseif d == 5
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).variables(:,d),tVar.allTouches.spikeMat,'equalX',numel(bounds{d})+1);
            elseif d == 6
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).variables(:,d),tVar.allTouches.spikeMat,'equalE',numel(bounds{d})+1,-1,1);
            else
                [sorted, sortedBy ,~]=binslin(tVar.(touchOrder{g}).variables(:,d),tVar.allTouches.spikeMat,'equalE',numel(bounds{d})+1,-100,100);
            end
            
            binrange = bounds{d};
            
            % Trimming unused ranges
            trims=[binrange' cell2mat(cellfun(@mean,sortedBy,'Uniformoutput',0))];
            indexHasTheta = ~isnan(trims(:,2));
            trims = trims(indexHasTheta, :);
            counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
            
            % Populating fields with values
            V.(fields{d}).bounds = bounds{d}; 
            V.(fields{d}).counts =counttmp(indexHasTheta==1,1);
            V.(fields{d}).range = trims(:,1);
            V.(fields{d}).spikes=cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
            V.(fields{d}).spikes=V.(fields{d}).spikes(indexHasTheta==1,:);
            V.(fields{d}).raw = sorted(indexHasTheta==1);
            
            %Trimming bins with touchcounts below thresholds
            mintouches=find(V.(fields{d}).counts<countThresh);
            V.(fields{d}).counts(mintouches,:)=[];
            V.(fields{d}).spikes(mintouches,:)=[];
            V.(fields{d}).range(mintouches,:)=[];
            V.(fields{d}).raw(mintouches,:) = [];
            
            popV{rec}=V;
            
        end
        
    end
end
