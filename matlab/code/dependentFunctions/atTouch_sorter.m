% This function will find the indices of touch, the spikes around that
% touch given by the RANGE (ie -25:50ms) you provide. Lastly it'll find the
% theta, phase, amplitude, setpoint, max kappa, and pre touch velocity.

% Simple organization tool for visualizing parameters at touch onset

% OUTPUT: First 6 columns will be values for the variables 
% 1) THETA 2) AMP 3) SETPOINT 4) PHASE 5) MAX KAPPA 6) PRE TOUCH VELOCITY 
% Last columns will be the spikes around your given window 
% Will output values for ALL TOUCHES, PREDECISION TOUCHES, and POST
% DECISION TOUCHES

%INPUT: vector with time points you want to view (ie [-25:50])

function [tVar] = atTouch_sorter(array,range,preDecisionMat)

%find touch indices 
touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
spikes = squeeze(array.R_ntk);

touchOnIdx = sort(touchOnIdx(touchOnIdx<(numel(spikes)-range(end))));
touchOffIdx = sort(touchOffIdx(1:length(touchOnIdx)));

pretOnIdx = find(preDecisionMat==1);
postOnIdx = setdiff(touchOnIdx,pretOnIdx);

pretOffIdx = touchOffIdx(ismember(touchOnIdx,pretOnIdx));
posttOffIdx =  setdiff(touchOffIdx,pretOffIdx);

tOnIndices = {touchOnIdx,pretOnIdx,postOnIdx};
tOffIndices = {touchOffIdx,pretOffIdx,posttOffIdx};

%R_ntk spikes all around touches 
spikesIdx = cellfun(@(x) x+repmat(range,length(x),1),tOnIndices,'uniformoutput',0);
spikesAligned = cellfun(@(x) spikes(x),spikesIdx,'uniformoutput',0);

%S_ctk stimulus mat for features around touches 
varDesign = cell(1,3);
for g = 1:length(tOnIndices)
    touchTnums = ceil(tOnIndices{g}./array.t);
    for i = 1:length(tOnIndices{g})
        varDesign{g}(i,1)=array.S_ctk(1,tOnIndices{g}(i)); %theta at touch
        varDesign{g}(i,3)=array.S_ctk(3,tOnIndices{g}(i)); %amp at at touch
        varDesign{g}(i,4)=array.S_ctk(4,tOnIndices{g}(i)); %setpoint at touch
        varDesign{g}(i,5)=array.S_ctk(5,tOnIndices{g}(i)); %phase at touch
        kwin=array.S_ctk(6,tOnIndices{g}(i):tOffIndices{g}(i)); %get values in touch window
        [~ ,maxidx] = max(abs(kwin)); %find idx of max kappa within each touch window, neg or pos
        varDesign{g}(i,6)=kwin(maxidx); %use idx to pull max kappa
        varDesign{g}(i,2)=mean(array.S_ctk(2,tOnIndices{g}(i)-5:tOnIndices{g}(i)-1)); %finds mean of velocity (-4:-1ms) before touch
    end
    varDesign{g}(:,7) = array.meta.motorPosition(touchTnums);
end

%composing output structure
tVar.allTouches.varNames = {'theta','velocity','amp','midpoint','phase','curvature','motorPosition'};
tVar.allTouches.S_ctk = varDesign{1}; 
tVar.allTouches.R_ntk = spikesAligned{1};

tVar.preDecisionTouches.varNames = {'theta','velocity','amp','midpoint','phase','curvature','motorPosition'};
tVar.preDecisionTouches.S_ctk = varDesign{2};
tVar.preDecisionTouches.R_ntk = spikesAligned{2}; 

tVar.postDecisionTouches.varNames = {'theta','velocity','amp','midpoint','phase','curvature','motorPosition'};
tVar.postDecisionTouches.S_ctk = varDesign{3};
tVar.postDecisionTouches.R_ntk = spikesAligned{3}; 



