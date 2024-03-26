function [stableMatching] = galeShapley(menPref, womenPref)
% GALESHAPLEY Solves the stable marriage problem using the Gale-Shapley algorithm.
% Handles cases with unequal numbers of men and women.

numMen   = size(menPref, 1);
numWomen = size(womenPref, 1);
largerGroupSize = max(numMen, numWomen);

% Pad the smaller group with dummy preferences
if numMen < numWomen 
    menPref = [menPref; zeros(numWomen - numMen, largerGroupSize)]; 
    for i = numMen+1 : largerGroupSize 
        menPref(i, 1:numWomen) = 1:numWomen;
        menPref(i, numWomen+1:largerGroupSize) = largerGroupSize:-1:numMen+1;
    end
elseif numWomen < numMen
    womenPref = [womenPref; zeros(numMen - numWomen, largerGroupSize)]; 
    for i = numWomen+1 : largerGroupSize 
        womenPref(i, 1:numMen) = 1:numMen;
        womenPref(i, numMen+1:largerGroupSize) = largerGroupSize:-1:numMen+1; % Corrected line
    end
end

% Initialize 
freeMen = 1:largerGroupSize;  % Note: might include dummy men now
engagedWomen = zeros(largerGroupSize, 1); 
women_suitor = zeros(largerGroupSize, largerGroupSize); 

% Precalculate women's ranking of men for efficiency
rank = zeros(largerGroupSize, largerGroupSize);
for i = 1:largerGroupSize
    for j = 1:largerGroupSize
        rank(i, womenPref(i, j)) = j; 
    end
end

% While there are still free men
while ~isempty(freeMen)
    m = freeMen(1);

    nextProposal = find(menPref(m, :) > 0, 1); 
    if ~isempty(nextProposal) 
        w = menPref(m, nextProposal);
        menPref(m, nextProposal) = 0; 
        women_suitor(w, m) = m; 

        if engagedWomen(w) == 0  
            engagedWomen(w) = m; 
        else  
            currentPartner = engagedWomen(w);
            if rank(w, m) < rank(w, currentPartner) 
                engagedWomen(w) = m;  
                freeMen = [freeMen currentPartner];
            end
        end
    else
        freeMen(1) = []; 
    end
end

% Construct the final matching (filter out dummy matches)
stableMatching = zeros(largerGroupSize, 2); 
matchCount = 0; 
for m = 1:largerGroupSize
    if engagedWomen(m) ~= 0 && m <= numMen % Filter out dummy men
        matchCount = matchCount + 1;
        stableMatching(matchCount, 1) = m;
        stableMatching(matchCount, 2) = engagedWomen(m);
    end
end

stableMatching = stableMatching(1:matchCount, :); 
end