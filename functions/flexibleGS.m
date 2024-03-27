function matched = flexibleGS(prefP,prefA)

    N = size(prefP,1);
    M = size(prefA,1);
    if N < M
        matched = zeros(M,1);
    else
        matched = zeros(N,1);
    end
    proposed = false(N,M);

    % while we still have free proposers ...
    while any(matched(1:N) == 0)
        % find them matches
        [matched,proposed] = findMatch(prefP,prefA,matched,proposed);
        % but stop if they run out of acceptors to propose to
        if all(all(proposed(matched(1:N) == 0,:),2))
            break
        end
    end
    if N < M
        matched = [(1:N)'  matched(matched>0)];
%         
%         (matched == 0) = setdiff(1:M,matched);
%         matched = [[(1:N)';0],matched];
    else
        matched = [(1:size(prefP,1))',matched];
    end
end