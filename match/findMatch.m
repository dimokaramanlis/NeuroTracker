function [matched,proposed] = findMatch(prefP,prefA,matched,proposed)
    p = find(matched == 0,1);
    % with |p| who still has acceptors to propose to
    acceptors = prefP(p,~proposed(p,:));

    for ii = 1:length(acceptors)
        % |a| is |p|'s highest ranked such acceptor to whom |p| has not yet
        % proposed
        a = acceptors(ii);
        proposed(p,prefP(p,:) == a) = true;

        % if |a| is free
        if ~any(matched == a)
            % (p, a) become engaged
            matched(p) = a;
            break
        else
            % else some pair (p', a) already exists
            p_ = find(matched == a);
            [ranking,~] = find(prefA(a,:) == [p;p_]);

            % if |a| prefers |p| to |p'|
            if  ranking(1) < ranking(2)
                % |p'| becomes free
                matched(p_) = 0;
                % (p, a) become engaged
                matched(p) = a;
                break
            else
                % else (p', a) remain engaged
            end
        end
    end
end