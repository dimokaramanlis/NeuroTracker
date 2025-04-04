
function [mcorr, imax]= maxSlidingCorr(temp1, temp2, maxLag)

[Nyx, Nt]= size(temp1);

lag = -maxLag:maxLag;
XCmat = zeros(2, 2, length(lag), 'single');
spike_templates = reshape(cat(3,temp1, temp2), Nyx, Nt, 2);

for ilag = 0:maxLag
    %shift matrices

    Y1 = reshape(spike_templates(:, 1+ilag:end, :), Nyx * (Nt - ilag), 2); 
    Y2 = reshape(spike_templates(:, 1:end-ilag, :), Nyx * (Nt - ilag), 2);
    resmat = corr(Y1,Y2); %core calculation
    
    XCmat(:,:, maxLag+1+ilag) = resmat; 
    XCmat(:,:, maxLag+1-ilag) = resmat';
end

[mcorr, im] = max(XCmat(1,2,:));
imax = lag(im);

end