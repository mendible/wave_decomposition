function [Cin] = initialize_C(x,t,W,k,n,lambda,library)
% initA initializes the SR3 optimization using an iteratively thresholded
% least squares

Cin = zeros(k,n);

for jj = 1:n
    ttmp = t(logical(W(:,jj)));
    xtmp = x(logical(W(:,jj)));
    Ttmp = library(ttmp);
    Cincol = Ttmp\xtmp; % initial guess: Least-squares
    % lambda is our sparsification knob.
    for iter = 1:10
        smallinds = (abs(Cincol)<lambda); % find small coefficients
        Cincol(smallinds)=0; % and threshold
        biginds = ~smallinds(:);
        % Regress dynamics onto remaining terms to find sparse A
        Cincol(biginds) = Ttmp(:,biginds)\xtmp(:);
    end
    Cin(:,jj) = Cincol;
end

end

