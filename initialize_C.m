function [Ain] = initA(x,t,W,k,n,lambda,library)
% initA initializes the SR3 optimization using an iteratively thresholded
% least squares

Ain = zeros(k,n);

for jj = 1:n
    ttmp = t(logical(W(:,jj)));
    xtmp = x(logical(W(:,jj)));
    Ttmp = library(ttmp);
    Aincol = Ttmp\xtmp; % initial guess: Least-squares
    % lambda is our sparsification knob.
    for iter = 1:10
        smallinds = (abs(Aincol)<lambda); % find small coefficients
        Aincol(smallinds)=0; % and threshold
        biginds = ~smallinds(:);
        % Regress dynamics onto remaining terms to find sparse A
        Aincol(biginds) = Ttmp(:,biginds)\xtmp(:);
    end
    Ain(:,jj) = Aincol;
end

end

