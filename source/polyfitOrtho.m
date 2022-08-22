function [resultsTable] = polyfitOrtho(x,f,k)
    %polyfitOrtho fits polynomial orthogonal
    %   resultsTable = polyfitOrtho(x,f,k)
    %
    %% Inputs
    % x - vector of x points for fit.
    % f - vector of y points for fit. Must be the same size as x.
    %
    %% Outputs
    % resultsTable - 
    %% Description
    %
    %% Example
    % Fit a 100th degree polynomial to the Runge function on 150 Chebyshev points. Compare to polyfit/Vandermonde
    % matrices. The monomial basis, x^(0:n), is exponentially ill-conditioned as n increases. polyfit/Vandermonde cannot
    % fit 100th degree polynomial.
    %   runge = @(x) 1./(1+25*x.^2);
    %   theta = linspace(pi,0,150);
    %   x = cos(theta); % Chebyshev points
    %   y = runge(x);
    %   k = 100;
    %   resultsTable = polyfitOrtho(x,y,k);
    %
    %
    %% References
    % The algorithm and nomenclature comes from Forsythe 1957.
    % Forsythe, George E. "Generation and use of orthogonal polynomials for data-fitting with a digital computer."
    % Journal of the Society for Industrial and Applied Mathematics 5.2 (1957): 74-88.
    %
    % See also polyvalOrtho, polyfit, polyval.
    %
    
    %% Input Checking
    arguments
        x (:,1) double;
        f (:,1) double;
        k (:,1) double {mustBeInteger, mustBePositive};
    end
    
    m = numel(x);
    assert(k<m-1,"k < m-1 is required for this algorithm.")
    
    assert(m==numel(f),"x and f must have the same number of elements.");
    
    %% Algorithm
    % This algorithm is 
    
    % initialize
    p_iMinus1 = zeros(m,1); % p_(-1)
    p_i = ones(m,1); % p_0
    deltaSquarediMinus1 = f'*f;
    resultsTable = table(nan(k+1,1), nan(k+1,1), nan(k+1,1), nan(k+1,1),'VariableNames',["s","alpha","beta","variance"]);
    resultsTable.beta(1) = 0;
    resultsTable.Row = "" + (0:k)';
    wii = m;
    
    for i=0:k
        omegai = f'*p_i;
        resultsTable.s(i+1) = omegai/wii;
        deltaSquaredi = deltaSquarediMinus1 - resultsTable.s(i+1)^2*wii;
        resultsTable.variance(i+1) = deltaSquaredi./(m-i-1);
        % TODO. We could test here if the variance is low enough.
        
        if i>=k
            break
        end
        
        resultsTable.alpha(i+2) = (x.*p_i)'*p_i/wii;
        
        % Prepare for next loop
        p_iPlus1 = (x-resultsTable.alpha(i+2)).*p_i - resultsTable.beta(i+1)*p_iMinus1;
        p_iMinus1 = p_i;
        p_i = p_iPlus1;
        wiPlus1iPlus1 = p_iPlus1'*p_iPlus1;
        resultsTable.beta(i+2) = wiPlus1iPlus1/wii;
        wii = wiPlus1iPlus1;
        deltaSquarediMinus1 = deltaSquaredi;
        
    end
    
end

