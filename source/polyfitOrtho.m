function [coefficientsAndResults] = polyfitOrtho(x_mu,f_mu,k)
    %polyfitOrtho fits orthogonal polynomials in a least squares sense
    %
    %   coefficientsAndResults = polyfitOrtho(x_mu,f_mu,k)
    %
    %% Inputs
    % x_mu - vector of x points for fit. m = numel(x_mu) by definition.
    % f_mu - vector of y points for fit. f_mu m elements.
    % k - polynomial fit order. This algorithm requires k<=m-1. 
    %
    %% Outputs
    % coefficientsAndResults - A table containing the coefficients, alphas, betas, and variances'. The table always has k+1 rows.
    %     s - is the cofficients of each polynomial.
    %     alpha - The 1st coefficient it the orthogonality procedure.
    %     beta - The 2nd coefficient in the orthogonality procedure.
    %     variance - The variance of the residuals for each polynomial fit. 
    %
    %% Description
    % This functions fits a polynomial of the form.
    %
    %   y_k(x) = s_0*p_0(x) + s_1*p_1(x) + ... + s_k*p_k(x)
    %
    %   where 
    %
    %   - p_i(x) is a polynomial of order i. By definition, p_i is orthogonal to p_j as long as i~=j on the inner 
    %     product, sum(p_j(x_mu).*p_i(x_mu)). 
    %   - s_i is the coefficient of p_i. - y_k(x) is then a polynomial of order k constructed from polynomials that are
    %     orthogonal on the inner product sum(p_j(x_mu).*p_i(x_mu)).
    %
    %   The orthogonal polynomials, p_i, are generated using the following definition:
    %
    %   p_0(x) = 1;
    %   p_1(x) = x*p_0(x) - alpha_1*p_0(x);
    %   p_2(x) = x*p_1(x) - alpha_2*p_1(x) - beta_1*p_0(x);
    %   ...      ...
    %   p_(i+1)(x) = x*p_i(x) - alpha_(i+1)*p_i(x) - beta_i*p_(i-1)(x);
    %
    %   Alpha and Beta are chosen to make the orthogonality relations hold. See Forsythe, [1] for derivation. The
    %   formulas for alpha's and beta's are
    %  
    %   alpha_(i+1) = sum(x_mu.*p_i(x_mu).^2)./sum(p_i(x_mu).^2); 
    %   beta(i) = sum(x_mu.*p_i(x_mu).*p_(i-1)(x_mu))./sum(p_(i-1)(x_mu).^2);
    %   
    %   % Note that more involved, less intuitive formulas are used in the code. These formulas can be found in 
    %     Forsythe, [1].
    %
    %  The beauty of this formulation is that the least squares problem is diagonal:
    %
    %  w_ii = sum(p_i(x_mu).*p_i(x_mu));
    %  omega_i = sum(f_mu.*p_i(x_mu);
    %  
    %  w_00*s_0                  = omega_0
    %         w11*s_1            = omega_1
    %                ...
    %                   w_kk*s_k = omega_k 
    %
    %  
    %   Other important properties:
    %   - The algorithm does not require forming the A matrix in the least squares problem Ax=b. 
    %   - The algorithm finds each coefficient s_i, in order from 0 to k. This has the property that you get all
    %     polynomial fits from 0 to k by requesting kth polynomial fit. Use 0 to k2 rows where k2<k.
    %   - The variance of the residuals is calculated for each i from 0 to k by an update formula.
    %   - The least square problem is better formed than Vandermonde matrices (i.e. x.^(0:k), aka monomial basis). This
    %     allows fitting of much higher order polynomials.
    %
    %% Example
    % % Fit a 100th degree polynomial to the Runge function on 150 Chebyshev points. Compare to polyfit/Vandermonde %
    % matrices. The monomial basis, x^(0:n), is exponentially ill-conditioned as n increases. polyfit/Vandermonde cannot
    % % fit a 100th degree polynomial. The Chebyshev points are chosen because high degree interpolating polynomials are
    % % known to converge to the Runge function on this set of points (this is not interpolation though but rather least
    % % sqaures fitting).
    %   runge = @(x) 1./(1+25*x.^2);
    %   theta = linspace(pi,0,150);
    %   x = cos(theta); % Chebyshev points
    %   y = runge(x);
    %   k = 100;
    %   coefficientsAndResults = polyfitOrtho(x,y,k);
    %
    %
    %% References
    % The algorithm and nomenclature comes from Forsythe 1957.
    % [1] Forsythe, George E. "Generation and use of orthogonal polynomials for data-fitting with a digital computer."
    % Journal of the Society for Industrial and Applied Mathematics 5.2 (1957): 74-88.
    %
    % See also polyvalOrtho, polyfit, polyval.
    %
    
    %% Input Checking
    arguments
        x_mu (:,1) double;
        f_mu (:,1) double;
        k (1,1) double {mustBeInteger, mustBePositive};
    end
    
    m = numel(x_mu);
    assert(m==numel(f_mu),"x and f must have the same number of elements.");
    
    %% Algorithm
    % This algorithm is 
    
    % initialize
    p_iMinus1 = zeros(m,1); % p_(-1)
    p_i = ones(m,1); % p_0
    deltaSquarediMinus1 = f_mu'*f_mu;
    coefficientsAndResults = table(nan(k+1,1), nan(k+1,1), nan(k+1,1), nan(k+1,1),'VariableNames',["s","alpha","beta","variance"]);
    coefficientsAndResults.beta(1) = 0; % beta_0
    coefficientsAndResults.Row = "" + (0:k)';
    wii = m;
    
    for i=0:k
        omegai = f_mu'*p_i;
        coefficientsAndResults.s(i+1) = omegai/wii; % s_i
        deltaSquaredi = deltaSquarediMinus1 - coefficientsAndResults.s(i+1)^2*wii;
        coefficientsAndResults.variance(i+1) = deltaSquaredi./(m-i-1);
        % TODO. We could test here if the variance is low enough.
        
        if i>=k
            break
        end
        
        coefficientsAndResults.alpha(i+2) = (x_mu.*p_i)'*p_i/wii; % alpha_(i+1)
        
        % Prepare for next loop
        p_iPlus1 = (x_mu-coefficientsAndResults.alpha(i+2)).*p_i - coefficientsAndResults.beta(i+1)*p_iMinus1;
        p_iMinus1 = p_i;
        p_i = p_iPlus1;
        wiPlus1iPlus1 = p_iPlus1'*p_iPlus1;
        coefficientsAndResults.beta(i+2) = wiPlus1iPlus1/wii; % beta_(i+1)
        wii = wiPlus1iPlus1;
        deltaSquarediMinus1 = deltaSquaredi;
        
    end
    
end

