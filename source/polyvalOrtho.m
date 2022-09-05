function y = polyvalOrtho(x,coefficients)
    % polyvalOrtho evaluates the orthogonal polynomial created from polyfitOrtho
    %
    %   y = polyvalOrtho(x,coefficients)
    %
    %% Inputz
    % x - double of any size. points to evaluate the polynomial. 
    % coefficients - Structure created from polyfitOrtho. 
    %
    %% Outputs
    % y - double. Same size as x. Values of polynomial evaluated at x.
    %
    %% Description
    % polyvalOrtho is complementary function to polyfitOrtho to evaluated the orthogonal polynomial created with
    % polyfitOrtho. 
    %
    % The time complexity of this function is O(k*m) where k is the order of the polynomial and m is the number of
    % points in x. 
    %
    %% Example
    % % Fit a 186th degree polynomial to the Runge function on 500 Chebyshev points. The Chebyshev points are chosen
    % % because high degree interpolating polynomials are known to converge to the Runge function.
    %   runge = @(x) 1./(1+25*x.^2);
    %   theta = linspace(pi,0,500);
    %   x = cos(theta); % Chebyshev points
    %   y = runge(x);
    %   k = 186;
    %   coefficientsAndResults = polyfitOrtho(x,y,k);
    %   fplot(@(x) runge(x)-polyvalOrtho(x,coefficientsAndResults),[-1 1])
    %   ylabel('Error')
    %   figure;
    %   fplot(@(x) runge(x),[-1 1])
    %   hold all;
    %   fplot(@(x) polyvalOrtho(x,coefficientsAndResults),[-1 1])
    %   legend(["Runge Function", "100th degree polynomial"],'location','best')
    %
    % See also polyfitOrtho, polyfit, polyval
    
    arguments
        x double;
        coefficients struct;
    end
    
    % initialize
    sz = size(x);
    x = x(:);
    m = numel(x);
    y = zeros(m,1);
    k = numel(coefficients.s)-1;
    p_iMinus1 = zeros(m,1);
    p_i = ones(m,1);
    
    for i =0:k
        y = y + coefficients.s(i+1)*p_i;
       
        if i==k
            continue
        end
        % prepare for next loop
        p_iPlus1 = (x-coefficients.alpha(i+2)).*p_i - coefficients.beta(i+1)*p_iMinus1;
        p_iMinus1 = p_i;
        p_i =p_iPlus1;
    end
    
    y = reshape(y,sz);
end

