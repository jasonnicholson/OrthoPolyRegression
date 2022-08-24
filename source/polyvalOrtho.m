function y = polyvalOrtho(x,coefficients)
    
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

