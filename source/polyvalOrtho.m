function y = polyvalOrtho(x,coeffTable)
    
    arguments
        x (:,1) double;
        coeffTable (:,4) table;
    end
    
    % initialize
    m = numel(x);
    y = zeros(m,1);
    k = numel(coeffTable.s)-1;
    p_iMinus1 = zeros(m,1);
    p_i = ones(m,1);
    
    for i =0:k
        y = y + coeffTable.s(i+1)*p_i;
       
        if i==k
            continue
        end
        % prepare for next loop
        p_iPlus1 = (x-coeffTable.alpha(i+2)).*p_i - coeffTable.beta(i+1)*p_iMinus1;
        p_iMinus1 = p_i;
        p_i =p_iPlus1;
    end
end

