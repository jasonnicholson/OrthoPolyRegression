function y = orthoPolyEval(x,coeffTable)
    
    arguments
        x (:,1) double;
        coeffTable (:,3) table;
    end
    
    % initialize
    m = numel(x);
    y = zeros(m,1);
    k = numel(coeffTable.s)-1;
    piMinus1 = zeros(m,1);
    pi = ones(m,1);
    
    for i =0:k
        y = y + coeffTable.s(i+1)*pi;
       
        if i==k
            continue
        end
        % prepare for next loop
        piPlus1 = (x-coeffTable.alpha(i+2)).*pi - coeffTable.beta(i+1)*piMinus1;
        piMinus1 = pi;
        pi =piPlus1;
    end
end

