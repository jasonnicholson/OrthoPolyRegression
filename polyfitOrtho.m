function [coeffTable] = orthoPolyFit(x,f,k)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    arguments
        x (:,1) double;
        f (:,1) double;
        k (:,1) double {mustBeInteger, mustBePositive};
    end
    
    m = numel(x);
    
    
    % initialize
    piMinus1 = zeros(m,1);
    pi = ones(m,1);
    wii = m;
    beta = nan(k+1,1);
    beta(1) = 0;
    deltaSquarediMinus1 = f'*f;
    alpha = nan(k+1,1);
    s = nan(k+1,1);
    
    coeffTable = table(s,alpha,beta,variance);
    coeffTable.Row = "" + (0:k)';
    
    for i=0:k
        omegai = f'*pi;
        coeffTable.s(i+1) = omegai/wii;
         deltaSquaredi = deltaSquarediMinus1 - coeffTable.s(i+1)^2*wii;
%         variancei = deltaSquaredi./(m-i-1);
        % TODO. We could test here if the variance is low enough.
        
        if i>=k
            continue
        end
                
        coeffTable.alpha(i+2) = (x.*pi)'*pi/wii;
        
        % Prepare for next loop
        piPlus1 = (x-coeffTable.alpha(i+2)).*pi - coeffTable.beta(i+1)*piMinus1;
        piMinus1 = pi;
        pi = piPlus1;
        wiPlus1iPlus1 = piPlus1'*piPlus1;
        coeffTable.beta(i+2) = wiPlus1iPlus1/wii;
        wii = wiPlus1iPlus1;
        deltaSquarediMinus1 = deltaSquaredi;
        
    end
    
end

