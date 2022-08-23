classdef orthoPoly1d
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        s (:,1);
        alpha (:,1);
        beta (:,1);
        variance (:,1);
    end
    
    methods
        function this = orthoPoly1d(k)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            this.s = nan(k+1,1);
            this.alpha = nan(k+1,1);
            this.beta = nan(k+1,1);
            this.beta(1) = 0;
            this.variance = nan(k+1,1);
        end
    end
end
