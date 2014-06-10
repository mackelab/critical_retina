classdef ProbDist < handle
    %PROBDISTRIBUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private);
        sup_x
        pdf_x
    end
    
    methods
        function obj = ProbDist(x, pdf_x)
            obj.sup_x = x(:);
            f = trapz(x, pdf_x);
            obj.pdf_x = pdf_x(:) / f;
        end
        
        function y = expectation(obj)
            y = trapz(obj.sup_x, obj.sup_x .* obj.pdf_x);
        end
        
        function y = variance(obj)
            err = (obj.sup_x - obj.expectation()).^2;
            y = trapz(obj.sup_x, err .* obj.pdf_x);
        end
        
        function y = stdev(obj)
            y = sqrt(obj.variance);
        end
        
        function y = eval_pdf(obj, x)
            y = interp1(obj.sup_x, obj.pdf_x, x, 'pchip');
        end
        
        function y = moment(obj, order)
            x = obj.sup_x;
            y = trapz(obj.sup_x, (x .^ order) .* obj.pdf_x);
        end
    end
    
end

