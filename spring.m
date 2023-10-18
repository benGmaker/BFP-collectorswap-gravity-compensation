classdef spring
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        k %stiffness [N/m]
        x0 = [0; 0] %position base [x y]
        x1 = [0; 0] %postion end effector 
        L0 %initial length
        r = [1; 0]%direction
        L %length
    end

    methods
        function obj = spring(k,x0,L0,r)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.k = k;
            obj.x0 = x0;
            obj.x1 = x0 + r*L0;
            obj.L0 = L0;
            obj.r = r;
            obj.L = obj.L0;
        end
        
        function obj = setpos(obj,x1)
            obj.x1 = x1;
        end

        function obj = setLength(obj)
            obj.L = sum(abs(obj.x1 - obj.x0));
        end

        function obj = elongation(obj, F)
            obj.L = sum(F) / obj.k + obj.L0;
            obj.x1 = obj.x0 + obj.r*obj.L;
        end

        function F = F(obj)
            %Output spring given the elenogation of the spring
            obj = obj.setLength();
            F = obj.k * (obj.L - obj.L0);
        end
    end
end