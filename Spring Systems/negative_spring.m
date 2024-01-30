classdef negative_spring
    %Negative_spring describes a negative spring system with 4 leaf springs
    %   Optimizes the parameters to create desired stiffness
    %   Detailed explanation goes here
    
    properties
        %Constant properties 
        E           %[Pa] Elastic modulus spring
        sigma_y     %[Pa] Yield stress spring
        t           %[m] Plate thickness spring
        b           %[m] Width spring
        n           %[] number of springs
        sm = 0.87   %[] Safety margin on stroke from JPE
        thread_pitch = 0.7 %[mm\revolution] 
        
        
        %Computed properties
        I           %[m^4] moment of inertia
        L           %[m] length of the leaf springs
        k           %[N/m] stiffness
        uz_max      %[m] compresion in the z direction fo the leaf spring
        S           %[m] maximal stroke of the spring
        n_rotation  %[] number the thread needs to make
        
    end
    
    methods
        function obj = negative_spring(E,sigma_y,t,b,n,k_des)
            %Constructing method
            obj.E = E;
            obj.sigma_y = sigma_y;
            obj.t = t;
            obj.b = b;
            obj.n = n;
            obj = obj.comp_properties(k_des);
        end
        
        function obj  = comp_properties(obj, k_des)
            %Computes the desired spring properties given a desired
            %stiffness
            obj.I = obj.t^3*obj.b/12; %[m^4]  moment of inertia leaf spring
            Length = ((44.4*obj.E*obj.I*obj.n)/k_des)^(1/3); %desired length given desired
            obj.L = round(Length,4); %ounding the result
            obj = obj.comp_stiff();
            obj.uz_max = obj.sigma_y*obj.L^2/(53*obj.E*obj.t); %[m] max compression in z direction leaf spring
            obj.n_rotation = obj.uz_max/(obj.thread_pitch*1e-3); %number of rotations
            obj = obj.comp_stroke();
        end
        
        function obj = comp_stiff(obj)
            %computes system stiffness given the system parameters
            obj.k = round(44.4*obj.E*obj.I/obj.L^3 * obj.n );
        end
        
        function obj = comp_stroke(obj)
            %computes system maximal stroke
            obj.S = obj.sm * sqrt(5/3*obj.uz_max*obj.L); %[m] max stroke of springs with safety margin
        end
        
    end
end

