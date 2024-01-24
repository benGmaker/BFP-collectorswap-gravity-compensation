classdef springsystem
    %Springsystem describes a system including construction with two loading conditions F1 and F2
    %It has functions to optimize the defining parameters for the given
    %desired input variables
    
    properties
        %INPUT VARIABLES
        %input spring variables
        F1 = 2          %[N] initial loading condition (HIGH)
        F2 = 1          %[N] second loading condition (LOW)
        SR = 2          %[] spring elongation ratio
        S               %[m] desired / build stroke 
        
        %input construction
        h_max           %[m] maximal build height
        h_mech          %[m] height needed for construction and 
        R1 = 0          %[m] inner diameter pulley
       
        %real spring input
        name            %name of the spring used
        n               %[] number of springs used
        Fn_tot          %[N] force maximal elongation spring
        fn              %[mm] maximal extension of the spring
        
        %RESULT
        %general results
        max_stroke = 0  %[m] maximal stroke of the system
        h_adjust        %[m] height needed to achieve desired stroke & mass transfer
        R2 = 0          %[m] outer diameter pulley 
        
        %spring results
        k = 0           %[N/m] stiffness
        L0              %[m] initiall spring length
        L1              %[m] initial spring elongation (HIGH) (excl. L0)
        L2              %[m] secondary spring elongation (LOW) (excl. L0)
        Lmax            %[m] maximal spring elongation incl L0
        springstroke    %[m] stroke possible by spring
    end
    
    methods
        function obj = springsystem(F1, F2, S, h_mech, h_max, SR, R1)
            %Constructing method
            %Set R1 to zero if no spiral pulley is desired in the system
            obj.F1 = F1;
            obj.F2 = F2;
            obj.S = S;
            obj.h_mech = h_mech;
            obj.h_max = h_max;
            obj.SR = SR;
            obj.R1 = R1;
        end
      
        function obj = comp_stiffness(obj,decimals)
            %Computes the lowest possible stiffness given system parameters 
            %the stiffnes is rounded by with n decimals
            %Uses SR the estimate the extension of the spring
            F = ((1+1/(1-obj.SR) )*obj.F1-obj.F2); 
            d = (obj.h_max-obj.R2-obj.h_mech-obj.S/(1-obj.SR)); 
            obj.k =  round(F/d,decimals); %results is rounded
        end
        
        function obj = comp_stiffnessAndPulley(obj,decimals)
            %optimizes system stiffness and pulley dimensions
            
            %setup
            obj = obj.comp_stiffness(decimals); %compute stiffness
            if obj.R1 == 0 %no inner diameter is set, meaning system has no pulley 
                return
            end
            
            %optimizing R2
            for i = 1:100 %doing at most 100 itteration
                R = obj.F1/obj.k; %computing new outer diameter
                if abs(R - obj.R2) < 1e-5 % 
                    %if the new value is close to previous one a optimimum
                    %is reached
                    obj.R2 = R;
                    return
                end
                obj.R2 = R;
                obj = obj.comp_stiffness(decimals); % compute new parameters 
            end
        end
        
        function obj = compute_lengths(obj)
            %The lengths can be computed if the stiffness is known
            obj.h_adjust = (obj.F1 - obj.F2)/obj.k; %adjustment height
            obj.Lmax = obj.h_max - obj.R2- obj.h_mech - obj.h_adjust; %maximum length of the spring
            obj.L1 = obj.F1/obj.k; %initial elongation (without initial length)
            obj.L2 = obj.F2/obj.k; %second elongation (without initial length)
        end
        
        function height_check(obj)
            %Validation function to check if heights of the system match
            %the max height defined
            heigh_sum = round(obj.R2 + obj.S + obj.L1 + obj.L0 + obj.h_adjust + obj.h_mech,4);
            if ((heigh_sum == obj.h_max)) == 0 
                display('max height and sum of heights do not match!')  
                display('maximal height = ' +  string(obj.h_max))
                display('sum of heights = ' + string(heigh_sum))
            end
        end
        
        function obj = comp_sys(obj,decimals)
            %Given the system parameters computes the optimal system
            %properties
            obj = obj.comp_stiffnessAndPulley(decimals);
            obj = obj.compute_lengths();
            obj.L0 = obj.Lmax*obj.SR; %initial length spring
            obj.height_check()
        end
        
        function desired_properties(obj)
            %Displays the desired spring properties for spring selection   
            fn = (obj.S + obj.L1)*1e3;
            Fn = (fn*obj.k)*1e-3;
            
            display("k = " + string(obj.k));
            display('l0 = ' + string(obj.L0*1e3));
            display('fn = ' + string(fn));
            display('Fn = ' + string(Fn));
        end
        
        function obj = real_spring_properties(obj, l0, fn, Fn, n)
            %computes the properties of the stage given the properties of a
            %real spring system with 
            %l0 initial length [mm]
            %fn extension length [mm]
            %Fn force at max elongation [N]
            %N numbert of springs
            obj.L0 = l0*1e-3;
            obj.Fn_tot = Fn*n;  %effective max spring output force
            obj.fn = fn;
            obj.k = round(obj.Fn_tot/fn*1e3);     %stiffness of spring system
            if obj.R1 > 0 %computing pulley size if there is a pulley
                obj.R2 = obj.F1/obj.k; %outer diameter spiral pulley
            end
            obj = obj.compute_lengths(); %computing other lengths
            obj.springstroke = (fn*1e-3 - obj.L1); %[mm] possible stroke with spring
            obj.S = (obj.h_max - obj.L0 - obj.L1 - obj.h_adjust - obj.h_mech - obj.R2); %[m] possible stroke within the construction
            obj.height_check(); %checking if the height is correct
            obj.max_stroke = min(obj.S, obj.springstroke);
        end
        
        function obj = spring_properties_for_k(obj)
            %Computes system properties given the given stiffness of the
            %spring
            obj.fn = obj.Fn_tot / obj.k * 1e3; %maximal extension spring [mm]
            if obj.R1 > 0 %computing pulley size if there is a pulley
                obj.R2 = obj.F1/obj.k; %outer diameter spiral pulley
            end
            obj = obj.compute_lengths();
            obj.springstroke = (obj.fn*1e-3 - obj.L1); %[mm] possible stroke with spring
            obj.S = (obj.h_max - obj.L0 - obj.L1 - obj.h_adjust - obj.h_mech - obj.R2); %[m] possible stroke within the construction
            obj.height_check(); %checking if the height is correct
            obj.max_stroke = min(obj.S, obj.springstroke);
        end
        
        function [obj, isbetter] = higher_stroke(obj, l0, fn, Fn, n)
            %given new real spring parameters returns the spring system with the
            %larger maximal stroke
            isbetter = false;
            new_sys = obj.real_spring_properties(l0, fn, Fn, n);
            if new_sys.max_stroke > obj.max_stroke
                obj = new_sys;
                display('Better spring found')
                isbetter = true;
            end 
        end
        
        function [obj, isbetter] = lower_stiff(obj, l0, fn, Fn, n, desired_stroke) 
            %given new real spring parameters returns a spring system with lower
            %stiffness and that reaches the desired stroke
            isbetter = false;
            new_sys = obj.real_spring_properties(l0, fn, Fn, n);
            if new_sys.max_stroke <= 0
                %if the spring has negative stroke it does not work
                return 
            end
            if new_sys.k <= 0 
                %if the new spring system has negative stiffness is does
                %not work
                return
            end
            
            if new_sys.max_stroke < desired_stroke
               %the new max stroke is smaller than the desired stroke
               if new_sys.max_stroke < obj.max_stroke
                   %the new max stroke is smaller than the maximal stroke
                   %and the desired stroke, thus it is not better
                   return
               end
            end
            
            if new_sys.k < obj.k
                obj = new_sys;
                display('Better spring found')
                isbetter = true;
            end   
        end
    end
    
 
end

