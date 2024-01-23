classdef spiral
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %MAIN PROPERTIES
        r1;         %inner radius
        r2;         %outer radius
        F0;         %initial load
        k;          %spring stiffness
        name %name of the spiral type
        
        %Constant Force variables
        S
        
        %SHAPE PROPERTIES
        r;          %input for theta
        theta;      %angle as a function of r
        x;
        y;
        spiral_max_stroke %maximal stroke allowed by the spiral
        curve_angle
        
        %SPRING AND MASS relation
        x_m2 %max extension mass
        x_m %positions of the mass
        x_s %position of spring
        rs %radius as a function of angle phi
        phi %input for rs
        spring_offset %offset spring has from equilibrium possition
        
        %FORCE BALANCE
        Fm %force on the mass from the spring
        Fs %force from the spring
        Fres %resultant force on the spring
        Fm_control %resultant force with pulley r1 to r2 (minimal expected behavior) 
        
        %POTENTIAL ENERGY
        Eg %gravitational potential energy 
        Es %spring potential energy
        Etot %total potential energy
    end
    
    
    methods
        function obj = spiral(r1,r2,F0,k)
            %initializing method
            obj.r1 = r1;
            obj.r2 = r2;
            obj.F0 = F0;
            obj.k = k;
        end
        
        %%SHAPE FUNCTIONS
        function obj = Spiral_constant_Force(obj, n_steps)
            %Returns theta for the given constant force 
            %size of the matrices is determined by n_steps
            obj.name = 'Constant force spiral';
            stepsize = (obj.r2- obj.r1)/n_steps;
            obj.r = [obj.r2:-stepsize:obj.r1]; %computing range of r
            obj.S = obj.F0*obj.r2/obj.k;  %computing variabble S 
            obj.theta =  real( sqrt(obj.S^2 - obj.r.^4) )./ (2*(obj.r.^2)) + 0.5 *real(asin((obj.r.^2)/obj.S));
            
            const = obj.theta(1); 
            obj.theta = (obj.theta); %flipping the shape
            obj.theta = obj.theta - const; %adding constant such that theta(1) = 0
            
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = archimedes(obj, n_steps, thetamax)
            %creates a archimedes spiral that takes n_steps
            obj.name = 'Archimedes spiral';
            stepsize = thetamax/n_steps; 
            obj.theta = [0:stepsize:thetamax]; %spiral will do one full rotation
           
            a = (obj.r2 -obj.r1)/thetamax;  %computing the slope variable such that only thetamax rotations are made
            
            obj.r = -a*obj.theta + obj.r2;
            obj.theta = (obj.theta); %flipping the shape
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = logarithmic(obj, n_steps, thetamax, curve_angle)
            obj.curve_angle = curve_angle;
            %creates a lorarithmic spiral
            % 0 < curve_angle < pi/2
            obj.name = 'Logarithmic spiral';
            stepsize = thetamax/n_steps; 
            obj.theta = [0:stepsize:thetamax]; %spiral will do one full rotation
            c = cot(obj.curve_angle);
            a = (obj.r2 -obj.r1) / exp(thetamax*c);
            obj.r = obj.r2 - a*exp(obj.theta*c);
            obj.theta = (obj.theta); %flipping the shape
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
            
        end
        
        function obj = linear_extend(obj, extension_angle,n_steps)
            %extends the pulley with constant radius
            th0 = max(obj.theta); %initial maximal angle
            step_size = extension_angle/n_steps;
            th_extend = [th0 : step_size : th0+extension_angle]; %e
            obj.theta = [flip(th_extend), obj.theta]; %extending the pulley
            
            obj.r = [ones(1,length(th_extend))*min(obj.r), obj.r]; %adding aditional radius
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = close_shape(obj, n_steps)
            %Can close the last 179 degrees of the spiral shape
            x_0 = obj.x(1);
            y_0 = obj.y(1);
            x_1 = obj.x(length(obj.x));
            y_1 = obj.y(length(obj.y));

            %computing in cartesian coordinates linear path
            x_fill = [x_1:(x_0 - x_1)/n_steps:x_0];
            y_fill = [y_1:(y_0 - y_1)/n_steps:y_0];
            
            %converting linear path in cylindrical coordinates
            r_fill = (x_fill.^2+y_fill.^2).^(1/2);
            theta_fill = 2*pi + atan(y_fill./x_fill); 
            %addjusting for loop over of atan
            theta_fill(theta_fill > 2*pi) = theta_fill(theta_fill > 2*pi) -pi; 
            
            %appending additional path
            L = length(theta_fill);
            obj.r = [obj.r,(r_fill(2:L))];
            obj.theta = [ obj.theta,(theta_fill(2:L)),];
            
            %recomputing the x and y coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        %%SYSTEM FUNCTIONS
        function obj = complete_analysis(obj, n_steps, spring_offset)
            %performs the compelte analysis
            obj = obj.mass_spring_movement(n_steps, spring_offset);
            obj = obj.force_balance();
            obj = obj.potential_energy();
        end
        
        function obj = mass_spring_movement(obj,n_steps, spring_offset)
            obj.spring_offset = spring_offset; %safing the spring ofset
            
            %computing max extension of the spring
            obj.x_m2 = max(obj.theta)*obj.r2; 
            
            stepsize = obj.x_m2/n_steps;
            obj.x_m = -[0:stepsize:obj.x_m2]; %computing linear position of the mass
            
            obj.x_s = zeros(1,length(obj.x_m)); %[m] creating empty spring position array
            obj.rs = zeros(1,length(obj.x_m)); %[rad] empty array for radius versus angle
            
            obj.phi = -obj.x_m / (obj.r2); % [rad] computing the angles that the system will go trough during the stroke of the system
            obj.x_s(1) = obj.F0 *(obj.r2/obj.r1)/ obj.k + spring_offset;  %initial extension spring
            obj.rs(1) = obj.r2;
            for i = 2 : length(obj.phi) %starting at the second position as the first angle position is zero
                validdata(i) = obj.phi(i) <= max(obj.theta); %if the current angle is greater than the maximal design angle the data is invalid
                [val, idx] = min(abs(obj.theta - obj.phi(i))); %finding the nearest index to the current angle
                obj.rs(i) = obj.r(idx); %finding the pulley radius for the given displacement of the mass
                obj.x_s(i) = obj.x_s(i-1) + obj.rs(i)/obj.r2* stepsize; %computing the next spring position with dx_m is equal to the stepsize
            end
            obj.spiral_max_stroke = max(-obj.x_m.*validdata); %maximal stroke
        end
        
        function obj = force_balance(obj)
            %computes the force balance resulting from the spiral pulley
            obj.Fs = (obj.x_s)*obj.k; %spring force versus mass displacement
            obj.Fm = obj.Fs.*(obj.rs/obj.r2); %force felt on the output
            obj.Fres = obj.Fm-obj.F0;
            %comparing to reduction r1/r2
            obj.Fm_control = -obj.x_m*obj.k* obj.r1/obj.r2 + obj.F0*(obj.r2/obj.r1) / obj.k - obj.F0; %reduced force spring : [force spring * reduction + initial extension spring - mass 
        end
        
        function obj = potential_energy(obj)
            obj.Eg = obj.F0*obj.x_m; %gravitational potential energy
            obj.Es = 0.5*obj.k*obj.x_s.^2; %spring potential energy intial position
            obj.Etot = obj.Eg + obj.Es;
        end
           
        %%PLOTTING FUNCTIONS
        function obj = plot_all(obj, n_steps, spring_offset)
           %computes everything and performs all plots
           obj = obj.mass_spring_movement(n_steps, spring_offset);
           obj.plot_r_theta();
           obj.plot_shape();
           obj = obj.plot_force_balance();
           obj = obj.plot_pot_energy();
        end
            
        function obj = plot_force_balance(obj)
            %computes the force balance and plots the normalised results
            obj = obj.force_balance();
            figure()
            hold on
            title(append(obj.name, ' normalised forces in the system'))
            x_axis = obj.x_m;
            plot(x_axis,obj.Fs/obj.F0) 
            plot(x_axis,obj.Fm/obj.F0)
            plot(x_axis,obj.Fres/obj.F0)
            plot(x_axis,obj.Fm_control/obj.F0) %reduced spring 
            %plot(x_axis,zeros(length(x_axis),1))
            legend('Spring force','output force','resulting force','Fres basic pulley system')
            xlabel('Displacement mass [mm]')
            ylabel('Force / F0 [ ]')
            hold off
         end
        
        function obj = plot_pot_energy(obj)
            obj = obj.potential_energy();
            
            figure()
            hold on
            x_axis = obj.x_m;
            plot(x_axis, obj.Eg)
            plot(x_axis, obj.Es)
            plot(x_axis, obj.Etot)
            %plot(x_axis,zeros(1,length(x_axis))) %plotting zero line
            title(append(obj.name,  ' energy versus position mass'))
            legend('gravity potential energy','spring potential energy','total potential energy')
            xlabel('Displacement mass [mm]')
            ylabel('Potential Energy [J]')
            hold off
        end
        
        function plot_r_theta(obj)
            %plots r versus theta
            %IMPORTANT the shape should be as follows
            % * at theta = 0 -> r is highest value
            %  * minimal theta is 0 and decreases with high values
            
            figure()
            plot(obj.theta,obj .r)
            title(append(obj.name, 'theta versus radius'))
            ylabel('radius [mm]'); xlabel('theta [rad]') %here we want to see the beginning be at theta == 0 
        end
        
        function fig = plot_force_energy(obj)
            x_axis = obj.x_m; %choice of x_axis 
            
            fig = figure();
            hold on
            title('Resulting force and Total potential energy') 
            plot(x_axis,obj.Fres/obj.F0)
            yyaxis left
            ylabel('Fres / F0 [N]') 
            
            yyaxis right
            plot(x_axis, obj.Etot)
            ylabel('Total Potential Energy [J]') 
            
            xlabel('Height Mass [mm]')
            %legend('Resulting force','Total potential energy')
        end
        
        function plot_shape(obj)
            figure()
            hold on
            plot(obj.x, obj.y);
            title(append(obj.name, ' pulley shape'))
            xlabel('x [mm]'); ylabel('y [mm]');
            hold off
            
        end

        
        
        
    end
end

