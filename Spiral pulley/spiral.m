classdef spiral
    %Spiral: describes a spirral pulley attached to a spring and a mass
    %attached to a pulley with r2

    
    properties
        %MAIN PROPERTIES
        r1;                 %[mm] inner radius
        r2;                 %[mm] outer radius
        F0;                 %[N] initial load
        k;                  %[N/mm] spring stiffness
        name %name of the spiral type
        
        %SHAPE PROPERTIES
        r;                  %[mm] radius spiral at given theta
        theta;              %[rad] angle spiral at given r
        x;                  %[mm] x-coordinates spiral
        y;                  %[mm] y-coordinates spiral

        %Constant force spiral variables
        S                   %[] calculation constant
        %logarithmic spiral
        curve_angle         %[rad] defines slope logarithmic spiral
        
        %SPRING AND MASS relation
        x_m2            %[mm] max position mass
        x_m             %[mm] positions of the mass
        x_s             %[mm] position of spring
        rs              %[mm] radius as a function of angle phi
        phi             %[rad] input for rs
        spring_offset   %[m]offset spring has from equilibrium possition
        
        %FORCE BALANCE
        Fm          %[N] force on the mass from the spring
        Fs          %[N] force from the spring
        Fres        %[N] resultant force on the spring
        Fm_control  %[N] resultant force with pulley r1 to r2 (minimal expected behavior) 
        
        %POTENTIAL ENERGY
        Eg          %[J] gravitational potential energy 
        Es          %[J] spring potential energy
        Etot        %[J] total potential energy
    end
    
    
    methods
        function obj = spiral(r1,r2,F0,k)
            %constructing method
            obj.r1 = r1; 
            obj.r2 = r2;
            obj.F0 = F0;
            obj.k = k;
        end
        
        %%SHAPE FUNCTIONS
        function obj = Spiral_constant_Force(obj, n_steps)
            %Creates a constant force spiral
            %The final shape a n_step coordinates
            obj.name = 'Constant force spiral'; %naming object
            
            %setup 
            stepsize = (obj.r2- obj.r1)/n_steps; 
            obj.r = [obj.r2:-stepsize:obj.r1]; %computing range of r
            
            %Computing spiral
            obj.S = obj.F0*obj.r2/obj.k;  %computing variable S 
            obj.theta =  real( sqrt(obj.S^2 - obj.r.^4) )./ (2*(obj.r.^2)) + 0.5 *real(asin((obj.r.^2)/obj.S));
            
            %adding constant such that theta(1) = 0
            const = obj.theta(1); 
            obj.theta = obj.theta - const; 
            
            %computing cartesian coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = archimedes(obj, n_steps, thetamax)
            %Creates a archimedes spiral 
            %it has n-steps coordinates and goes until thetamax
            obj.name = 'Archimedes spiral';
            
            %setup
            stepsize = thetamax/n_steps; 
            obj.theta = [0:stepsize:thetamax];
           
            %creating shape
            a = (obj.r2 -obj.r1)/thetamax;  %slope variable
            obj.r = -a*obj.theta + obj.r2;
            
            %computing cartesian coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = logarithmic(obj, n_steps, thetamax, curve_angle)
            %creates a lorarithmic spiral
            %It has n_steps coordinates and rotates goes unitl thetamax
            %The slope of the system is proportional to curve_angle
            %IMPORTANT: 0 < curve_angle < pi/2
            obj.name = 'Logarithmic spiral';
            obj.curve_angle = curve_angle; %saving used curve_angle
            
            %setup            
            stepsize = thetamax/n_steps; 
            obj.theta = [0:stepsize:thetamax]; 
            
            %creating shape
            c = cot(obj.curve_angle);
            a = (obj.r2 -obj.r1) / exp(thetamax*c);
            obj.r = obj.r2 - a*exp(obj.theta*c);
            
            %computing cartesian coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
            
        end
        
        function obj = linear_extend(obj, extension_angle,n_steps)
            %Adds a constant radius shape for extension_angle radians of
            %the shape
            %the radius is the last used radius
            
            %setup
            th0 = max(obj.theta);
            step_size = extension_angle/n_steps;
            th_extend = [th0 : step_size : th0+extension_angle]; 
            radius = obj.r(length(obj.r));
            
            %shape creation
            obj.theta = [obj.theta, th_extend ]; 
            obj.r = [ obj.r, ones(1,length(th_extend))*radius]; 
            
            %converting to cartesian coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        function obj = close_shape(obj, n_steps)
            %Connects to two end points of the shape with a straight line
            %in cartensian coordinates
            %IMPORTANT: can only close the last 179 degrees of the spiral
            %shape
            
            %Finding coordinates
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
            obj.r = [obj.r,(r_fill(1:L))];
            obj.theta = [ obj.theta,(theta_fill(1:L)),];
            
            %recomputing the x and y coordinates
            obj.x = real(cos(obj.theta).*obj.r);
            obj.y = real(sin(obj.theta).*obj.r);
        end
        
        %%SYSTEM FUNCTIONS
        function obj = complete_analysis(obj, n_steps, spring_offset)
            %Performs all analysis of the system
            obj = obj.mass_spring_movement(n_steps, spring_offset);
            obj = obj.force_balance();
            obj = obj.potential_energy();
        end
        
        function obj = mass_spring_movement(obj,n_steps, spring_offset)
            %Computes the relation between x_m and x_s due to pulley
            obj.spring_offset = spring_offset; %safing the spring ofset
            
            %computing max extension of the spring
            obj.x_m2 = max(obj.theta)*obj.r2; 
            
            %setup
            stepsize = obj.x_m2/n_steps;
            obj.x_m = -[0:stepsize:obj.x_m2]; 
            obj.x_s = zeros(1,length(obj.x_m)); 
            obj.rs = zeros(1,length(obj.x_m)); 
            
            %movement calculation
            obj.phi = -obj.x_m / (obj.r2); % angles system goes trought
            obj.x_s(1) = obj.F0 / obj.k + spring_offset;  %initial extension spring
            
            %integration by parts
            obj.rs(1) = obj.r2;
            for i = 2 : length(obj.phi) %starting at second to be able to integrate
                [val, idx] = min(abs(obj.theta - obj.phi(i))); %finding the nearest index to the current angle
                obj.rs(i) = obj.r(idx); %finding the pulley radius for the given displacement of the mass
                
                %computing the next spring position with dx_m is equal to the stepsize
                obj.x_s(i) = obj.x_s(i-1) + obj.rs(i)/obj.r2* stepsize; 
            end
        end
        
        function obj = force_balance(obj)
            %computes the force balance resulting from the spiral pulley
            
            obj.Fs = (obj.x_s)*obj.k; %spring force versus mass displacement
            obj.Fm = obj.Fs.*(obj.rs/obj.r2); %force felt on the output
            obj.Fres = obj.Fm-obj.F0; %resulting force balance
            
            %Comparison result
            obj.Fm_control = -obj.x_m*obj.k* obj.r1/obj.r2 + obj.F0*(obj.r2/obj.r1) / obj.k - obj.F0;
            %reduced force spring : [force spring * reduction + initial extension spring - mass 
        end
        
        function obj = potential_energy(obj)
            obj.Eg = obj.F0*obj.x_m; %gravitational potential energy
            obj.Es = 0.5*obj.k*obj.x_s.^2; %spring potential energy intial position
            obj.Etot = obj.Eg + obj.Es; %total potential energy
        end
           
        %%PLOTTING FUNCTIONS
        function obj = plot_all(obj, n_steps, spring_offset)
           %Peforms complete analysis and plots all plotting functions
           obj = obj.mass_spring_movement(n_steps, spring_offset);
           obj.plot_r_theta();
           obj.plot_shape();
           obj = obj.plot_force_balance();
           obj = obj.plot_pot_energy();
        end
            
        function obj = plot_force_balance(obj)
            %computes the force balance and plots the normalised results
            obj = obj.force_balance(); %computing force balance
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
            %computes potential energy and plots results
            obj = obj.potential_energy(); %computing potential energy
            
            figure()
            hold on
            x_axis = obj.x_m;
            plot(x_axis, obj.Eg)
            plot(x_axis, obj.Es)
            plot(x_axis, obj.Etot)
            
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
            % * minimal theta is 0 and decreases with high values
            
            figure()
            plot(obj.theta,obj .r)
            title(append(obj.name, 'theta versus radius'))
            ylabel('radius [mm]'); xlabel('theta [rad]') 
        end
        
        function fig = plot_force_energy(obj)
            %plots the resultant force and total potential energy
            x_axis = obj.x_m; %choice of x_axis 
            
            fig = figure();
            hold on
            title(append(obj.name,' resulting force and total potential energy'))
            
            yyaxis left
            plot(x_axis,obj.Fres/obj.F0) 
            ylabel('Fres / F0 [N]') 
            ylim([-0.5 0.5]) %limits the plot shape
            
            yyaxis right
            plot(x_axis, obj.Etot)
            ylabel('Total Potential Energy [J]') 
            
            xlabel('Height Mass [mm]')
            xlim([-obj.r2*pi*1.5 0])
        end
        
        function fig = plot_shape(obj)
            %Plots the shape of the spiral in cartensian coordinates
            fig = figure();
            hold on
            plot(obj.x, obj.y);
            title(append(obj.name, ' pulley shape'))
            xlabel('x [mm]'); ylabel('y [mm]');
        end
        
        %OTHER FUNCTIONS
        function save(obj,name)
            %Saves coordinates into excel such that it can be imported
            output = table(transpose(obj.x),transpose(obj.y)); 
            writetable(output,name);
        end
   
    end
end

