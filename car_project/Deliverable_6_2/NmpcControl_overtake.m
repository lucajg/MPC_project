classdef NmpcControl_overtake < handle

    properties
        % The NMPC problem
        opti

        % Problem parameters
        x0, ref, x0other

        % Most recent problem solution
        sol

        % The input that you want to apply to the system
        u0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add any variables you would like to read to debug here
        % and then store them in the NmpcControl function below.
        % e.g., you could place X here and then add obj.X = X
        % in the NmpcControl function below.
        % 
        % After solving the problem, you can then read these variables 
        % to debug via
        %   nmpc.sol.value(nmpc.X)
        % 
        X, U
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    methods
        function obj = NmpcControl_overtake(car, H)

            import casadi.*
            


            N_segs = ceil(H/car.Ts); % Horizon steps
            N = N_segs + 1;          % Last index in 1-based Matlab indexing

            f_discrete = @(x,u) RK4(x,u,car.Ts,(@ car.f));

            nx = 4;
            nu = 2;

            % Define the NMPC optimization problem
            opti = casadi.Opti();
            
            % Parameters (symbolic)
            obj.x0 = opti.parameter(nx, 1);       % initial state
            obj.ref = opti.parameter(2, 1);       % target y, velocity
            obj.x0other = opti.parameter(nx, 1);  % initial state of other car

            % SET THIS VALUE TO BE YOUR CONTROL INPUT
            obj.u0 = opti.variable(nu, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Define your problem using the opti object created above


            X = opti.variable(nx,N+1); % state trajectory variables
            X_other = opti.variable(nx, N+1);
           

            U = opti.variable(nu,N);   % control trajectory (throttle, brake)
            A = [1 0 0 car.Ts
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1];

            H = [0.005 0
                 0 0.11];

            pos_x0 = obj.x0(1);
            pos_y0 = obj.x0(2);
            theta_0 = obj.x0(3) ;  
            v_0 = obj.x0(4); 

            pos_x   = X(1,:);
            pos_y = X(2,:);
            theta = X(3,:);
            speed = X(4,:);

            theta_max =  deg2rad(5);
            theta_min = -theta_max;

            throttle = U(2,:);
            steering = U(1,:);

            opti.subject_to(pos_x(1)==pos_x0);   % use initial position
            opti.subject_to(pos_y(1)==pos_y0);   % use initial position
            opti.subject_to(theta(1)==theta_0);   % use initial position
            opti.subject_to(speed(1)==v_0); % use initial speed

            V_err = speed - obj.ref(2)*ones(1,N+1);
            y_err = pos_y - obj.ref(1)*ones(1,N+1);
            
            cost = 10*(V_err*V_err')  + ... 
                1*(y_err*y_err') + ... 
                10 * (throttle*throttle')  + ...  
                10 * (theta)*(theta)'  + ...  
                10 * (steering*steering'); 
           
            opti.subject_to(X_other(:,1)==obj.x0other);
            
            p = [X(1,:); X(2,:)];
            pL = [X_other(1,:); X_other(2,:)];

            for k=1:N % loop over control intervals
                opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
                opti.subject_to(X_other(:,k+1) == A*X_other(:,k));
                opti.subject_to(((p(:,k)-pL(:,k))'*H*(p(:,k)-pL(:,k)))>=1);
            end
            
            opti.subject_to(-1 <= throttle <= 1); 
            opti.subject_to(-0.5 <= pos_y <= 3.5); 
            opti.subject_to(-0.0872 <= theta <= 0.0872);
            opti.subject_to(-0.5236 <= steering <= 0.5236);  

            opti.subject_to( obj.u0 == U(:,1) );

            opti.minimize(cost);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Store the defined problem to solve in get_u
            obj.opti = opti;

            % Setup solver
            options = struct;
            options.ipopt.print_level = 0;
            options.print_time = 0;
            options.expand = true;
            %options.ipopt.max_iter = 10000000000;
            obj.opti.solver('ipopt', options);
        end

        function u = get_u(obj, x0, ref, x0other)

            if nargin < 4
                x0other = zeros(4, 1);
            end

            % Compute solution from x0
            obj.solve(x0(1:4), ref, x0other(1:4));
            
            u = obj.sol.value(obj.u0);
        end

        function solve(obj, x0, ref, x0other)

            % Pass parameter values
            obj.opti.set_value(obj.x0, x0);
            obj.opti.set_value(obj.ref, ref);
            obj.opti.set_value(obj.x0other, x0other);
            obj.sol = obj.opti.solve();   % actual solve
            
            % Set warm start for next solve
            obj.opti.set_initial(obj.sol.value_variables());
            obj.opti.set_initial(obj.opti.lam_g, obj.sol.value(obj.opti.lam_g));
        end
    end
end
