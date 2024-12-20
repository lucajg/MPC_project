classdef NmpcControl < handle

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
        function obj = NmpcControl(car, H)

            import casadi.*

            N_segs = ceil(H/car.Ts); % Horizon steps
            N = N_segs + 1;          % Last index in 1-based Matlab indexing

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
            % t = 0:car.Ts:H; % Sample times
            % ode.name = 'ODE45';
            % ode.f_discrete = @(X,U) ode45(@(t,x) car.f(x,U),[0 h], X');
            % ode.X = obj.x0;
            
            % Define your problem using the opti object created above
            X = opti.variable(nx,N+1);
            U = opti.variable(nu,N  );

            V_err = X(4,:) - obj.ref(2)*ones(1,N+1);
            y_err = X(2,:) - obj.ref(1)*ones(1,N+1);
            theta = X(3, :);
            delta = U(1,:);
            cost =  1 * (V_err*V_err') + ...
                    1 * (y_err*y_err') + ...
                    0.1 * (theta * theta') + ...
                    0.1* (delta * delta')+ ...
                    -10* (U(2,:)*U(2,:)');
            
            opti.subject_to(X(:,1) == obj.x0);

            for k = 1:N
                opti.subject_to(X(:,k+1) == RK4(X(:,k),U(:,k),car.Ts,(@car.f)));
            end

            opti.subject_to(-1 <= U(2,:) <= 1);
            opti.subject_to(-0.5 <= X(2,:) <= 3.5);
            opti.subject_to( -0.0873 <= X(3,:) <= 0.0873);
            opti.subject_to(-0.5236 <= U(1,:) <= 0.5236);

            % change this line accordingly
            opti.subject_to( obj.u0 == U(:,1));

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
