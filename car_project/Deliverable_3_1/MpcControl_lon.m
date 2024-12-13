classdef MpcControl_lon < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   V_ref, u_ref - reference state/input
            %   d_est        - disturbance estimate
            %   x0other      - initial state of other car
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            V_ref = sdpvar(1);
            u_ref = sdpvar(1);

            % Disturbance estimate (Ignore this before Todo 4.1)
            d_est = sdpvar(1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this before Todo 5.1)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system.
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            obj = 0;
            con = [];

            % Decision variables over the horizon
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);

            % Weights for the cost function
            % Assume we only care about velocity tracking and input usage
            Q = diag([0, 1]);     % penalize deviation in velocity state only (x(2))
            R = 1;  %0.25           % penalize input deviation from reference
            
            [K, Qf, ~] = dlqr(mpc.A(2,2),mpc.B(2,1),Q(2,2),R);
            K = -K;

            % Linearization points
            xs = mpc.xs;
            us = mpc.us;

            % Input constraits
            umin = -1;
            umax =  1;
           
            %% Set up the MPC cost and constraints using the computed set-point
            
            % Initial condition
            con = con + (X(:,1) == x0);

            % Build the prediction model over the horizon
            for k = 1:N-1
                % Dynamics
                con = con + (X(:,k+1) == mpc.A*X(:,k) + mpc.B*U(:,k) + ...
                               (mpc.f_xs_us - mpc.A*xs - mpc.B*us)); 
                
                % Cost accumulation: 
                % Track velocity (X(2,k)) to V_ref and input U(k) to u_ref
                x_err = X(:,k) - [xs(1); V_ref];  % Error in state
                u_err = U(:,k) - u_ref;           % Error in input
                obj = obj + x_err'*Q*x_err + u_err'*R*u_err;
            end
            
            % P22 = dlyap(mpc.A(2,2),Q(2,2));
            % P = diag([0, P22]);
            % Terminal cost
            x_err_terminal = X(:,N) - [xs(1); V_ref];
            obj = obj + x_err_terminal'*diag([0,Qf])*x_err_terminal;

            % input constraints
            con = con + (umin <= U <= umax);
            % note: there are no constraints on the state for this
            % subsystem
            
            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = con + (u0 == U(:,1));

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {X, U};
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, V_ref, u_ref, d_est, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [Vs_ref, us_ref] = compute_steady_state_target(mpc, ref, d_est)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate (Ignore before Todo 4.1)
            % OUTPUTS
            %   Vs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state subsystem
            A = mpc.A(2, 2);
            B = mpc.B(2, 1);

            % Subsystem linearization steady-state
            xs = mpc.xs(2);
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % We want to find us_ref such that V = ref is an equilibrium:
            %
            % At steady state: Vs_ref = A(Vs_ref - xs) + B(us_ref - us) + xs
            %
            % Rearrange:
            % Vs_ref - xs = A(Vs_ref - xs) + B(us_ref - us)
            % (I - A)(Vs_ref - xs) = B(us_ref - us)
            % us_ref = us + ( (1 - A)/B )*( Vs_ref - xs )
            
            Vs_ref = ref;
            us_ref = us + ((1 - A)/B)*(Vs_ref - xs);
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
