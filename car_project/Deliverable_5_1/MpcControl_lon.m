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

            obj = 0;
            con = [];

            load('tube_mpc_data.mat', 'Xf', 'x_safe_pos', 'K', 'U_tilde', 'X_tilde', 'E', 'Q', 'R')
            Ff = Xf.A;
            ff = Xf.b;
            x_safe = [x_safe_pos; 0];
            %DeltaX0 = sdpvar(nx,1);
            DeltaX0 = x0other - x0 - x_safe;

            DELTAZ = sdpvar(nx,N);

            U_T = sdpvar(1,N);

            xs = mpc.xs;
            us = mpc.us;

            umin = -1;
            umax =  1;

            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            
            %con = con + (DELTAZ(:,1) == DeltaX0);
            FXt = X_tilde.A;
            fXt = X_tilde.b;

            FUt = U_tilde.A;
            fUt = U_tilde.b;
            
           

            FE = E.A;
            fe = E.b;
            % con = con + (DELTAZ(:,1) == DeltaX0);
            % Minkowski sum constraint: DeltaX0 in DELTAZ(:,1) + E
            % ==> DeltaX0 - DELTAZ(:,1) in E
            con = con + (FE * (DeltaX0 - DELTAZ(:,1)) <= fe);
            for k = 1:N-1
                % Dynamics
                
                con = con + (DELTAZ(:,k+1) == mpc.A*DELTAZ(:,k) - mpc.B*U_T(:,k)); 
                
                con = con + (FXt*DELTAZ(:,k) <= fXt);
                con = con + (FUt*U_T(:,k) <= fUt);
                
                %con = con + (umin <= K*(DeltaX0 - DELTAZ(:,k)) + U_T(:,k) <= umax);
                
                % Cost accumulation: 
                
               % u_err = U_T(:,k) - u_ref;           % Error in input
                obj = obj + DELTAZ(:,k)'*Q*DELTAZ(:,k) + U_T(:,k)'*R*U_T(:,k);
            end
            Acl = mpc.A-mpc.B*K;
            P = dlyap(Acl,Q);
            obj = obj + DELTAZ(:,N)'*P*DELTAZ(:,N);

            con = con + (Ff*DELTAZ(:,N) <= ff);
            
            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            %con = con + (u0 == K*(DeltaX0 - DELTAZ(:,1)) + U_T(:,1));
            con = con + (u0 == K*(DeltaX0 - DELTAZ(:,1)) + U_T(:,1));

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {};
            
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
            Vs_ref = ref;
            us_ref = us;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
