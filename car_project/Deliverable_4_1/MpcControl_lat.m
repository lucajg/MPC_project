classdef MpcControl_lat < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this, not used)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE

            obj = 0;
            con = [];

            % Decision variables over the horizon
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);

            % Q = eye(2);
            % R = eye(1);
            % Weights for the cost function
            % Assume we primarily care about velocity tracking and input usage
            % Q = diag([0.01,1]) ;   
            % R = 0.01;             
            Q = diag([0.001385, 1]);
            R = 0.1;  

            % Linearization points
            xs = mpc.xs;
            us = mpc.us;
       
            % Input constraits

            ymin = -0.5;
            ymax =  3.5;

            thetamin = -deg2rad(5);
            thetamax =  deg2rad(5);

            deltamin =  -deg2rad(30);
            deltamax =   deg2rad(30);

            %% Compute maximal invariant set
          
            [K,Qf,~] = dlqr(mpc.A,mpc.B,Q,R);
            K = -K;
           
            A_cl = mpc.A + mpc.B*K;
            
            Fx = [ 1  0
                  -1  0
                   0  1
                   0 -1];
            fx = [ ymax
                  -ymin
                   thetamax
                  -thetamin];
            
            Fu = [ 1
                  -1];
            fu = [ deltamax
                  -deltamin];
            
            % Combine constraints for closed-loop system
            F_cl = [Fx; Fu*K];
            f_cl = [fx; fu  ];
            
           
            % Compute the invariant set:
            Xf = polytope(F_cl,f_cl);
            while 1
                prevXf = Xf;
                [T,t] = double(Xf);
                presetXf = polytope(T*A_cl,t);
                Xf = intersect(Xf, presetXf);
                if isequal(prevXf, Xf)
                    break
                end
            end
            [Ff,ff] = double(Xf);

            %ff-Ff*[3;0];
            
            % Plot
            figure(3);
            hold on;
            
            plot(polytope(Fx,fx),'g');
            plot(polytope(Ff,ff),'r');
            title('Maximal invariant set for lateral subsystem');
            scatter([0,3],[0,0]);
            xlabel('y');
            ylabel('\theta');
            grid on;
            hold off;
            %% Set up the MPC cost and constraints using the computed set-point
            
            % Initial condition
            con = con + (X(:,1) == x0 - xs);

            % Build the prediction model over the horizon
            for k = 1:N-1
                % Dynamics
                con = con + (X(:,k+1) == mpc.A*X(:,k) + mpc.B*U(:,k));
                con = con + (Fx*(X(:,k) + xs) <= fx);
                con = con + (Fu*(U(:,k) + us) <= fu);
                % Cost accumulation: 
                % Track position y X(1) to x_ref(1), 
                % angle theta X(2) to x_ref(2) 
                % and input U(k) to u_ref
                x_err = X(:,k) - (x_ref-xs);           % Error in state
                u_err = U(:,k) - (u_ref-us);           % Error in input
                obj = obj + x_err'*Q*x_err + u_err'*R*u_err;
            end
            
            % Terminal cost
            x_err_terminal = X(:,N) - (x_ref-xs);
            obj = obj + x_err_terminal'*Qf*x_err_terminal;
            
            % terminal constraint
            con = con + (Ff*(X(:,N) + xs - x_ref) <= ff);

            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = con + (u0 == U(:,1)+us); % set u0 to be the first input

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lat.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {X, U};
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, x_ref, u_ref, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [xs_ref, us_ref] = compute_steady_state_target(mpc, ref)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Steady-state system
            A = mpc.A;
            B = mpc.B;

            % Linearization steady-state
            xs = mpc.xs;
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % The y state has to be equal to the reference and theta must
            % be 0 at steady state
            C = [1 0];
            sizeC =  size(C);
            ny = sizeC(1);
            sizeB = size(B);
            nu = sizeB(2);
            nx = sizeB(1);
            aug_mat = [eye(size(A))-A,         - B
                                    C,  zeros(ny,nu)];
            r = [-A*xs - B*us + xs ; ref];
            xu = aug_mat\r;
            xs_ref = [xu(1); xu(2)];
            % we want the input delta to be 0 at steady state
            us_ref = xu(3);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end