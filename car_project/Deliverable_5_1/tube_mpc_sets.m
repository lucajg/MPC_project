function  tube_mpc_sets(sys, Ts)
    % Get discrete dynamics
    [~, Ad, Bd, ~, ~] = Car.c2d_with_offset(sys, Ts);
    
    % Set Q and R
    Q = eye(2)*10;
    R = 1;

    % Use DLQR controller
    K = dlqr(Ad,-Bd,Q,R);
    K = -K;
    
    % Compute closed loop controller
    Acl = Ad - Bd*K;
    
    % Define W, the disturbance set in the state space
    u_tilde_min = -0.5 + sys.UserData.us;
    u_tilde_max =  0.5 + sys.UserData.us;
    F = [1
         -1];
    f = [u_tilde_max
         -u_tilde_min];
    W = Bd*Polyhedron(F,f);

    %% Compute the minimal robust invariant set E
    E = W;
    tol = 1e-2; % tolerance for convergence
    maxIter = 100;
    AW_i = W ;
    for i = 1:maxIter
        % Minkowski sum of Acl^i*W and E_prev
        AW_i = Acl*AW_i;
        E = AW_i + E;
        
        % Check for convergence to stop computation
        if norm(Acl^i) <= tol 
            break;
        end
    end
    
    %% Compute tightened constraints
    % Define state constraints
    x_safe_pos = 6.5;
    Fx = [ 0, 0 ; -1, 0 ; 0, 0 ; 0, 0];
    fx = [             0
          -(6 - x_safe_pos)
                       0
                       0   ];
    X = Polyhedron(Fx,fx);
    
    % Define input constraints
    umin = -1;
    umax = 1;
    Fu = [    1;    -1];
    fu = [ umax; -umin];
    U = Polyhedron(Fu,fu);
    
    % Tighten state and input constraints
    X_tilde = X-E;
    U_tilde = U-K*E;
    
    %% Maximal Positive invariant set (terminal set)
    % Combine tightened constraints for closed-loop system
    Fcl = [X_tilde.A; U_tilde.A*K];
    fcl = [X_tilde.b; U_tilde.b  ];
    
    % Compute maximal positive invariant set from tightened constraints
    Xf = polytope(Fcl,fcl);
    while 1
        prevXf = Xf;
        [T,t] = double(Xf);
        Xf = polytope([T;T*Acl],[t;t]);
        %Xf = intersect(Xf, presetXf);
        if isequal(prevXf, Xf)
            break
        end
    end
    [Ff,ff] = double(Xf);
    Xf = Polyhedron('A',Ff,'b',ff);
    
    %% Save useful variables
    save('tube_mpc_data.mat', 'Xf', 'x_safe_pos', 'K', 'U', 'U_tilde', ...
         'X', 'X_tilde', 'E', 'Q', 'R', 'W')
end