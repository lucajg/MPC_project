function  tube_mpc_sets(sys, Ts)
    
    close all 
    
    [~, Ad1, Bd2, ~, ~] = Car.c2d_with_offset(sys, Ts);
    Ad = [1 0.0999
         0 0.9976];
    Bd = [0.0083
        0.1665];

    Ad = Ad1;
    Bd = Bd2;
    
    K = dlqr(Ad,-Bd,eye(2),1);
    K = -K;
    %K = -place(Ad,-Bd,[0.5,0.4]);
    
    
    Acl = Ad - Bd*K;
    
    P = dlyap(Acl,eye(2));
    
    eig(Acl);
    
    u_tilde_min = -0.5 + sys.UserData.us;
    u_tilde_max =  0.5 + sys.UserData.us;
    
    F = [1
         -1];
    
    f = [u_tilde_max
         -u_tilde_min];
    
    W = Bd*Polyhedron(F,f); %Polyhedron(F*Bd',f); %
    
    % Acl = A + B*K is your closed-loop matrix
    % W is a Polyhedron representing the disturbance set
    E_prev = W;
    E = W;
    tol = 1e-2; % some tolerance for convergence
    maxIter = 100;
    AW_i = W ;
    figure; 
    W.plot('color','r','alpha',0.1)
    figure; 
    hold on;
    for i = 1:maxIter
        % Compute next set: Acl*E_prev ⊕ W
        AW_i = Acl*AW_i;
        E_next = AW_i + E_prev;
        
        % Check for convergence: E_next is a subset of E
        if norm(Acl^i) <= tol  % If Polyhedron E_next is contained in E
            E = E_next;
            break;
        else
            % Update: take intersection or form union to approach fixed point
            E_prev = E_next;
        end
        % Plot the initial polyhedron W in a separate figure
               % Create a new figure for W
        E_prev.plot('color','r','alpha',0.1);
        xlabel('x_1');
        ylabel('x_2');
        title('E_i');
    end
    hold off;
    
    
    
    % % Plot the final polyhedron E in another separate figure
    % figure;        % Create another new figure for E
    % W.plot('color','b','alpha',0.5);
    % xlabel('x_1');
    % ylabel('x_2');
    % title('Initial Polyhedron W');
    x_safe_pos = 15;
    Fx = [ 0 0
          -1 0
           0 0
           0 0];
    fx = [             0
          -(6 - x_safe_pos)
                       0
                       0];
    polyX = Polyhedron(Fx,fx);
    umin = -1;
    umax = 1;
    Fu = [1
          -1];
    fu = [ umax
          -umin];
    polyU = Polyhedron(Fu,fu);
    polyU_tilde = polyU-K*E;
    polyX_tilde = polyX-E;
    figure;
    polyU_tilde.plot('color','g','alpha',0.1);
    title('PolyU tilde');
    figure;
    polyU.plot('color','g','alpha',0.1);
    
    figure;
    hold on;
    polyX_tilde.plot('color','blue','alpha',0.1);
    polyX.plot('color','blue','alpha',0.1);
    hold off;
    
    KE = K*E;
    figure;
    KE.plot('color','blue','alpha',0.1);
    
    % Combine constraints for closed-loop system
    Ff = [Fx; Fu*K]
    ff = [fx; fu  ]
    
    % Plot
    figure;
    hold on;
    % Compute the invariant set:
    Xf = polytope(Ff,ff);
    % Xfplot = Polyhedron(double(Xf));
    plot(Xf,'r');
    while 1
        prevXf = Xf;
        [T,t] = double(Xf);
        Xf = polytope([T;T*Acl],[t;t]);
        %Xf = intersect(Xf, presetXf);
        if isequal(prevXf, Xf)
            break
        end
    end
    % Xfplot = Polyhedron(double(Xf));
    plot(Xf,'r');
    %Xmpi.plot('color','g','alpha',0.1);
    %[Ff,ff] = double(Xf);
    
    %ff-Ff*[3;0];
    
    
    
    title('Maximal invariant set X_f');
    xlabel('x');
    ylabel('V');
    grid on;
    hold off;
    
    save('tube_mpc_data.mat', 'Ff', 'ff', 'x_safe_pos', 'K', 'P', 'polyU_tilde', 'polyX_tilde', 'E')
end