function plots()
    load('tube_mpc_data.mat', 'Xf', 'x_safe_pos', 'K', 'U', 'U_tilde', 'X', ...
         'X_tilde', 'E', 'Q', 'R', 'W')
    close all 
    
    %% Plot the polyhedron E
    figure;        
    hold on;
    E.plot('color','b','alpha',0.5);
    xlabel('$x$','Interpreter', 'latex');
    ylabel('$V$','Interpreter', 'latex');
    title('Minimal robust positive invariant set $\mathcal{E}$','Interpreter', 'latex');
    plot(0, 0, 'blackx', 'MarkerSize', 7, 'LineWidth', 1); 
    axis square;
    hold off;
    
    %% Plot Maximal invaritant set Xf
    
    figure;
    hold on;
    Xf.plot('color','r','alpha',0.3);
    title('Maximal robust positive invariant set $\mathcal{X}_f$','Interpreter', 'latex');
    xlabel('$x$','Interpreter', 'latex');
    ylabel('$V$','Interpreter', 'latex');
    grid on;
    plot(0, 0, 'blackx', 'MarkerSize', 7, 'LineWidth', 1); 
    axis square;
    hold off;
end