function simulate_box()
% %define system parameters (these can be changed)
% box_params = struct();
% box_params.m = 1;                 
% box_params.I = 1/12*(10^2 + 1^2);  
% box_params.g = 9.81;    
% % spring stiffnesses
% box_params.k_list = [10 10 20 40]; 
% box_params.l0_list = [1.5 1.5 1.5 1.5];
% % spring endpoints
% box_params.P_world = [ ...
%     -2  2  2 -2;
%     -2 -2  2  2];
% % box-mounted points
% box_params.P_box = [ ...           
%     -1  1  1 -1;
%     -1 -1  1  1];

    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;
    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];

    %define system parameters
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;
%load the system parameters into the rate function
%via an anonymous function

my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

% test equlibirum

temp_function = @(V_in) my_rate_func(0,V_in);

solver_params = struct();
solver_params.dxtol = 1e-10;
solver_params.ftol = 1e-10;
solver_params.max_iter = 200;
solver_params.dxmax = 1e8;
solver_params.numerical_diff = 1;

x_guess = [.1 -2 pi/6 0 0 0];

[x_root] = multi_newton_solver_5(temp_function,x_guess,solver_params);

disp('x_root = ');
disp(x_root);

V_eq = [x_root(1,1); x_root(2,1); x_root(3,1); 0; 0; 0];

x0 = 2;%.0182; %0;
y0 = .1478; %0;
theta0 = -0.0321;%.2;
vx0 = 0;%.5;
vy0 =0; %.2;
vtheta0 = 0;%.35;

%small number used to scale initial perturbation
epsilon = .4;%your code here
V0 = V_eq + epsilon*[x0;y0;theta0;vx0;vy0;vtheta0];

% V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
tspan = [0,5];
% t_range = linspace(tspan(1),tspan(2),100);

% "Original" fourth-order Runge-Kutta
ogRunge = struct();
ogRunge.C = [0; 1/2; 1/2; 1];
ogRunge.B = [1/6, 1/3, 1/3, 1/6];
ogRunge.A = [0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];

expMethod = ogRunge;

h_ref = 0.05;
% [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration_global(my_rate, tspan, V0, h_ref, expMethod);
% Run the integration
[tlist,Vlist] = explicit_RK_fixed_step_integration_global(my_rate_func, tspan, V0, h_ref, expMethod);
[tlist_test,Vlist_test] = explicit_RK_fixed_step_integration_global(my_rate_func, tspan, V_eq, h_ref, expMethod);


% LINEARIZATION ---------
f = @(V) box_rate_func(0, V, box_params);
J_approx = approximate_jacobian(f, V_eq);

% linear rate function
my_linear_rate = @(t_in, V_in) J_approx * (V_in - V_eq);

% linearized simulation
[t_linear, V_linear] = explicit_RK_fixed_step_integration_global(my_linear_rate, tspan, V0, h_ref, expMethod);

% linear vs nonlinear
%figure(1); hold on;
%plot(tlist,   Vlist(:,1), 'b', 'LineWidth', 1.5);
%plot(t_linear,V_linear(:,1), 'r--', 'LineWidth', 1.5);
%legend('Nonlinear', 'Linearized');
%title('x(t): Nonlinear vs Linear');
%xlabel('time'); ylabel('x'); hold off;

% Spring plotting example
% [spring_plot_struct, P1, P2] = spring_plotting_example();

% num_zigs = 5;
% w = .1;
% hold on;
% spring_plot_struct = initialize_spring_plot(num_zigs,w);

% Update spring plot
% update_spring_plot(spring_plot_struct, P1, P2);

% num_frames = 5000;
% figure(2);
% % animate_box(tlist,Vlist,box_params,num_frames);
% figure(3);
% animate_box(tlist_test,Vlist_test,box_params,num_frames);


% MODAL ANALYSIS
%[Jacobian_box] = approximate_jacobian(my_rate_func, V0);

Q = -[vx0;vy0;vtheta0]/[x0;y0;theta0];
[U_modes, D] = eig(Q);
omega_n = sqrt(abs(diag(D))); % Extract natural frequencies



%run the integration of nonlinear system
for i = 1:3
    U_mode = U_modes(:, i);

    %small number
    V0 = V_eq + epsilon*[U_mode;0;0;0];
    tspan = [0,5];

    [t_list, V_list, ~, ~] = explicit_RK_fixed_step_integration_global(my_rate_func, tspan, V0, h_ref, ogRunge);  
        
    % Predicted behavior from the modal decomposition
    x_modal = zeros(size(t_list));
    y_modal = zeros(size(t_list));
    theta_modal = zeros(size(t_list));
    t_list
    for j = 1: t_list(end)
        x_modal(j) = V_eq(1)+epsilon*U_mode(1)*cos(omega_n(i)*t_list(j));
        y_modal(j) = V_eq(2)+epsilon*U_mode(2)*cos(omega_n(i)*t_list(j));
        theta_modal(j) = V_eq(3)+epsilon*U_mode(3)*cos(omega_n(i)*t_list(j));
    end

    % Store modal data
    Vlist_modal = [x_modal, y_modal, theta_modal];

    % Plots comparing the predicted vibration mode to the simulated nonlinear
    % behavior
    figure(4);
    plot(t_list,V_list(:, 1), 'r-'); hold on
    plot(t_list, x_modal, 'b-');
    title("Predicted vibration mode vs. simulated nonlinear behavior");
    legend("Predicted vibration mode", "Simulated nonlinear behavior");

    figure(5);
    plot(t_list,V_list(:, 1), 'r-'); hold on
    plot(t_list, y_modal, 'b-');
    title("Predicted vibration mode vs. simulated nonlinear behavior");
    legend("Predicted vibration mode", "Simulated nonlinear behavior");
    
    figure(6);
    plot(t_list,V_list(:, 1), 'r-'); hold on
    plot(t_list, theta_modal, 'b-');
    title("Predicted vibration mode vs. simulated nonlinear behavior");
    legend("Predicted vibration mode", "Simulated nonlinear behavior");
end


end

function animate_box(tlist,Vlist,box_params,num_frames)
    T_total = tlist(end)-tlist(1);
    dt = T_total/num_frames;

    t = tlist(1);
    current_frame = 1;

    figure(1);
    hold on
    l = 20;
    axis equal
    axis([-l,l,-l,l]);
    
    
    box_plot = plot(0,0,'k','linewidth',1.5);
    spring_list = {};

    num_springs = size(box_params.P_world,2);

    num_zigs = 5;
    w = .1;
    
    
    for n = 1:num_springs
        spring_list{n} = initialize_spring_plot(num_zigs,w);
    end

    while t<=tlist(end)
        while tlist(current_frame)<t && current_frame < length(tlist)
            current_frame = current_frame+1;
        end
    
        V_current = Vlist(current_frame,:)';


        x_current = V_current(1);
        y_current = V_current(2);
        theta_current = V_current(3);

        Plist_boundary = compute_rbt(x_current,y_current,theta_current,box_params.boundary_pts);
        P1_list = box_params.P_world;
        P2_list = compute_rbt(x_current,y_current,theta_current,box_params.P_box);

        for n = 1:num_springs
            P1 = P1_list(:,n);
            P2 = P2_list(:,n);
            update_spring_plot(spring_list{n},P1,P2);
        end

        set(box_plot,'xdata',Plist_boundary(1,:),'ydata',Plist_boundary(2,:));
        drawnow;



        t = t+dt;
    end
    


end