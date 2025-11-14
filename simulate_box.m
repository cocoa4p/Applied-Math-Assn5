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
x0 = 0;
y0 = 0;
theta0 = deg2rad(35);
vx0 = 5;
vy0 =2;
vtheta0 = 35;
V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
% tspan = [0 100];
tspan = [0,30];
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

% Spring plotting example
[spring_plot_struct, P1, P2] = spring_plotting_example();

% Update spring plot
update_spring_plot(spring_plot_struct, P1, P2);

end