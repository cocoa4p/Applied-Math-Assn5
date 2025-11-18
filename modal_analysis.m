function modal_analysis()
    [Jacobian_box, num_evals] = approximate_jacobian(my_rate_func, V0);

    U_mode = eig(Q); %your code here
    omega_n = %your code here
    %small number
    epsilon = 0.1; %your code here
    V0 = Veq + epsilon*[Umode;0;0;0];
    V0_1 = Veq + epsilon*[dx_p;dy_p;dtheta_p;vx_p;vy_p;vtheta_p];
    tspan = %your code here
    %run the integration of nonlinear system
    % [tlist_nonlinear,Vlist_nonlinear] =...
    % your_integrator(my_rate_func,tspan,V0,...);

    % Predicted behavior from the modal decomposition
    x_modal = Veq(1)+epsilon*Umode(1)*cos(omega_n*tlist);
    y_modal = Veq(2)+epsilon*Umode(2)*cos(omega_n*tlist);
    theta_modal = Veq(3)+epsilon*Umode(3)*cos(omega_n*tlist);

    % Plots comparing the predicted vibration mode to the simulated nonlinear
    % behavior
end