function assignment_5_main
    % Parameters
    k = 5; % spring stiffness
    l0 = 5; % initial length of spring (when nothing is attached)
    PA = 1; % position of one end of spring
    PB = 2; % position of other end of spring
    x = 0; % x position of centroid of box
    y = 0; % y position of centroid of box
    theta = 30; % orientation of box
    Plist_box = [-1, 1; 1, 1; 1, -1; -1, -1]; % Points in the box frame

    % Compute force the spring exerts at one of its ends
    F = compute_spring_force(k,l0,PA,PB);
    disp(F);
    % Map box-frame coordinates to corresponding world-frame coordinates
    Plist_world = compute_rbt(x,y,theta,Plist_box);
    disp(Plist_world);
end