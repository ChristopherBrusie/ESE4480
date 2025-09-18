phi = pi/6;

d = 14;
l = 1;
theta1 = pi/8;
theta2 = pi/5;

base = [0 0 0 0]';


original = [1,0,0]; % object space


% link 2 (in its own frame) - forearm
p1 = [0,0,0];
p2 = [0,12,0];
link2_frame2 = [p1' p2'; ones(1,2)];

% link 1 (in its own frame) - bicep
pts = [0, 0, 0, 0, 0;
       0, 0, 6, 6, 8;
       0, -1, -1, 0, 0];
link1_frame1 = [pts; ones(1,5)];



for i = 0:0.1:2*pi
    for j = 0:0.1:2*pi
        % all points in world frame
        full_model_frame0 = [TRANS0_1(link1_frame1, i), TRANS0_1(TRANS1_2(link2_frame2, j, d), i)];
        
        
        % PLOT ROBOT
        pointsorg = full_model_frame0;
        s = 1;
        g = 7; % total num of points
        plot3(pointsorg(1,s:g), pointsorg(2,s:g), pointsorg(3,s:g), 'LineWidth',3);
        xlim([-20, 20]);
        ylim([-20, 20]);
        zlim([-20, 20]);
        hold on
        pause(0.001);
    end
end


%%
simulate_pendulum(pi/5, pi/8);

function points_frame1 = TRANS1_2(link2_frame2, theta2,d) % this is task 5

    t2_1 = [cos(theta2) -sin(theta2)       0       0;
              0                   0        -1       d;
            sin(theta2)  cos(theta2)        0       0;
                0                0            0       1];

    points_frame1 = t2_1*link2_frame2;



end

function points_frame0 = TRANS0_1(link1_frame1, theta1) % this is task 3

    t1_0 = [cos(theta1) -sin(theta1)       0       0;
     sin(theta1) cos(theta1)        0        0;
     0           0            1            0;
     0             0        0              1];

    points_frame0 = t1_0*link1_frame1;


end

function simulate_pendulum(theta1, theta2)

    d = 14;
    % link 2 (in its own frame) - forearm
    p1 = [0,0,0];
    p2 = [0,12,0];
    link2_frame2 = [p1' p2'; ones(1,2)];
    
    % link 1 (in its own frame) - bicep
    pts = [0, 0, 0, 0, 0;
           0, 0, 6, 6, 8;
           0, -1, -1, 0, 0];
    link1_frame1 = [pts; ones(1,5)];
    full_model_frame0 = [TRANS0_1(link1_frame1, theta1), TRANS0_1(TRANS1_2(link2_frame2, theta2, d), theta1)];


    % PLOT ROBOT
    pointsorg = full_model_frame0;
    s = 1;
    g = 7; % total num of points
    plot3(pointsorg(1,s:g), pointsorg(2,s:g), pointsorg(3,s:g), 'LineWidth',3);


end




