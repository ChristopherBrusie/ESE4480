
function t01 = T01(t1)
    t01 = [cos(t1) -sin(t1) 0 0;
           sin(t1) cos(t1)  0 0;
           0       0        1 0;
           0       0        0 1];
end
function t12 = T12(t2)
    t12 = [cos(t2) -sin(t2) 0 0;
           0       0        -1 8;
           sin(t2) cos(t2)  1 0;
           0       0        0 1];
end
function trans1 = TRANS0_1(t1, link1)
    
    trans1 = zeros(size(link1,1), 4);
    tm1 = T01(t1);
    v1 = [link1,ones(5, 1)];
    for n = 1:size(link1, 1)
        trans1(n, :) = (tm1*(v1(n, :))')';
    end
end

function trans2 = TRANS1_2(t1, t2, link2)

    trans2 = zeros(size(link2, 1), 4);
    tm2 = T01(t1)*T12(t2);
    v2 = [link2, ones(2, 1)];
    for n = 1:size(link2, 1)
        trans2(n, :) = (tm2*(v2(n, :))')';
    end
end

function plot_robot(t1, t2)
    base = [0, 0, -10;
            0, 0, 0];
    link1 = [0,0,0; 
         0,0,-1; 
         0,6,-1; 
         0,6,0; 
         0,8,0];
    link2 = [0, 0, 0;
         0, 12, 0];
    trans_link1 = TRANS0_1(t1, link1);
    trans_link2 = TRANS1_2(t1, t2, link2);
    robot_points = [base; trans_link1(:, 1:3); trans_link2(:, 1:3)];
    plot3(robot_points(:, 1), robot_points(:, 2), robot_points(:,3))
    hold on
    axis equal
end
% Test
figure;
plot_robot(0, 0)
axis([-12 12 -12 12 -12 12])
% Work envelope
function plot_end_effector(t1, t2)
base = [0, 0, -10;
            0, 0, 0];
    link1 = [0,0,0; 
         0,0,-1; 
         0,6,-1; 
         0,6,0; 
         0,8,0];
    link2 = [0, 0, 0;
         0, 12, 0];
    trans_link1 = TRANS0_1(t1, link1);
    trans_link2 = TRANS1_2(t1, t2, link2);
    robot_points = [base; trans_link1(:, 1:3); trans_link2(:, 1:3)];
    scatter3(robot_points(end, 1), robot_points(end, 2), robot_points(end,3), "rx")
end


hold on
for i = 0:pi/8:2*pi
    for j = 0:pi/8:2*pi
        plot_end_effector(i, j)
    end
end
axis equal;
axis([-12 12 -12 12 -12 12])
hold off
