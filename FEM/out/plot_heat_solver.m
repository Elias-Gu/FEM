clearvars; close all; clc;

numm = dir(strcat('heat_data/', "*.txt"));
step_number = length(numm) - 1;         % Number of steps

figure()
title_var = {'$\dot{u}(\mathbf{x},t) - \Delta u(\mathbf{x}) = 4\pi\sin(\pi x)\sin(\pi y), \mathbf{x} \in \Omega$',...
    '$\nabla u(\mathbf{x},t) \cdot \mathbf{n} = 8\pi \pmatrix{\sin(2\pi y) \cr \sin(2\pi y)}, \mathbf{x} \in \{\mathbf{x}|x=1\}\cup t \in [0,2]$', ...
    '$u(\mathbf{x},t)= 4\cos(2\pi t), \mathbf{x} \in \{ \mathbf{x}|x=0 \}$',...
    '$u(\mathbf{x},0) = 0$'};
output_file = 'heat_data/nodes_coo_0.txt';
data = dlmread(output_file);
dx = min(data(:,1)):1/sqrt(length(data(:,1))):max(data(:,1));
[X, Y] = meshgrid(dx, dx);

Z = griddata(data(:,1), data(:,2), data(:,3), X, Y);
surf(X, Y, Z);
xlabel('$x$', 'interpreter', 'latex');
ylabel('$y$', 'interpreter', 'latex');
zlabel('$u(x,y)$', 'interpreter', 'latex');
title(title_var, 'interpreter', 'latex')
zlim([-10 10]);
gif('heat_solution.gif', 'frame', gcf); %https://www.mathworks.com/matlabcentral/fileexchange/63239-gif

for step = 1:step_number
    
    output_file = ['heat_data/nodes_coo_', num2str(step), '.txt'];
    data = dlmread(output_file);

    Z = griddata(data(:,1), data(:,2), data(:,3), X, Y);
    surf(X, Y, Z);
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$y$', 'interpreter', 'latex');
    zlabel('$u(x,y)$', 'interpreter', 'latex');
    title(title_var, 'interpreter', 'latex')
    zlim([-10 10]);
    
    pause(0.01);
    gif

end




