clearvars; close all; clc;

output_file = 'poisson_data.txt';
data = dlmread(output_file);
dx = min(data(:,1)):1/sqrt(length(data(:,1))):max(data(:,1));
[X, Y] = meshgrid(dx, dx);

Z = griddata(data(:,1), data(:,2), data(:,3), X, Y);
surf(X, Y, Z);

xlabel('$x$', 'interpreter', 'latex');
ylabel('$y$', 'interpreter', 'latex');
zlabel('$u(x,y)$', 'interpreter', 'latex');
title({'$\Delta u(\mathbf{x}) = 4\pi\sin(2\pi x)\cos(2\pi y), \mathbf{x} \in \Omega$',...
    '$\nabla u(\mathbf{x}) \cdot \mathbf{n} =  \pmatrix{\frac{\pi}{5}\cos(2\pi x) \cr \sin(2\pi y)}, \mathbf{x} \in \{\mathbf{x}|x=1\cup\mathbf{x}|y=0,1\}$', ...
    '$u(\mathbf{x})= \cos(2\pi x) \cos(2\pi y), \mathbf{x} \in \{ \mathbf{x}|x=0 \}$'}, 'interpreter', 'latex')
