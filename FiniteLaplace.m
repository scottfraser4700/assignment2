%% Laplace's Equation - Part 1: One Dimensional Case 
% The first part of this problem is a solution to Laplace's equation along
% a 1-D section of a material with one side connected to a voltage source
% and the other end grounded. We expect to have a linear solution here. We
% are not using boundary conditions where the y-derivative of the potential
% is smooth at the upper boundary, so this is strictly a 1-D problem.

%initializing matrices
nx = 150;
ny = 100;
G = sparse(nx);
B = zeros(1,nx);
for i = 1:nx
    %Set up G-Matrix
     if i == 1
        G(i,i) = 1;
        B(i) = 1;
     elseif i == nx
        G(i,i) = 1;
     else
        G(i,i+1) = 1;
        G(i,i-1) = 1;
        G(i,i) = -2;
    end
end
%solve for V and plot
V = G\B';
X = linspace(0,150,length(V));
figure(1)
plot(X,V)
title('1D Voltage Plot')
xlabel('x dimension')
ylabel('potential')





