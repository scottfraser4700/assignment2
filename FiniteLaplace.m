%% Laplace's Equation - Part 1: One Dimensional Case 

%TODO Set up boundary conditions with conductivity maps
%for part 1 only doing a 1d solution from 0 to L
%for part 1&2 the conductivity map will use resistors of value 1
%for part 3 the boxes will have very low but finite conductivity

%TODO Set up G-Matrix similar to PA
%TODO implement and compare to analytical solution
%required plots: 2d plot of V(x), surface plot of V(x,y)

%TODO with altered conductivity map, plot conductivity, potential, electric
%field, and current density. play around with mesh density, bottleneck
%width, and conductivity and produce graphs

nx = 50;
ny = 50;
G = sparse(nx);
B = zeros(1,nx);
cMap = ones(nx*ny);
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
X = linspace(0,1,length(V));
figure(1)
plot(X,V)
title('1D Voltage Plot')


%setting conductivity map for 3rd part
% for q = 1:nx
%     for w = 1:ny
%         if (q<0.6*nx && q>0.4*nx && w>0.6*ny) || (q<0.6*nx && q>0.4*nx && w<0.4*ny)
%             cMap(q,w) = 1e-2;
%         end
%     end
% end
