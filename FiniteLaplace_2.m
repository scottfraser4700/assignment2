%% Part Two: Two Dimensional Case
% Here we are solving a 2-dimensional version of Laplace's equation in what
% is essentially a rectangular pipe. Our goal is to compare our numerical
% solution to the analytical solution that can be obtained with a
% separation of variables and fourier analysis. Solving the differential
% equations on paper can be annoying, so a numerical solution is
% appropriate here. This approach would be significantly harder in a
% different coordinate system, but for rectangular coordinates, its OK.

% There is an inherent error involved with both methods. The analytical
% solution involves an infinite sum, but since that's not physically
% realizable, I cut it off at 100. On my plot of the analytical solution,
% the top edges are still wavy, as its comprised of a large quantity of
% sinusoids. I chose to end it at 100 because when watching it as a movie,
% it would stop visually changing around 100 but it would slow down
% dramatically afterwards.

% The finite difference solution I obtained looks really good, but with
% a wider mesh it would look much less smooth. 

nx = 150;
ny = 100;
G = sparse(nx*ny, nx*ny);
B = zeros(1,nx*ny);
cMap = ones(nx,ny);

for i = 1:nx
    for j = 1:ny
        %Setting up the G-Matrix
        n = j + (i-1)*ny;
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,n) = 1;
            B(n) = 1;
        elseif j == 1
            G(n,n) = 1;
        elseif j == ny
            G(n,n) = 1;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
            G(n,nxm) = 1;
            G(n,nym) = 1;
            G(n,n) = -4;
        end
    end
end
%solving for matrix of potentials
V = G\B';
%remapping back to i,j
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end
%plotting V(x,y)
 figure(2)
 surf(Vmap)
 title('Surface plot of potential (Matrix solution)')
 ylabel('X dimension (L=150)')
 xlabel('Y dimension (W=100)')

%analytical solution
[x, y] = meshgrid(-75:2:75, 0:2:100);
a = 100;
b = 75;
V_an = 0;

for N = 1:100
    if rem(N,2) == 1
        V_an = V_an + (4/pi)*(cosh(N*pi*x/a).*sin(N*pi*y/a)) ./ (N*cosh(N*pi*b/a));
        figure(3)
        surf(y,x,V_an)
        title('Surface plot of potential (analytical solution)')
        hold on
        pause(0.01)
    end
end


