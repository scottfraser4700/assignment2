%% Part Two: Two Dimensional Case
nx = 50;
ny = 50;
G = sparse(nx*ny, nx*ny);
B = zeros(1,nx*ny);
cMap = ones(nx*ny);
for i = 1:nx
    for j = 1:ny
        %Setting up the G-Matrix
        n = j + (i-1)*ny;
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,n) = 1;
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

%analytical solution
[x, y] = meshgrid(0:0.02:1 , 0:0.02:1);
a = 1;
b = 1;
N = 1;
V_an = (4/pi).*(cosh(N.*pi.*x./a).*sin(N.*pi.*y./a)) ./ (N.*cosh(N.*pi.*b./a));

surf(x,y,V_an)
