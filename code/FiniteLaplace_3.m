%% Part 3 - Bottleneck Analysis
% In this part, a bottleneck of a different material is introduced into the
% medium. The bottleneck has a significantly lower conductivity than the
% surrounding medium, so its expected that less current should flow through
% it. Graphs for conductivity, electric potential, electric field, and
% current density are produced.
% By altering the values of nx and ny, we can play around with the mesh
% density. normally, the program is run with nx = 150, and ny = 100. By
% reducing both by a factor of 10, the total current flow increased by
% over 50% from 0.314A to 0.452A. Naturally when the mesh is made more
% dense, the current is lower.
% Narrowing the bottleneck will also have an effect on the total current
% output. I reduced the width of the opening by a factor of 2, and the
% current went from 0.314A to 0.220A. The current continues to drop off as
% the bottleneck is made thinner, as expected. This essentially is
% increasing the resistance of the material.
% As the conductivity of the bottleneck is changed, then the total current
% will change as well. As the conductivity approaches 1, it essentially
% behaves as if there is no bottleneck whatsoever since the surrounding
% medium has a conductivity of 1. It follows naturally that as the
% conductivity of the bottleneck approaches zero, the current density will
% approach a steady state value where it can't really interact with the
% bottleneck at all.
nx = 150;
ny = 100;
G = sparse(nx*ny, nx*ny);
B = zeros(1,nx*ny);
%setting conductivity map for 3rd part
cMap = ones(nx,ny);
 for q = 1:nx
     for w = 1:ny
         if ((q<(0.6*nx)&&q>(0.4*nx)&&w>(0.6*ny)) || (q<(0.6*nx)&&q>(0.4*nx)&&w<(0.4*ny)))
             cMap(q,w) = 1e-2;
         end
     end
 end

for i = 1:nx
    for j = 1:ny
        %Setting up the G-Matrix
        n = j + (i-1)*ny;
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,n) = 1;
            B(n) = 0;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G(n,n) = -(rxm + rxp + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        elseif j == ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            G(n,nxm) = rxm;
            G(n,nym) = rym;
            G(n,n) = -(rxm + rxp + rym + ryp);
        end
    end
end
%solving for matrix of potentials
V = G\B';
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end

%conductivity map plot
figure(4)
surf(cMap)
title('Conductivity')
ylabel('X dimension (L=150)')
xlabel('Y dimension (W=100)')
 
%plotting potential
figure(5)
surf(Vmap)
title('Potential')
ylabel('X dimension (L=150)')
xlabel('Y dimension (W=100)')

%plotting electric field
[Ex,Ey] = gradient(-Vmap);
figure(6)
quiver(Ex,Ey)
title('Electric field')
ylabel('X dimension (L=150)')
xlabel('Y dimension (W=100)')

%plotting current density
Jx = Ex .* cMap;
Jy = Ey .* cMap;
figure(7)
quiver(Jx,Jy)
title('Current Density')
ylabel('X dimension (L=150)')
xlabel('Y dimension (W=100)')

 %calculating current flow
 J = sqrt(Jx.*Jx+Jy.*Jy);
 C0 = sum(J(1,:));
 Cnx = sum(J(nx,:));
 Current = (C0 + Cnx)*0.5;
 fprintf('the total current flow is %f amps\n',Current);