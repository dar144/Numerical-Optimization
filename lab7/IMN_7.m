% clear; clc; close all;

delta = 0.01;
ro = 1;
mu = 1;
nx = 200;
ny = 90;
i1 = 50;
j1 = 55;
IT_MAX = 20000;

x_max = nx * delta;
y_max = ny * delta;

x = 0:delta:x_max;
y = 0:delta:y_max;

Q_wy = @(Q_we,y) Q_we*(y(ny+1)^3-y(j1+1)^3-3*y(j1+1)*y(ny+1)^2+...
                 3*y(j1+1)^2*y(ny+1))/(y(ny+1)^3);

q = [-1000 -4000 4000];

for k = 1:length(q)

    q_for = q(k);

    psi = zeros(ny+1,nx+1);
    zeta = zeros(ny+1,nx+1);

    [psi,zeta,u,v] = NS(nx, ny, Q_wy, q_for, psi, y, zeta, j1, i1, IT_MAX, delta, ro, mu);
    
    psi(1:j1,1) = NaN;
    psi(1,1:i1) = NaN;

    u(isnan(u))=0;
    v(isnan(v))=0;

    figure;
    contour(x, y, psi, 100);
    colorbar;
    xlabel('x');
    ylabel('y');
    title("\Omega = " + q_for +", \psi(x,y)");

    if(k ~= 3)
        figure;
        contour(x, y, zeta, 100);
        colorbar;
        xlabel('x');
        ylabel('y');
        title("\Omega = " + q_for + ", \zeta(x,y)");

        figure;
        s = pcolor(x, y, u);
        set(s, 'EdgeColor', 'none');
        colorbar;
        colormap(jet);
        xlabel('x');
        ylabel('y');
        title("\Omega = " + q_for + ", u(x,y)");
        
        figure;
        s = pcolor(x, y, v);
        set(s, 'EdgeColor', 'none');
        colorbar;
        colormap(jet);
        xlabel('x');
        ylabel('y');
        title("\Omega = " + q_for + ", v(x,y)");
    end

end











 

 





