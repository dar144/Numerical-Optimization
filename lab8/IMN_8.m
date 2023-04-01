% clc; clear; close all;
format long;

nx = 400;
ny = 90;
i1 = 200;
i2 = 210;
j1 = 50;
delta = 0.01;
xa = 0.45;
ya = 0.45;

x = 0:delta:nx*delta;
y = 0:delta:ny*delta;

fileID = fopen('psi.dat','r');
f = fscanf(fileID,'%d %d %f', [3 Inf]);
f = f';
f = f(:,3);
psi = reshape(f,ny+1,[]);

[Vx,Vy] = speedField(psi, nx, ny, delta, i1, i2, j1);
drawSpeedField(x,y,Vx,Vy);


v = sqrt(Vx.^2+Vy.^2);
v_max = max(max(v));
delta_t = delta/(4*v_max);


[c1, t1, x_sr1, m1, m2, m3, m4, m5] = AD(nx, ny, 0, delta, Vx, Vy, xa, ya, x, y, 5000, i1, i2, j1, delta_t);
[c2, t2, x_sr2, n1, n2, n3, n4, n5] = AD(nx, ny, 0.1, delta, Vx, Vy, xa, ya, x, y, 5000, i1, i2, j1, delta_t);

figure;
plot(0:delta_t:delta_t*(length(c1)-1), c1, 0:delta_t:delta_t*(length(c2)-1), c2);
xlabel('t_n');
ylabel('c(t_n)');
title("c(t_n)");
legend('D = 0','D = 0.1');

figure;
plot(0:delta_t:delta_t*(length(c1)-1), x_sr1, 0:delta_t:delta_t*(length(c2)-1), x_sr2);
xlabel('t_n');
ylabel('x_{sr}(t_n)');
title("x_{sr}(t_n)");
legend('D = 0','D = 0.1');

% ____________ D = 0  _____________

figure;
h = 2; w = 3; i = 1;

subplot(h, w, i);
s = pcolor(x, y, m1);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 1000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, m2);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 2000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, m3);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 3000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, m4);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 4000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, m5);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 5000");
i = i + 1;


% ____________ D = 0.1  _____________

figure;
h = 2; w = 3; i = 1;

subplot(h, w, i);
s = pcolor(x, y, n1);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 1000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, n2);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 2000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, n3);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 3000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, n4);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 4000");
i = i + 1;

subplot(h, w, i);
s = pcolor(x, y, n5);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 5000");
i = i + 1;