clc; clear; close all;

nx = 40;
ny = 40;
N = (nx+1)*(ny+1);

% ____________ dt = 1 _____________

figure;
h = 2; w = 3; i = 1;

subplot(h, w, i);
fileID = fopen('7_100.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 100");
i = i + 1;

subplot(h, w, i);
fileID = fopen('7_200.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 200");
i = i + 1;

subplot(h, w, i);
fileID = fopen('7_500.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 500");
i = i + 1;


subplot(h, w, i);
fileID = fopen('7_1000.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 1000");
i = i + 1;

subplot(h, w, i);
fileID = fopen('7_2000.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 2000");
i = i + 1;


% ____________ dt = 0 _____________

figure;
h = 2; w = 3; i = 1;

subplot(h, w, i);
fileID = fopen('8_100.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 100");
i = i + 1;

subplot(h, w, i);
fileID = fopen('8_200.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 200");
i = i + 1;

subplot(h, w, i);
fileID = fopen('8_500.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');

colorbar;
xlabel('x');
ylabel('y');
title("it = 500");
i = i + 1;


subplot(h, w, i);
fileID = fopen('8_1000.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 1000");
i = i + 1;

subplot(h, w, i);
fileID = fopen('8_2000.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);
A = reshape(A,sqrt(N),[])';
s = pcolor(0:40,0:40,A);
set(s, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title("it = 2000");
i = i + 1;

