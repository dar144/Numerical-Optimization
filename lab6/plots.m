clc; clear; close all;
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});

% ------------------- 5a ------------------- 

fileID = fopen('5a.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
xlabel('x');
ylabel('y');
title("nx = ny = " + 50 + " \epsilon_1 = \epsilon_2 = 1");

% ------------------- 5b ------------------- 

fileID = fopen('5b.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
xlabel('x');
ylabel('y');
title("nx = ny = " + 100 + " \epsilon_1 = \epsilon_2 = 1");

% ------------------- 5c ------------------- 

fileID = fopen('5c.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
xlabel('x');
ylabel('y');
title("nx = ny = " + 200 + " \epsilon_1 = \epsilon_2 = 1");

% ------------------- 6a ------------------- 

fileID = fopen('6a.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
clim([-0.8 0.8])
xlabel('x');
ylabel('y');
title("nx = ny = 100 \epsilon_1 = \epsilon_2 = 1");

% ------------------- 6b ------------------- 

fileID = fopen('6b.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
clim([-0.8 0.8])
xlabel('x');
ylabel('y');
title("nx = ny = 100 \epsilon_1 = 1 \epsilon_2 = 2");


% ------------------- 6c ------------------- 

fileID = fopen('6c.txt','r');
A = fscanf(fileID,'%f');
fclose(fileID);

A = reshape(A,3,[])';
len = A(:,1);
x = reshape(A(:,1),sqrt(length(len)),[]);
y = reshape(A(:,2),sqrt(length(len)),[]);
V = reshape(A(:,3),sqrt(length(len)),[]);

figure;
s = pcolor(x,y,V);
set(s, 'EdgeColor', 'none');
colormap(mycolormap);
colorbar;
clim([-0.8 0.8])
xlabel('x');
ylabel('y');
title("nx = ny = 100 \epsilon_1 = 1 \epsilon_2 = 10");









 


