clc; close all; clear;
format long;

epsilon = 1;
delta = 0.1;
nx = 150;
ny = 100;
V1 = 10;
V2 = 0;
V = [V1.* ones(1,nx+1); zeros(ny, nx+1)];
TOL = 1e-8;
x_max = delta * nx;
y_max = delta * ny;
sigma_x = 0.1 * x_max;
sigma_y = 0.1 * y_max;
iter_max  = 1e5;

x = 0:delta:x_max;
y = 0:delta:y_max;
[X,Y] = meshgrid(x,y);

sigma_x2 = sigma_x^2;
sigma_y2 = sigma_y^2;

f1 = @(x,y) exp(-((x-0.35*x_max)^2/sigma_x2)-((y-0.5*y_max)^2/sigma_y2));
f2 = @(x,y) (-1)*exp(-((x-0.65*x_max)^2/sigma_x2)-((y-0.5*y_max)^2/sigma_y2));
f = @(x,y) f1(x,y)+f2(x,y);
ro = [];

for i = 1:nx+1
    for j = 1:ny+1
        ro(j,i) = f(x(i), y(j));
    end
end 

omega_g = [0.6, 1.0];
omega_l = [1.0, 1.4, 1.8, 1.9];


% omega_g = 0.6
[it, sumTab, it_final, error] = relaksacjaGlobalna(omega_g(1), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it1 = it;
sumTab1 = sumTab;
it_final1 = it_final;
% omega_g = 1.0
[it, sumTab, it_final, error] = relaksacjaGlobalna(omega_g(2), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it2 = it;
sumTab2 = sumTab;
it_final2 = it_final;


fc = contourf(X, Y, error);
% fc.LevelList = [0.003 0.0025 0.002 0.0015 0.001 0.0005 0];
% colorbar

figure;
semilogx(it1, sumTab1, it2, sumTab2,'LineWidth',3);
title("Relaksacja globalna");
xlabel('nr iteracji');
ylabel('S');
legend("\Omega = 0.6, "+it_final1+" it","\Omega = 1.0, "+it_final2+" it");  


% omega_l = 1.0
[it, sumTab, it_final] = relaksacjaLokalna(omega_l(1), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it1 = it;
sumTab1 = sumTab;
it_final1 = it_final;
% omega_l = 1.4
[it, sumTab, it_final] = relaksacjaLokalna(omega_l(2), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it2 = it;
sumTab2 = sumTab;
it_final2 = it_final;
% omega_l = 1.8
[it, sumTab, it_final] = relaksacjaLokalna(omega_l(3), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it3 = it;
sumTab3 = sumTab;
it_final3 = it_final;
% omega_l = 1.9
[it, sumTab, it_final] = relaksacjaLokalna(omega_l(4), V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL);
it4 = it;
sumTab4 = sumTab;
it_final4 = it_final;

figure;
semilogx(it1, sumTab1, it2, sumTab2, it3, sumTab3, it4, sumTab4,'LineWidth',3);
title("Relaksacja lokalna");
xlabel('nr iteracji');
ylabel('S');
legend("\Omega = 1.0, "+it_final1+" it","\Omega = 1.4, "+it_final2+" it",...
"\Omega = 1.8, "+it_final3+" it","\Omega = 1.9, "+it_final4+" it");


function [it, sumTab, it_final, error] = relaksacjaGlobalna(omega, V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL)    
    delta2 = delta^2;
    tmpIt = [];
    tmpSumTab = [];

    iteracja = 0;
    nextS = 0;
    currVs = V; 
    Vn = V;
    
    
    while 1
        iteracja = iteracja + 1;
        
        for i = 2:nx
            for j = 2:ny
                Vn(j,i) = (currVs(j,i+1)+currVs(j,i-1)+currVs(j+1,i)+currVs(j-1,i)+delta2/epsilon*ro(j,i))/4;
            end
        end
        
        Vn(:,1) = Vn(:,2);
        Vn(:,nx+1) = Vn(:,nx);
        
        nextVs = (1-omega)*currVs+omega*Vn;
        
        currS = nextS;
        nextS = S(nextVs, delta, x, y, ro, nx, ny);

        tmpIt(iteracja) = iteracja;
        tmpSumTab(iteracja) = nextS;


        currVs = nextVs;   
        
        if iteracja > iter_max || abs((nextS - currS)/currS) < TOL
            break;
        end
    
    end

    for i = 1:nx+1
    	for j = 1:ny+1
        	tmpErrTab(j,i) = delta2*currVs(j,i)+ro(j,i)/epsilon;
        end
    end   


    
    it_final = iteracja;
    it = tmpIt;
    sumTab = tmpSumTab;
    error = tmpErrTab;
    
end

function [it, sumTab, it_final] = relaksacjaLokalna(omega, V, delta, nx, ny, epsilon, iter_max, x, y, ro, TOL)    
    delta2 = delta^2;
    tmpIt = [];
    tmpSumTab = [];

    iteracja = 0;
    nextS = 0;
    Vs = V; 
    
    while 1
        iteracja = iteracja + 1;
        
        for i = 2:nx
            for j = 2:ny
                Vs(j,i) = (1-omega)*Vs(j,i)+omega/4*(Vs(j,i+1)+Vs(j,i-1)+Vs(j+1,i)+Vs(j-1,i)+delta2/epsilon*ro(j,i));
            end
        end
        
        Vs(:,1) = Vs(:,2);
        Vs(:,nx+1) = Vs(:,nx);
        
        currS = nextS;
        nextS = S(Vs, delta, x, y, ro, nx, ny);
        
        tmpIt(iteracja) = iteracja;
        tmpSumTab(iteracja) = nextS; 
        
        if iteracja > iter_max || abs((nextS - currS)/currS) < TOL
            break;
        end
    
    end   
    
    it_final = iteracja;
    it = tmpIt;
    sumTab = tmpSumTab;
    
end

function sum = S(V, delta, x, y, ro, nx, ny)
    delta2 = delta^2;
    tmpSum = 0;
    for i = 1:nx
        for j = 1:ny
            tmpSum = tmpSum + delta2*((((V(j,i+1)-V(j,i))/delta)^2)/2+(((V(j+1,i)-V(j,i))/delta)^2)/2 - ro(j,i)*V(j,i));
        end
    end
   
    sum = tmpSum;
end
