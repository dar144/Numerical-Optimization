clc; clear; close all;
format long;

delta = 0.2;
nx = 128;
ny = 128;
x_max = delta*nx;
y_max = delta*ny;
TOL = 1e-8;

x = 0:delta:x_max;
y = 0:delta:y_max;
[X,Y] = meshgrid(x,y);

Vb1 = sin(pi*y/y_max);
Vb2 = -sin(2*pi*x/x_max);
Vb3 = Vb1;
Vb4 = -Vb2;

V = zeros(ny+1, nx+1);
V(:,1) = Vb1;
V(ny+1,:) = Vb2;
V(:,nx+1) = Vb3;
V(1,:) = Vb4;

k_max = 16;
k = 16;
it = 0;    
nextS = 0;
s_it = zeros(sqrt(k)+1, 1000);


for q = 1:sqrt(k_max)+1
    while 1
        
        kh = k/2;
        it = it + 1;
        
        for i = k+1:k:nx-k+1
            for j = k+1:k:nx-k+1
                V(j, i) = (V(j,i+k)+V(j,i-k)+V(j+k,i)+V(j-k,i))/4;
            end
        end
      
        currS = nextS;
        nextS = S(V, k, delta, nx, ny);
        s_it(q,it) = nextS;

        if abs((nextS - currS)/currS) < TOL
            break;
        end
  
    end
    
    if k > 1
        for i = 1:k:nx-k+1
            for j = 1:k:nx-k+1
                V(j+kh, i+kh) = (V(j,i)+V(j,i+k)+V(j+k,i)+V(j+k,i+k))/4;
                V(j+kh, i+k) = (V(j,i+k)+V(j+k,i+k))/2;
                V(j+k, i+kh) = (V(j+k,i)+V(j+k,i+k))/2;
                V(j, i+kh) = (V(j,i)+V(j,i+k))/2;
                V(j+kh, i) = (V(j,i)+V(j+k,i))/2;
            end
        end
    end


    V(:,1) = Vb1;
    V(ny+1,:) = Vb2;
    V(:,nx+1) = Vb3;
    V(1,:) = Vb4;


    if k > 1
        V_draw = V(2:end-1,2:end-1);
        V_draw = V_draw(any(V_draw),any(V_draw));
    else
        V_draw = V;
    end
    
    figure;
    s = pcolor( V_draw);
    set(s, 'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    xlabel('x');
    ylabel('y');
    title("k = " + k); 
    
    k = k/2;
end

s_it = s_it(:,any(s_it));
begin_plot = 1;

figure;
hold on
for q = 1:sqrt(k_max)+1
    y = nonzeros(s_it(q,:)');
    x = begin_plot:begin_plot+length(y)-1; 
    begin_plot = begin_plot+length(y)+1; 
    plot(x,y); 
    xlabel('it');
    ylabel('S');
    title("S(t), it = " + it);
    legend('k = 16','k = 8','k = 4','k = 2', 'k = 1');
end

