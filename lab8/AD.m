function [c_out, t_out, x_sr_out, m1, m2, m3, m4, m5] = AD(nx, ny, D, delta, Vx, Vy, xa, ya, x, y, it_max, i1, i2, j1, delta_t)
    sigma = 10 * delta;
    
    c1 = 1/(2*pi*sigma^2);
    c2 = 2*sigma^2;

    u0 = zeros(ny+1,nx+1);
    u1 = zeros(ny+1,nx+1);
    c = zeros(1, it_max);
    x_sr = zeros(1, it_max);
    
    for j = 1:ny
    	u0(j,:) = c1*exp(-((x-xa).^2+(y(j)-ya)^2)/c2);
    end
    
    c3 = 1/(1+(2*D*delta_t)/delta^2);

    for it = 1:it_max
        it
        u1 = u0;
        for k = 1:20
            for i = 1:nx+1
                for j = 2:ny
                    if i >= i1+1 && i <= i2+1 && j <= j1+1 
                        continue;
                    else
                        if i == 1 || i == nx+1
                            u1(j,i) = relaxation_edge(D, delta_t, delta, j, i, u0, u1, Vx, Vy, c3, nx);  
                        else
                            u1(j,i) = relaxation(D, delta_t, delta, j, i, u0, u1, Vx, Vy, c3);  
                        end
                    end
                end
            end
        end
        u0 = u1;
        
        c_tmp = 0;
        x_sr_tmp = 0;
        for i = 1:nx+1
            for j = 1:ny+1
                c_tmp = c_tmp + u0(j,i)*delta^2;
                x_sr_tmp = x_sr_tmp + x(i)*u0(j,i)*delta^2;
            end
        end

        if(it == 1000)
            m1 = u0;
        end
        if(it == 2000)
            m2 = u0;  
        end
        if(it == 3000)
            m3 = u0;  
        end
        if(it == 4000)
            m4 = u0;
        end
        if(it == 5000)
            m5 = u0; 
        end
        
        c(it) = c_tmp;
        x_sr(it) = x_sr_tmp;
    end
    
   
    c_out = c;
    x_sr_out = x_sr;
    t_out = 0:delta:it_max*delta-delta;
end
