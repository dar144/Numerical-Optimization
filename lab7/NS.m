function [psi_out,zeta_out, u_out, v_out] = NS(nx, ny, Q_wy, Q_we, psi, y, zeta, j1, i1, IT_MAX, delta, ro, mu)
    tmp_psi = psi_function(psi, nx, ny, Q_wy, Q_we, y, mu, j1, i1);
    tmp_zeta = zeta;
    u = zeros(ny+1,nx+1);
    v = zeros(ny+1,nx+1);
     
    for it = 1:IT_MAX
        it
        if it < 2000
            omega = 0;
        else
            omega = 1;
        end
        
        for i = 2:nx
            for j = 2:ny
                if is_edge(i,j,nx,ny,i1,j1) == false
                    tmp_psi(j,i) = (tmp_psi(j,i+1)+tmp_psi(j,i-1)+tmp_psi(j+1,i)+tmp_psi(j-1,i)-delta^2*tmp_zeta(j,i))/4;
                    tmp_zeta(j,i) = (tmp_zeta(j,i+1)+tmp_zeta(j,i-1)+tmp_zeta(j+1,i)+tmp_zeta(j-1,i))/4 - omega*ro/(16*mu)*...
                    ((tmp_psi(j+1,i)-tmp_psi(j-1,i))*(tmp_zeta(j,i+1)-tmp_zeta(j,i-1))-(tmp_psi(j,i+1)-tmp_psi(j,i-1))*(tmp_zeta(j+1,i)-tmp_zeta(j-1,i)));
                end	

                u(j,i) = (tmp_psi(j+1,i) - tmp_psi(j-1,i))/(2*delta);
                v(j,i) = -(tmp_psi(j,i+1) - tmp_psi(j,i-1))/(2*delta);     
            end
        end
        tmp_zeta = zeta_function(tmp_zeta, tmp_psi, nx, ny, Q_wy, Q_we, y, mu, j1, i1, delta);       
    end
    
    psi_out = tmp_psi;
    zeta_out = tmp_zeta;
    v_out = v;    
    u_out = u;
end