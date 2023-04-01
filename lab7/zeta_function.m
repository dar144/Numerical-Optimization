function zeta = zeta_function(tmp_zeta, psi, nx, ny, Q_wy, Q_we, y, mu, j1, i1, delta)
    c = 2/delta^2;
    tmp_zeta(j1+1:ny+1,1) = Q_we/(2*mu)*(2.*y(j1+1:ny+1)-y(j1+1)-y(ny+1));  
    tmp_zeta(:,nx+1) = Q_wy(Q_we,y)/(2*mu)*(2.*y-y(ny+1));
    tmp_zeta(ny+1,2:nx) = c*(psi(ny,2:nx)-psi(ny+1,2:nx));
    tmp_zeta(1,i1+2:nx) = c*(psi(2,i1+2:nx)-psi(1,i1+2:nx));
    tmp_zeta(2:j1,i1+1) = c*(psi(2:j1,i1+2)-psi(2:j1,i1+1));
    tmp_zeta(j1+1,2:i1+1) = c*(psi(j1+2,2:i1+1)-psi(j1+1,2:i1+1));
    tmp_zeta(j1+1,i1+1) = (tmp_zeta(j1+1,i1)+tmp_zeta(j1,i1+1))/2;
    
    zeta = tmp_zeta;
end