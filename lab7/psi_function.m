function psi = psi_function(tmp_psi, nx, ny, Q_wy, Q_we, y, mu, j1, i1)
    tmp_psi(j1+1:ny+1,1) = Q_we/(2*mu)*(y(j1+1:ny+1).^3/3-y(j1+1:ny+1).^2/2*(y(j1+1)+y(ny+1))+y(j1+1:ny+1)*y(j1+1)*y(ny+1));
    tmp_psi(:,nx+1) = Q_wy(Q_we,y)/(2*mu)*(y.^3/3-y.^2/2*y(ny+1))+(Q_we*y(j1+1)^2*(-y(j1+1)+3*y(ny+1)))/(12*mu);
    tmp_psi(ny+1,2:nx) = tmp_psi(ny+1,1);
    tmp_psi(1,i1+1:nx) = tmp_psi(j1+1,1);
    tmp_psi(2:j1+1,i1+1) = tmp_psi(j1+1,1);
    tmp_psi(j1+1,2:i1+1) = tmp_psi(j1+1,1);

    psi = tmp_psi;
end