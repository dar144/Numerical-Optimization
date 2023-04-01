function [vx_out,vy_out] = speedField(psi, nx, ny, delta, i1, i2, j1)
    vx = zeros(ny+1,nx+1);
    vy = zeros(ny+1,nx+1);

    vx(2:ny,2:nx) = (psi(3:ny+1,2:nx)-psi(1:ny-1,2:nx))/(2*delta);
    vy(2:ny,2:nx) = -(psi(2:ny,3:nx+1)-psi(2:ny,1:nx-1))/(2*delta);
    
    vx(1:j1+1,i1+1:i2+1) = 0;
    vy(1:j1+1,i1+1:i2+1) = 0;

    vx(1,2:nx) = 0;
    vy(ny+1,2:nx) = 0;

%     vx(1:ny,1) = vx(1:ny,2);
%     vy(1:ny+1,nx+1) = vy(1:ny+1,nx);

    vx(1:ny,1) = vx(1:ny,2);
    vx(1:ny+1,nx+1) = vx(1:ny+1,nx);

    vy(1:ny,1) = vy(1:ny,2);
    vy(1:ny+1,nx+1) = vy(1:ny+1,nx);
    
    vx_out = vx;
    vy_out = vy;
end