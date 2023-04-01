function drawSpeedField(x,y,vx,vy)
    figure;
    s = pcolor(x, y, vx);
    set(s, 'EdgeColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('x');
    ylabel('y');
    title("V_x");

    figure;
    s = pcolor(x, y, vy);
    set(s, 'EdgeColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('x');
    ylabel('y');
    title("V_y");
end