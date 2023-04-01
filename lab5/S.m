function sum = S(V, k, delta, nx, ny)
    tmpSum = 0;
    divider = 2*k*delta;
    for i = 1:k:nx+1-k
        for j = 1:k:ny+1-k
            tmpSum = tmpSum + (k*delta)^2/2*(((V(j,i+k)-V(j,i)+V(j+k,i+k)-V(j+k,i))/divider)^2+((V(j+k,i)-V(j,i)+V(j+k,i+k)-V(j,i+k))/divider)^2);
        end
    end
    sum = tmpSum;
end