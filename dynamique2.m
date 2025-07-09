function dx = dynamique2(t,x)
    
    global N;
    global d;
    global Krep;
    global Katt;
    global Kal;

    global x_min;
    global x_max;
    global y_min;
    global y_max;

    L = x_max - x_min;
    W = y_max - y_min;
    
    x_tmp = x;

    x = zeros(N,1);
    y = zeros(N,1);
    vx = zeros(N,1);
    vy = zeros(N,1);
    
    for i = 1:N
        vy(i) = x_tmp((i-1)*4+4);
        vx(i) = x_tmp((i-1)*4+3);
        y(i) = x_tmp((i-1)*4+2);
        x(i) = x_tmp((i-1)*4+1);
    end
    
    bar_x = zeros(N,1);
    bar_x(:,:,2) = zeros(N,1);
    bar_y = zeros(N,1);
    bar_y(:,:,2) = zeros(N,1);
    bar_vx = zeros(N,1);
    bar_vx(:,:,2) = zeros(N,1);
    bar_vy = zeros(N,1);
    bar_vy(:,:,2) = zeros(N,1);


    dUrep_x = zeros(N,1);
    dUrep_y = zeros(N,1);
    dUatt_x = zeros(N,1);
    dUatt_y = zeros(N,1);
    dUal_vx = zeros(N,1);
    dUal_vy = zeros(N,1);

    
    dist = zeros(N,N);
    dist(:,:,2) = zeros(N,N);

    % Calcul de la matrice des distances
    for i = 1:N
        for j = i:N
            if i ~= j
                dist_n = norm([x(i);y(i)]-[x(j);y(j)]);
                if x(i) <= x(j)
                    dist_x = norm([x(i);y(i)]-[x(j)-L;y(j)]);
                else
                    dist_x = norm([x(i);y(i)]-[x(j)+L;y(j)]);
                end
                if y(i) <= y(j)
                    dist_y = norm([x(i);y(i)]-[x(j);y(j)-W]);
                else
                    dist_y = norm([x(i);y(i)]-[x(j);y(j)+W]);
                end
                if dist_n <= dist_x
                    if dist_n <= dist_y
                        dist(i,j,1) = dist_n;
                        dist(i,j,2) = 1;
                    else
                        dist(i,j,1) = dist_y;
                        dist(i,j,2) = 3;
                    end
                else
                    if dist_x <= dist_y
                        dist(i,j,1) = dist_x;
                        dist(i,j,2) = 2;
                    else
                        dist(i,j,1) = dist_y;
                        dist(i,j,2) = 3;
                    end
                end
                dist(j,i,1) = dist(i,j,1);
                dist(j,i,2) = dist(i,j,2);
            end
        end
    end

    for i = 1:N
        for j = i+1:N
            if dist(i,j,1) <= d
                if dist(i,j,2) == 1
                    dUrep_x(i) = dUrep_x(i) + Krep*(1/dist(j,i,1)-1/d)*((x(j)-x(i))/(dist(j,i,1)^3));
                    dUrep_y(i) = dUrep_y(i) + Krep*(1/dist(j,i,1)-1/d)*((y(j)-y(i))/(dist(j,i,1)^3));
                    dUrep_x(j) = dUrep_x(j) + Krep*(1/dist(i,j,1)-1/d)*((x(i)-x(j))/(dist(i,j,1)^3));
                    dUrep_y(j) = dUrep_y(j) + Krep*(1/dist(i,j,1)-1/d)*((y(i)-y(j))/(dist(i,j,1)^3));
                    bar_x(i,1,1) = bar_x(i,1,1) + x(j);
                    bar_x(i,1,2) = bar_x(i,1,2) + 1;
                    bar_y(i,1,1) = bar_y(i,1,1) + y(j);
                    bar_y(i,1,2) = bar_y(i,1,2) + 1;
                    bar_x(j,1,1) = bar_x(j,1,1) + x(i);
                    bar_x(j,1,2) = bar_x(j,1,2) + 1;
                    bar_y(j,1,1) = bar_y(j,1,1) + y(i);
                    bar_y(j,1,2) = bar_y(j,1,2) + 1;
                elseif dist(i,j,2) == 2
                    if x(i) <= x(j)
                        dUrep_x(i) = dUrep_x(i) + Krep*(1/dist(j,i,1)-1/d)*((x(j)-L-x(i))/(dist(j,i,1)^3));
                        dUrep_x(j) = dUrep_x(j) + Krep*(1/dist(i,j,1)-1/d)*((x(i)+L-x(j))/(dist(i,j,1)^3));
                        bar_x(i,1,1) = bar_x(i,1,1) + x(j)-L;
                        bar_x(i,1,2) = bar_x(i,1,2) + 1;
                        bar_x(j,1,1) = bar_x(j,1,1) + x(i)+L;
                        bar_x(j,1,2) = bar_x(j,1,2) + 1;
                    else
                        dUrep_x(i) = dUrep_x(i) + Krep*(1/dist(j,i,1)-1/d)*((x(j)+L-x(i))/(dist(j,i,1)^3));
                        dUrep_x(j) = dUrep_x(j) + Krep*(1/dist(i,j,1)-1/d)*((x(i)-L-x(j))/(dist(i,j,1)^3));
                        bar_x(i,1,1) = bar_x(i,1,1) + x(j)+L;
                        bar_x(i,1,2) = bar_x(i,1,2) + 1;
                        bar_x(j,1,1) = bar_x(j,1,1) + x(i)-L;
                        bar_x(j,1,2) = bar_x(j,1,2) + 1;
                    end
                    dUrep_y(i) = dUrep_y(i) + Krep*(1/dist(j,i,1)-1/d)*((y(j)-y(i))/(dist(j,i,1)^3));
                    dUrep_y(j) = dUrep_y(j) + Krep*(1/dist(i,j,1)-1/d)*((y(i)-y(j))/(dist(i,j,1)^3));
                    bar_y(i,1,1) = bar_y(i,1,1) + y(j);
                    bar_y(i,1,2) = bar_y(i,1,2) + 1;
                    bar_y(j,1,1) = bar_y(j,1,1) + y(i);
                    bar_y(j,1,2) = bar_y(j,1,2) + 1;
                else
                    if y(i) <= y(j)
                        dUrep_y(i) = dUrep_y(i) + Krep*(1/dist(j,i,1)-1/d)*((y(j)-W-y(i))/(dist(j,i,1)^3));
                        dUrep_y(j) = dUrep_y(j) + Krep*(1/dist(i,j,1)-1/d)*((y(i)+W-y(j))/(dist(i,j,1)^3));
                        bar_y(i,1,1) = bar_y(i,1,1) + y(j)-L;
                        bar_y(i,1,2) = bar_y(i,1,2) + 1;
                        bar_y(j,1,1) = bar_y(j,1,1) + y(i)+L;
                        bar_y(j,1,2) = bar_y(j,1,2) + 1;
                    else
                        dUrep_y(i) = dUrep_y(i) + Krep*(1/dist(j,i,1)-1/d)*((y(j)+W-y(i))/(dist(j,i,1)^3));
                        dUrep_y(j) = dUrep_y(j) + Krep*(1/dist(i,j,1)-1/d)*((y(i)-W-y(j))/(dist(i,j,1)^3));
                        bar_y(i,1,1) = bar_y(i,1,1) + y(j)+L;
                        bar_y(i,1,2) = bar_y(i,1,2) + 1;
                        bar_y(j,1,1) = bar_y(j,1,1) + y(i)-L;
                        bar_y(j,1,2) = bar_y(j,1,2) + 1;
                    end
                    dUrep_x(i) = dUrep_x(i) + Krep*(1/dist(j,i,1)-1/d)*((x(j)-x(i))/(dist(j,i,1)^3));
                    dUrep_x(j) = dUrep_x(j) + Krep*(1/dist(i,j,1)-1/d)*((x(i)-x(j))/(dist(i,j,1)^3));
                    bar_x(i,1,1) = bar_x(i,1,1) + x(j);
                    bar_x(i,1,2) = bar_x(i,1,2) + 1;
                    bar_x(j,1,1) = bar_x(j,1,1) + x(i);
                    bar_x(j,1,2) = bar_x(j,1,2) + 1;
                end
                bar_vx(i,1,1) = bar_vx(i,1,1) + vx(j);
                bar_vx(i,1,2) = bar_vx(i,1,2) + 1;
                bar_vy(i,1,1) = bar_vy(i,1,1) + vy(j);
                bar_vy(i,1,2) = bar_vy(i,1,2) + 1;
                bar_vx(j,1,1) = bar_vx(j,1,1) + vx(i);
                bar_vx(j,1,2) = bar_vx(j,1,2) + 1;
                bar_vy(j,1,1) = bar_vy(j,1,1) + vy(i);
                bar_vy(j,1,2) = bar_vy(j,1,2) + 1;
            end
        end
    end

    for i = 1:N
        if bar_x(i,1,2) ~= 0
            dUatt_x(i) = -Katt*(bar_x(i,1,1)/bar_x(i,1,2)-x(i));
        end
        if bar_y(i,1,2) ~= 0
            dUatt_y(i) = -Katt*(bar_y(i,1,1)/bar_y(i,1,2)-y(i));
        end
        if bar_vx(i,1,2) ~= 0
            dUal_vx(i) = -Kal*(bar_vx(i,1,1)/bar_vx(i,1,2)-vx(i));
        end
        if bar_vy(i,1,2) ~= 0
            dUal_vy(i) = -Kal*(bar_vy(i,1,1)/bar_vy(i,1,2)-vy(i));
        end
    end

    dx = zeros(N,1);
    dy = zeros(N,1);
    dvx = zeros(N,1);
    dvy = zeros(N,1);
    
    for i = 1:N
        dx(i) = vx(i);
        dy(i) = vy(i);
        dvx(i) = -dUrep_x(i)-dUatt_x(i)-dUal_vx(i); %-eps*vx(i);
        dvy(i) = -dUrep_y(i)-dUatt_y(i)-dUal_vy(i); %-eps*vy(i);
    end
    
    dx_tmp = zeros(4*N,1);
    for i = 1:N
        dx_tmp(((i-1)*4+1):((i-1)*4+4)) = [dx(i);dy(i);dvx(i);dvy(i)];
    end

    dx = dx_tmp;

end