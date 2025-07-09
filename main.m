% while(1)

clear all
close all


global N;
global d;
global Krep;
global Katt;
global Kal;

global x_min;
global x_max;
global y_min;
global y_max;

x_min = -5;
x_max = 5;
y_min = -5;
y_max = 5;

N = 50;  % Nb d'oiseaux
d = 0.5;
Krep = 0.01; %0.01
Katt = 0.5; %0.5
Kal = 2; %2


x_init = 4*rand(N,1)-2;
y_init = 4*rand(N,1)-2;
vx_init = 8*rand(N,1)-4;
vy_init = 8*rand(N,1)-4;

% x_init = [-4;4];
% y_init = [-4.8;4.8];
% vx_init = [1;-0.5];
% vy_init = [1;0];

tstart = 0;
tfinal = 50;


X_init = zeros(4*N,1);
for i = 1:N
    X_init(((i-1)*4+1):((i-1)*4+4)) = [x_init(i);y_init(i);vx_init(i);vy_init(i)];
end

options = odeset('Events',@events,'MaxStep',1e-2);

tout = tstart;
yout = X_init.';
teout = [];
yeout = [];
ieout = [];

while tout(length(tout)) < tfinal
    tspan = [tstart tfinal];

    [ts,ys,te,ye,ie] = ode45(@dynamique2,tspan,X_init,options);

    tout = [tout; ts(2:length(ts))];
    yout = [yout; ys(2:length(ts),:)];
    teout = [teout; te];                % Events at tstart are never reported.
    yeout = [yeout; ye];
    ieout = [ieout; ie];

    for i = 1:N
        if ys(length(ts),(i-1)*4+1) <= x_min && ys(length(ts),(i-1)*4+3) < 0
            ys(length(ts),(i-1)*4+1) = x_max;
        end
        if ys(length(ts),(i-1)*4+1) >= x_max && ys(length(ts),(i-1)*4+3) > 0
            ys(length(ts),(i-1)*4+1) = x_min;
        end
        if ys(length(ts),(i-1)*4+2) <= y_min && ys(length(ts),(i-1)*4+4) < 0
            ys(length(ts),(i-1)*4+2) = y_max;
        end
        if ys(length(ts),(i-1)*4+2) >= y_max && ys(length(ts),(i-1)*4+4) > 0
            ys(length(ts),(i-1)*4+2) = y_min;
        end
    end
    
    X_init = ys(length(ts),:);

    tstart = ts(length(ts));
end


figure('WindowState','maximized')
hold on

cmap = hsv(N);
rectangle('Position',[x_min y_min x_max-x_min y_max-y_min]);

h = gobjects(N,1);
for i = 1:N
    h(i) = animatedline('Color',cmap(i,:),'LineWidth',1.5,'MaximumNumPoints',1,'Marker','.','MarkerSize',25);
end

axis([x_min x_max y_min y_max])
i = 1;
for k = 0:0.02:tout(length(tout))
    while tout(i) < k && i+1 <= length(tout)
        i = i+1;
    end
    for j = 1:N
        addpoints(h(j),yout(i,(j-1)*4+1),yout(i,(j-1)*4+2));
    end
    drawnow
    title(sprintf("t = %g",tout(i)))
    % pause(0.01)
end
for j = 1:N
    addpoints(h(j),yout(length(tout),(j-1)*4+1),yout(length(tout),(j-1)*4+2));
end
drawnow
title(sprintf("Simu finie (t = %g)",tout(length(tout))))

% figure()
% hold on
% plot(ts,ys(:,3))
% plot(ts,ys(:,7))
% plot(ts,ys(:,11))
% title('Vx')
% 
% figure()
% hold on
% plot(ts,ys(:,4))
% plot(ts,ys(:,8))
% plot(ts,ys(:,12))
% title('Vy')

% end

function [value,isterminal,direction] = events(t,y)
    global N;
    global x_min;
    global x_max;
    global y_min;
    global y_max;
    value = zeros(4*N,1);
    isterminal = zeros(4*N,1);
    direction = zeros(4*N,1);
    for i = 1:N
        value(4*(i-1)+1) = y((i-1)*4+1)-x_min;    % x
        value(4*(i-1)+2) = y((i-1)*4+1)-x_max;
        value(4*(i-1)+3) = y((i-1)*4+2)-y_min;    % y
        value(4*(i-1)+4) = y((i-1)*4+2)-y_max;
        isterminal(4*(i-1)+1:4*(i-1)+4) = [1;1;1;1];
        direction(4*(i-1)+1:4*(i-1)+4) = [-1;1;-1;1];
    end
end

