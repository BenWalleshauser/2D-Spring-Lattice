%B.W.
%2/17/22

%Simulating 2-D spring lattice w/. forward euler
%%
clear 
clc
close all
%Time
t_total = 1200;
h = 0.01;
time_step = 0:h:t_total;

%Spring constant
k = 0.1;

%Ln between nodes
Ln = 10;

%Mass of particles
mass = 0.1;

%Number of lattice points
n = 4;
m = 4;

%Displacements in the x-direction of the atoms (adding two in each
%direction as a solid boundary)
X = zeros(n+2,m+2,length(t_total));

%Displacement in the y-direction of the atoms
Y = zeros(n+2,m+2,length(t_total));

%Perturb a node
X(n/2,m/2) = 15;
Y(n/2,m/2) = 15;

%% Forward Euler

%New notation helps uncouple 2nd order ODES 
Z = zeros(n+2,m+2,length(t_total));
W = zeros(n+2,m+2,length(t_total));

for t = 1:length(time_step)
    for i = 2:1:n+1
        for j = 2:1:m+1
            
            %right
            dx = X(i,j+1,t)-X(i,j,t)+Ln;
            dy = Y(i,j+1,t)-Y(i,j,t);
            Ls = sqrt(dx^2+dy^2);
            ds = Ls-Ln;
            F = k*ds;
            Fx_right = F*dx/Ls;
            Fy_right = F*dy/Ls;
            
            %left
            dx = X(i,j-1,t)-X(i,j,t)-Ln;
            dy = Y(i,j-1,t)-Y(i,j,t);
            Ls = sqrt(dx^2+dy^2);
            ds = Ls-Ln;
            F = k*ds;
            Fx_left = F*dx/Ls;
            Fy_left = F*dy/Ls;
            
            %top
            dx = X(i-1,j,t)-X(i,j,t);
            dy = Y(i-1,j,t)-Y(i,j,t)+Ln;
            Ls = sqrt(dx^2+dy^2);
            ds = Ls-Ln;
            F = k*ds;
            Fx_top = F*dx/Ls;
            Fy_top = F*dy/Ls;
            
            %bottom
            dx = X(i+1,j,t)-X(i,j,t);
            dy = Y(i+1,j,t)-Y(i,j,t)-Ln;
            Ls = sqrt(dx^2+dy^2);
            ds = Ls-Ln;
            F = k*ds;
            Fx_bottom = F*dx/Ls;
            Fy_bottom = F*dy/Ls;
            
            Fx = Fx_right + Fx_left + Fx_top + Fx_bottom;
            Fy = Fy_right + Fy_left + Fy_top + Fy_bottom;
            
            Z(i,j,t+1) = Z(i,j,t) + h*Fx/mass;
            X(i,j,t+1) = X(i,j,t) + h*Z(i,j,t+1);
            
            W(i,j,t+1) = W(i,j,t) + h*Fy/mass;
            Y(i,j,t+1) = Y(i,j,t) + h*W(i,j,t+1);
            
        end
    end
end

%% Plotting

x_vec = zeros((n-2)*(m-2),length(time_step));
y_vec = zeros((n-2)*(m-2),length(time_step));
it = 1;

for t = 1:length(time_step)  
    for i = 2:1:n+1
        for j = 2:1:m+1
            x_vec(it,t) = X(i,j,t)+j*Ln;
            y_vec(it,t) = Y(i,j,t) + m*Ln - i*Ln;
            it = it+1;
        end
    end
    it = 1;
end

close all
for i = 1:1/(5*h):length(time_step)
   figure(1)
   grid on
   plot(x_vec(:,i),y_vec(:,i),'o') 
   xlim([min(min(x_vec)) max(max(x_vec))])
   ylim([min(min(y_vec)) max(max(y_vec))])
   pause(0.01); 
end
