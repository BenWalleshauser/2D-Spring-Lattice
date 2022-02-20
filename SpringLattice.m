%B.W.
%2/17/22

%Simulating 2-D spring lattice w/. forward Euler
%%

%Time
t_total = 100;
h = 0.001;
time_step = 0:h:t_total;

%Spring constant
k = 10;

%Spacing between nodes
spacing = 10;

%Mass of particles
m = 0.001;

%Number of lattice points
n = 10;
m = 10;

%Displacements in the x-direction of the atoms (adding two in each
%direction as a solid boundary)
X = zeros(n+2,m+2,length(t_total));

%Displacement in the y-direction of the atoms
Y = zeros(n+2,m+2,length(t_total));

%Perturb a node
X(n/2,m/2) = 5;
Y(n/2,m/2) = 4;

%% Forward Euler

%New notation helps uncouple 2nd order ODES 
Z = zeros(n+2,m+2,length(t_total));
W = zeros(n+2,m+2,length(t_total));

for t = 1:length(time_step)
    for i = 2:1:n-1
        for j = 2:1:m-1
            %-------X Component---------
            del_top = (X(i-1,j,t) - X(i,j,t))*(1-spacing/sqrt(spacing^2 + (X(i-1,j,t) - X(i,j,t))^2));
            del_bottom = (X(i+1,j,t) - X(i,j,t))*(1-spacing/sqrt(spacing^2 + (X(i+1,j,t) - X(i,j,t))^2));
            Z(i,j,t+1) = Z(i,j,t) + h*k/m*(X(i,j+1,t) - 2*X(i,j,t) + X(i,j-1,t) + del_top + del_bottom);
            X(i,j,t+1) = X(i,j,t) + h*Z(i,j,t);
            
            %------Y Component-----------
            del_right = (Y(i,j+1,t) - Y(i,j,t))*(1-spacing/sqrt(spacing^2 + (Y(i+1,j,t) - Y(i,j,t))^2));
            del_left = (Y(i,j-1,t) - Y(i,j,t))*(1-spacing/sqrt(spacing^2 + (Y(i-1,j,t) - Y(i,j,t))^2));
            W(i,j,t+1) = W(i,j,t) + h*k/m*(Y(i+1,j,t) - 2*Y(i,j,t) + Y(i-1,j,t));
            Y(i,j,t+1) = Y(i,j,t) + h*W(i,j,t);
        end
    end
end

%% Plotting

x_vec = zeros((n-2)*(m-2),length(time_step));
y_vec = zeros((n-2)*(m-2),length(time_step));
it = 1;

for t = 1:length(time_step)  
    for i = 2:1:n-1
        for j = 2:1:m-1
            x_vec(it,t) = X(i,j,t)+j*spacing;
            it = it+1;
        end
    end
    it = 1;
end

for t = 1:length(time_step)  
    for i = 2:1:n-1
        for j = 2:1:m-1
            y_vec(it,t) = Y(i,j,t) + m*spacing - i*spacing;
            it = it+1;
        end
    end
    it = 1;
end

close all
for i = 1:1/(5*h):length(time_step)
   figure(1)
   plot(x_vec(:,i),y_vec(:,i),'o') 
   xlim([spacing n*spacing])
   ylim([0 m*spacing-spacing])
   pause(0.01); 
end
























