%ADA SU SAATCI - 511231179
%% Cleaning
clc
clear all
close all
%% Parameter Decleration
theta = 41;
radius = 81;

t = theta+1;              
r = radius+1; 
number_of_grids = r*t;                

r_in = 1;  
r_out = 2;                

theta_min = 0;           
theta_max = 2*pi;         

d_radi = (r_out-r_in)/(radius);   
d_theta = (2*pi)/(theta); 

constant = 15; %1.2/0.08
     
radial_grid = r_in:d_radi:r_out;            
theta_grid = theta_min:d_theta:theta_max; 

%% Domain Discritization
% Initialize grids
cartesian_radial_grid = zeros(r, t);
cartesian_theta_grid = zeros(r, t);

x = zeros(1, r);
y = zeros(1, t);

analytic = zeros(r, t);

% Populate the grids and calculate values
for j = 1:t
    for i = 1:r
        % Compute radial and angular coordinates
        rr = r_in + (i-1) * d_radi;
        tr = theta_min + (j-1) * d_theta;

        % Assign values to radial and theta grids
        cartesian_radial_grid(i, j) = rr;      
        cartesian_theta_grid(i, j) = tr;

        % Convert polar coordinates to Cartesian
        x(i) = rr * cos(tr);
        y(j) = rr * sin(tr);

        % Calculate the analytic function
        analytic(i, j) = x(i)^2 + y(j)^2;  
    end
end

cartesian_x = zeros(r,t);
cartesian_y = zeros(r,t);

for i=1:r
    for j=1:t
        cartesian_x(i,j) = cartesian_radial_grid(i,j)*cos(cartesian_theta_grid(i,j));
        cartesian_y(i,j) = cartesian_radial_grid(i,j)*sin(cartesian_theta_grid(i,j));
    end 
end

x=reshape(cartesian_x, [r*t,1]);
y=reshape(cartesian_y, [r*t,1]);


for i=1:r
    for j=1:t
        Ua(i,j)=(radial_grid(i)-1/radial_grid(i))*sin(theta_grid(j));
        Ua(1,j)=0;
        Ua(r,j)=(2-1/2)*sin(theta_grid(j));
    end
end

%% Mesh Generation
% Initialize node number matrix
node_numbers = zeros(radius * theta, 4);     

% Set initial values
node_numbers(1, 1) = 1;
index = 1;

% First column of node numbers
for i = 2:(radius * theta)
    index = index + 1;
    if mod(i - 1, radius) ~= 0
        node_numbers(i, 1) = index;                    
    else          
        node_numbers(i, 1) = index + 1;
        index = index + 1;   
    end
end

% Set initial values for second column
node_numbers(1, 2) = 2;
index = 2;

% Second column of node numbers
for i = 2:(radius * theta)
    index = index + 1;
    if mod(i - 1, radius) ~= 0
        node_numbers(i, 2) = index;                     
    else          
        node_numbers(i, 2) = index + 1;
        index = index + 1;
    end
end

% Set initial values for third column
node_numbers(1, 3) = radius + 2;
index = radius + 2;

% Third column of node numbers
for i = 2:(radius * theta)
    index = index + 1;
    if mod(i - 1, radius) ~= 0
        node_numbers(i, 3) = index;                     
    else         
        node_numbers(i, 3) = index + 1;
        index = index + 1;    
    end
end

% Set initial values for fourth column
node_numbers(1, 4) = radius + 1;
index = radius + 1;

% Fourth column of node numbers
for i = 2:(radius * theta)
    index = index + 1;
    if mod(i - 1, radius) ~= 0
        node_numbers(i, 4) = index;                    
    else        
        node_numbers(i, 4) = index + 1;
        index = index + 1;   
    end
end

% Adjust third column for the last rows
index = 1;
for i = (radius * theta - radius + 2):(radius * theta)
    index = index + 1;       
    node_numbers(i, 3) = index;
end

% Adjust fourth column for the last rows
index = 0;
for i = (radius * theta - radius + 2):(radius * theta)
    index = index + 1;  
    node_numbers(i, 4) = index;
end


%% Finding Local Stiffness Matrix and the Global Stiffness Matrix

% Initialize sparse matrix and element counter
A = sparse(number_of_grids, number_of_grids); 
current_element = 0;


while current_element < (radius * theta)
    % Node numbers and coordinates for the current element
    nodes = node_numbers(current_element + 1, :); 
    x_coords = x(nodes);
    y_coords = y(nodes);
    
    % Intermediate variables
    xx = x_coords;
    yy = y_coords;
    weights = [1, 1];
    integrals = cell(2, 2);

    % Gauss points
    for m = 1:2
        for n = 1:2
            % Initialize stiffness matrix
            stiffness_matrix = zeros(4);
            ksi = (-1)^(m) / sqrt(3);
            eta = (-1)^(n) / sqrt(3); 

            % Shape function derivatives 
            dN_dksi = [-1-eta, -1+eta, 1-eta, 1+eta] * 0.25;
            dN_deta = [1-ksi, -1+ksi, -1-ksi, 1+ksi] * 0.25;
            
            % Jacobian components
            dx_dksi = 0; 
            dx_deta = 0; 
            dy_dksi = 0; 
            dy_deta = 0; 

            % Calculate Jacobian components
            for node = 1:4
                dx_dksi = dx_dksi + xx(node) * dN_dksi(node);
                dx_deta = dx_deta + xx(node) * dN_deta(node);
                dy_dksi = dy_dksi + yy(node) * dN_dksi(node);
                dy_deta = dy_deta + yy(node) * dN_deta(node);
            end
            
            % Compute Jacobian determinant
            Jacobian = abs(dx_dksi * dy_deta - dx_deta * dy_dksi);

            % Compute inverses of Jacobian matrix components
            dksi_dx = dy_deta / Jacobian;
            dksi_dy = -dx_deta / Jacobian;
            deta_dx = -dy_dksi / Jacobian;
            deta_dy = dx_dksi / Jacobian;

            % Calculate local stiffness matrix
            for i = 1:4
                for j = 1:4
                    stiffness_matrix(i, j) = ((dN_dksi(i) * dksi_dx + dN_deta(i) * deta_dx) * (dN_dksi(j) * dksi_dx + dN_deta(j) * deta_dx) + ...
                                              (dN_dksi(i) * dksi_dy + dN_deta(i) * deta_dy) * (dN_dksi(j) * dksi_dy + dN_deta(j) * deta_dy)) * Jacobian;
                end
            end

            % Store weighted local stiffness matrix
            integrals{m, n} = stiffness_matrix * weights(m) * weights(n);
        end
    end

    % Combine contributions 
    local_stiffness_matrix = integrals{1, 1} + integrals{1, 2} + integrals{2, 1} + integrals{2, 2};  
    
    % Assemble global stiffness matrix
    for i = 1:4
        for j = 1:4
            A(nodes(i), nodes(j)) = A(nodes(i), nodes(j)) + local_stiffness_matrix(i, j);
        end
    end

    % Increment element counter
    current_element = current_element + 1;
end


%%  Boundary Conditions
% Calculate the total number of nodes
num_nodes = r * t; 

% Initialize the RHS vector
rhs_vector = zeros(num_nodes, 1); 

% Populate matrix A with boundary conditions
for row = 1:r
    for col = 1:t
        index = row + (col - 1) * r;   
        if row == 1      
            A(index, index) = 1;  % Boundary at the inner radius
        elseif row == r    
            A(index, index) = (1.5) * sin(theta_grid(col));  % Boundary at the outer radius  
        elseif col == 1     
            A(index, index) = 1;  % Boundary at the start of theta 
        elseif col == t   
            A(index, index) = 1;  % Boundary at the end of theta   
        end
    end
end

% Apply boundary conditions to the last row of the grid
for col = 1:t
    index = r + (col - 1) * r;   
    A(index, index) = 1;
    rhs_vector(index, 1) = 1.5 * sin(theta_grid(col));
end

%% LU Decomposition

% Perform LU decomposition with partial pivoting
[L, U, P] = lu(A);

% Solve the system of equations
intermediate_solution = L \ (P * rhs_vector);
psi = U \ intermediate_solution;

% Initialize the solution matrix
solution_grid = zeros(r, t);

% The solution grid
for row = 1:r
    for col = 1:t
        index = row + (col - 1) * r;
        solution_grid(row, col) = psi(index);
    end
end


% USER MUST CHANGE BELOW ACCORDING TO THE MESH SIZE!!
rearranged_solution = zeros(82, 41);

for row = 1:82
    for col = 1:21
        rearranged_solution(row, col) = solution_grid(row, 43 - col);
        rearranged_solution(row, 43 - col) = solution_grid(row, col);
    end
end
 
rearranged_solution = constant * rearranged_solution;

%% Plots
radial_grid = linspace(1,2,r); 
z2 = linspace(0,360,t); 

figure
polarPcolor(radial_grid,z2,rearranged_solution); 
title('LU Decomposition')

figure
polarPcolor(radial_grid,z2,Ua); 
title('Analytical Solution')

%% Error
% Initialize error matrix
error_matrix = zeros(r, t);

% Compute error matrix
for i = 1:r
    for j = 1:t
        error_matrix(i, j) = norm(Ua(i, j) - rearranged_solution(i, j), 2) / sqrt(r * t);
    end
end

% Calculate the overall error norm
overall_error_norm = norm(error_matrix, 2) * 100;

error_diff_1 = (error_matrix(3, 3) - error_matrix(2, 2)) / d_radi;

% Initialize arrays for error and radial differences
errors = zeros(2, 1);
radial_diffs = zeros(2, 1);


errors(1) = error_matrix(19, 19);
errors(2) = error_matrix(18, 18);

radial_diffs(1) = cartesian_radial_grid(29, 29);
radial_diffs(2) = cartesian_radial_grid(28, 28);

slope = log(errors(2) - errors(1)) / log(radial_diffs(2) - radial_diffs(1));
fprintf("The slope is %e \n", slope);

figure
loglog(radial_diffs,errors); 
title('Error Graph')
xlabel('r')
ylabel('Error')