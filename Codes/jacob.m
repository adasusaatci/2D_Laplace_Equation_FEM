%% Mehmet Emin Güney-110170031
% LU Factorization
clc
clear all
% Inputs
inner_radius=1;     %inner radius
outer_radius=2;    %outer radius
radius_diff=outer_radius-inner_radius;
theta_min=0;
theta_max=2*pi;

% Mesh Size
r=81;  
theta=41;  
delta_r=(radius_diff)/(r-1); %Grid
delta_theta=(theta_max-theta_min)/(theta-1); % Grid
R=inner_radius:delta_r:outer_radius;
T=theta_min:delta_theta:theta_max;
[r,theta]=meshgrid(r,theta);
% A matrix
Node=r*theta;          
A=sparse(Node);   
RHS=sparse(Node,1); % Boundary Conditions
for i=2:r-1    
    for j=2:theta-1     
        k=i+(j-1)*r;    
        A(k,k)=(-2/delta_r^2) - (2/((R(i)^2)*delta_theta^2));       %  Uij
        A(k,k-1)=(1/delta_r^2) - (1/(R(i)*2*delta_r));          %  U_(i-1,j)
        A(k,k+1)=(1/delta_r^2) + (1/(R(i)*2*delta_r));          %  U_(i+1,j)           
        A(k,k-r)=1/((R(i)^2)*delta_theta^2);             % U_(i,j-1)
        A(k,k+r)=1/((R(i)^2)*delta_theta^2);             % U_(i,j+1)     
        RHS(k,1)=1;                  
    end
end
for i=1:r
    for j=1:theta
           k=i+(j-1)*r;   
           if i==1      
               A(k,k)=1;  
           elseif i==r    
               A(k,k)=1;    
           elseif j==1     
               A(k,k)=1;   
           elseif j==theta    
               A(k,k)=1;    
           end
    end
end

% Drihlet Boundary Conditions
for j=1:theta
        k=r+(j-1)*r;   
        A(k,k)=1;
        RHS(k,1)=(2-1/2)*sin(T(j));
end

%  Jacobi

tolerance = 10^-5;

D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);


TJ = D \ (L+U);        
CJ = D \ RHS;	

iteration = 1;
B = sparse(Node,1);

while iteration <= 50000
	
	Psi2 = B;
	B = TJ * Psi2 + CJ;
    told = norm( B - Psi2 );
	if told < tolerance
		fprintf('Jacobbi converged at %d iterations', iteration)
		break
	end
	iteration = iteration + 1;
    convergeplot(1,iteration)=told;
end

for i=1:r
    for j=1:theta
        m=i+(j-1)*r;
        Uij(i,j)=B(m);
    end
end


%Analytical Solution

for i=1:r
    for j=1:theta
        Ua(i,j)=(R(i)-1/R(i))*sin(T(j));
    end
end

aaa=r;
bbb=theta;
R = linspace(1,2,aaa); % (distance in km)
Az = linspace(0,360,bbb); % in degrees
figure
polarPcolor(R,Az,Uij); % col
title('Jacobbi')
figure
polarPcolor(R,Az,Ua); % col
title('Analytical Solution')





