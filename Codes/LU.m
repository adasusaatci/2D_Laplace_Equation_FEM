%% Mehmet Emin Güney-110170031
% LU Factorization
clc
clear all
close all
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

% LU Factorization
[L,U,P] = lu(A);
B = U \ ( L \ (P*RHS) );

for i=1:r
    for j=1:theta
        k=i+(j-1)*r;
        Uij(i,j)=B(k);
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
title('LU')
figure
polarPcolor(R,Az,Ua); % col
title('Analytical Solution')


% error1=abs(Uij(20,20)-Ua(20,20))/sqrt(20*20);
% error2=abs(Uij(40,40)-Ua(40,40))/sqrt(40*40);
% error3=abs(Uij(81,41)-Ua(81,41))/sqrt(81*41);
% error=[error1 error2 error3];
% figure
% loglog(error)