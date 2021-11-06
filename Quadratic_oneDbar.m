clc
clear all
disp('RN')
disp('***\n FEM ANALYSIS OF 1-D Cantilever BAR WITH UDL and 3 nodes per element ==>>>>>>>');
fprintf('\n \n');

L = input('Enter the length of the element(in mm) :');
e = input('\n Enter the number of elements: ');           
Ar = input('\n Enter the area of the element(in mm^2) :');
E = input('\n Enter the value of E (in N/mm^2):');        % Modulus of Elasticity         

% nodes = degree of freedom = size of global stiffness matrix
nodes = (e*2)+1 ;                         % number of nodes
le = L/e;                            % length of each element   


C_M = zeros(e,3);                   % Matrix Connectivity
for i=1:e
    for j=2
    C_M(i,1)= i + (i-1);
    C_M(i,j)= C_M(i,1) + 2 ;
    C_M(i,j+1) = C_M(i,j) - 1;
    end
end

% Generating force matrix
q = input("Enter the value of force in N/mm: ");
Fg = zeros(3,1);
%Global elemental load matrix 1 3 2
for s=1:3
   Fg(s,1) = Fg(s,1) + le*q*(1/6);
end
Fg(3,1) = Fg(3,1)*4;
disp('Global Elemental Force Matrix with index 1 3 2 is:')
disp(Fg);                            % elemental load matrix

%Global Elimental Stiffness Matrix with Index as 1 3 2
c = E*Ar/(3*le);
Kg = c*[7 1 -8;1 7 -8;-8 -8 16];
disp('Global elemental stiffness Matrix with index 1 3 2');
disp(Kg);

%Global Stiffness Matrix with Index as 1 2 3
Ke = zeros(3,3);
for j = 1:3
    for k=1:3
    Ke(j,k) = Ke(j,k) + Kg(C_M(1,j),C_M(1,k));
    end
end
disp('Elemental Stiffness Matrix Ke with index 1 2 3:');
disp(Ke);

%Global to Local Force Matrix with Index 1 2 3
Fe = zeros(3,1);
for r=1:3
    Fe(r,1) = Fe(r,1) + Fg(C_M(1,r));
end
disp('Global Elemental Force Matrix Fe with index 1 2 3');
disp(Fe)

%Assembling Stiffness Matrix
K = zeros(nodes,nodes);
for s=1:e
    value = [C_M(s):C_M(s,2)];
    K(value,value) = K(value,value) + Ke;
end
disp('Assembeled Stiffness Matrix:');
disp(K);

%Assembling Force Matrx Force Matrix
F = zeros(nodes,1);
for t=1:e
    value = [C_M(t):C_M(t,2)];
    F(value,1) = F(value,1) + Fe;
end
disp('Assembled Force Matrix:');
disp(F);


fix_n = [1];                 %declaration of fix node 1
Nodes = 1:nodes; 

% Matrix of free nodes
free_n = setxor(Nodes,fix_n);
K1 = K(free_n,free_n);
Fp1 = F(free_n,1);
disp('Assembled Stiffness Matrix at Free nodes');
disp(K1);
disp('Assembled Load Matrix at free nodes');
disp(Fp1);

 
%Displacement Matrix
U = zeros(nodes,1);             %Displacement defination
Up = K1\Fp1;
disp('Displacement at each node:')
U(free_n) = Up;
disp(U);
   

%Stress Matrix
B = zeros(1,3);                         %defining B matrix
S = zeros(nodes,1);                     % Stress Matrix defination
for m=1:e 
    for n=1:3
    xi = n-2;                           %variation of xi with respect to local nodes -1<xi<1
     for p=1
        B(p)=B(1,p)+(-(1-2*xi)/2);
        B(p+1)=B(p+1)+(-2*xi);
        B(p+2)=B(p+2)+((1+2*xi)/2);
     end
    S(n+(m-1)*2,1)=S(n+(m-1)*2,1) + E*B*U(C_M(m):C_M(m,2),1);       %Filling in stress values in S matrix
    B = zeros(1,3);
    end
end
for h=2:nodes-1                                                     %Dividing stress values by 2 at similar nodes as they were added twice
    if rem(h,2) ~= 0
        S(h)=S(h)/2;
    end
end
disp('Nodal Stress matrix is as follows:');
disp(S);

disp('Nodal Strain Matrix is as follows:');
Str = S/E;
disp(Str);

% Plotting displacement at each node generated using FEM
xx = 0:le/2:L;      %dividing L in nodes
uu = U;
figure;
hold on;
xlabel('Nodes')
ylabel('Displacement')
title('Displacement at each node')
plot(xx,uu,'b','LineWidth', 2);
plot( xx,uu ,'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('Displacement','Nodes','Location','East');

%analytical method
x = 0:0.1:L;
ux = (q/(2*E*Ar))*(2*L.*x - x.^2);                     %exact solution
figure;
subplot(1,2,1) ; plot(x,ux,'k','LineWidth', 1);        %analytical plot
title('Analytical')
subplot(1,2,2) ; plot(xx,uu,'b','LineWidth', 1);       %FEM plot
axis([min(xx) max(xx) min(uu) max(uu)]);
title('FEM')


% FEM vs Analytical
figure;
hold on;
title("FEM vs Analytical")
plot(x,ux,'k','LineWidth', 1); 
plot(xx,uu,'b');
plot( xx,uu ,'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('Analytical','FEM','Nodes','Location','East');

% Plotting Stress wrt x
figure;
hold on;
title("Stress Plot");
ss = S;
plot(xx,ss,'g','LineWidth',1);
plot( xx,ss ,'ro', 'MarkerSize', 10, 'LineWidth', 2);
title("Stress vs x");
xlabel('X');
ylabel('Stress');

% Plotting Strain
figure;
hold on;
sr = Str;
plot(xx,sr,'m','LineWidth',1);
plot( xx,sr ,'ro', 'MarkerSize', 10, 'LineWidth', 2);
title("Strain vs x");
xlabel('X');
ylabel('Strain');

%Strss vs Strain
%plot(ss,sr);


