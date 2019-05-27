clear all;
clc;
% close all;
N=100;

[X,Y,Z] = sphere(N);
U=zeros((N+1)^2,3);
fval=zeros((N+1)^2,1);

for i=1:(N+1)^2
e = [X(i);Y(i);Z(i)];
U(i,:)=crossover_point([0;0;0],e);
end
U1=zeros(N+1,N+1);
U2=zeros(N+1,N+1);
U3=zeros(N+1,N+1);
for i=1:(N+1)^2
    U1(i)=U(i,1);
    U2(i)=U(i,2);
    U3(i)=U(i,3);
end
% figure,
% plot3(U(:,1),U(:,2),U(:,3),'r.');hold on;

figure,
C=sqrt(U1.*U1+U2.*U2+U3.*U3);
surf(U1,U2,U3,C);% 'FaceColor','b','FaceAlpha',0.5,,'EdgeColor','none'
% hold on;
% plot3(0,0,0.3043,'r.');hold on;
% plot3(0.3088,0,0.1521,'r.');hold on;
% plot3(0.3088,0.3088,0,'r.');hold on;
% plot3(0,0.3088,0.1521,'r.');hold on;
% plot3(0.3088,0,0,'r.');hold on;
% plot3(0,0.3088,0,'r.');hold on;
% plot3(0,0.3088,0.1521,'r.');
% x1=linspace(0,1,10);
% y1=linspace(0,1,10);
% [X1,Y1]=meshgrid(x1,y1);
% Z1=0.3043 - X1.*0.3893 - Y1.*0.3893;
% mesh(X1,Y1,Z1)
% y1=linspace(0,1,10);
% z1=linspace(0,1,10);
% [Y1,Z1]=meshgrid(y1,z1);
% X1=ones(10).*0.39088;
% mesh(X1,Y1,Z1)
% x1=linspace(0,1,10);
% z1=linspace(0,1,10);
% [X1,Z1]=meshgrid(x1,z1);
% Y1=ones(10).*0.39088;
% mesh(X1,Y1,Z1)

colorbar
% C1=sqrt(X.*X+Y.*Y+Z.*Z);
% figure,
% surf(X,Y,Z,C1);
% colorbar
