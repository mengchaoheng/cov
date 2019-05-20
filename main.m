clear all;
clc;
% close all;
N=20;
K=1*[-0.56   0       0.56    0;
          0  -0.56    0       0.56;
      0.218   0.218   0.218   0.218];
a=0.3490000;
  
[X,Y,Z] = sphere(N);
% X=X';
% Y=Y';
% Z=Z';
x=zeros(4,(N+1)^2);
fval=zeros((N+1)^2,1);
 
for i=1:(N+1)^2
e = [X(i);Y(i);Z(i)];
% e = [xyz(i,1);xyz(i,2);xyz(i,3)];
e=e/norm(e);
cf = @(x)confun(x,e,K);
of = @(x)objfun(x,K);
% A=[];b=[];Aeq=[];beq=[];
A=-e'*K;b=0;Aeq=[0 -e(3) e(2);e(3) 0 -e(1);-e(2) e(1) 0]*K;beq=[0;0;0];
lb=[-a,-a,-a,-a];ub=[a,a,a,a];
options = optimoptions(@ fmincon,'Algorithm','interior-point','Display','iter');
% x0=[0;0;0;0]; %对解决方案进行初步猜测
c1=-e(1)+0.5*e(3);c2=-e(2)+0.5*e(3);c3= e(1)+0.5*e(3);c4= e(2)+0.5*e(3);
x0=0.05*[c1;c2;c3;c4];
[x(:,i),fval(i),exitflag,output,lambda,grad,hessian] = fmincon(of,x0,A,b,Aeq,beq,lb,ub,[],options);
end
U=K*x;
U1=zeros(N+1,N+1);
U2=zeros(N+1,N+1);
U3=zeros(N+1,N+1);
for i=1:(N+1)^2
    U1(i)=U(1,i);
    U2(i)=U(2,i);
    U3(i)=U(3,i);
end
save data.txt U -ascii;
xlswrite('data1.xls',U1);
xlswrite('data2.xls',U2);
xlswrite('data3.xls',U3);

figure,
plot3(U(1,:),U(2,:),U(3,:),'r.');hold on;
% plot3(X,Y,Z,'b*');
C=sqrt(U1.*U1+U2.*U2+U3.*U3);
figure,
surf(U1,U2,U3,C);% 'FaceColor','b','FaceAlpha',0.5,,'EdgeColor','none'
hold on;
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

% shading interp;
colorbar
C1=sqrt(X.*X+Y.*Y+Z.*Z);
figure,
surf(X,Y,Z,C1);
% shading interp;
colorbar
% figure,
% CC = gradient(U3);
% meshz(U1,U2,U3,CC)


[row,col,v] = find(U(3,:)>=-0.00001);
u=zeros(3,length(col));
for i=1:length(col)
        u(1,i)=U(1,col(i));
        u(2,i)=U(2,col(i));
        u(3,i)=U(3,col(i));
end
% u=u';
% save data1.txt u -ascii;
[row1,col1,v1] = find(u(2,:)>=-0.00001);%U(3,:)u(3,:)<0.3043
u1=zeros(3,length(col1));
for i=1:length(col1)
        u1(1,i)=u(1,col1(i));
        u1(2,i)=u(2,col1(i));
        u1(3,i)=u(3,col1(i));
end

[row2,col2,v2] = find(u1(1,:)>=-0.00001);
u2=zeros(3,length(col2));
for i=1:length(col2)
        u2(1,i)=u1(1,col2(i));
        u2(2,i)=u1(2,col2(i));
        u2(3,i)=u1(3,col2(i));
end
[row3,col3,v3] = find(u2(3,:)<0.3043);
u3=zeros(3,length(col3));
for i=1:length(col3)
        u3(1,i)=u2(1,col3(i));
        u3(2,i)=u2(2,col3(i));
        u3(3,i)=u2(3,col3(i));
end

u3=u3';
save datau3.txt u3 -ascii;
% xx=u3(:,1);
% yy=u3(:,2);
% zz=u3(:,3);

% u1=u1';
[azimuth,elevation,r] = cart2sph(u3(:,1),u3(:,2),u3(:,3));%1379  50 101
xx=cos(azimuth);
yy=cos(elevation);
zz=r;
[X2,Y2]=meshgrid( linspace(min(azimuth),max(azimuth))',linspace(min(elevation),max(elevation)) );
Z2=griddata(azimuth,elevation,r,X2,Y2);%插值
figure
surf(X2,Y2,Z2);grid on;hold on;
plot3(azimuth,elevation,r,'r.');
% stem3(azimuth(1395:1401),elevation(1395:1401),r(1395:1401));
%  stem3(azimuth,elevation,r);
 
 
p_start=zeros(length(u3(:,1)),3);%向量起点
x=p_start(:,1);
y=p_start(:,2);
z=p_start(:,3);
u=u3(:,1);
v=u3(:,2);
w=u3(:,3);
 figure
quiver3(x,y,z,u,v,w,1);

% [C,ia,ic] = unique(U(1:2,:)','rows');
% C=C';
% [X2,Y2]=meshgrid( linspace(min(C(1,:)),max(C(1,:)))',linspace(min(C(2,:)),max(C(2,:))) );
% Z2=griddata(u(1,:),u(2,:),u(3,:),X2,Y2);%插值
% figure
% % surf(X2,Y2,Z2);
% C1 = gradient(Z2);

% figure
% meshz(X2,Y2,Z2,C1)
% fill3(X,Y,Z)
% waterfall(X2,Y2,Z2);
% hold on;
% scatter3(u2(1,:),u2(2,:),u2(3,:),4,'r')
% figure
% patch(surf2patch(X2,Y2,Z2,Z2)); 
% shading faceted; 
% view(3)
% figure
% [X1,Y1]=meshgrid( linspace(min(u2(1,:)),max(u2(1,:)))',linspace(min(u2(2,:)),max(u2(2,:))) );
% Z1=griddata(u2(1,:),u2(2,:),u2(3,:),X1,Y1,'v4');
% surf(X1,Y1,Z1);
% hold on;
% scatter3(u2(1,:),u2(2,:),u2(3,:),4,'r');
% figure,
% tri = delaunay(u(1,:),u(2,:));
% trisurf(tri,u(1,:),u(2,:),u(3,:));
% shading interp

% figure
% k = 5;
% n = 2^k-1;
% theta = pi*(-n:2:n)/n;
% phi = (pi/2)*(-n:2:n)'/n;
% X3 = cos(phi)*cos(theta);
% Y3 = cos(phi)*sin(theta);
% Z3 = sin(phi)*ones(size(theta));
% colormap([0 0 0;1 1 1])
% C3 = hadamard(2^k); 
% surf(X3,Y3,Z3,C3)
% axis square


% X4=U(1,1:1:length(U(1,:)));
% Y4=U(2,1:1:length(U(1,:)));
% Z4=U(3,1:1:length(U(1,:)))'*ones(size(U(3,1:1:length(U(1,:)))));
% figure
% colormap([0 0 0;1 1 1])
% C5 = hadamard(2^(N/20)); 
% surf(X4,Y4,Z4)
% axis square

% n=1000;
% r=1;
% Phi=rand(N,1);
% Theta=rand(N,1);
% X=r.*sin(Theta).*cos(Phi);
% Y=r.*sin(Theta).*sin(Phi);
% Z=r.*cos(Theta); 




% xyz = zeros(N, 3);
% for i = 1:N
% theta = 2*pi*rand;
% beta = 2*pi*rand;
% xyz(i,:) = [cos(beta)*[cos(theta) sin(theta)] sin(beta)];
% end
% figure,
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.');

% plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.');
%  plot3(U(1,:),U(2,:),U(3,:),'.','markersize',12)
% n=1000;
% r=1;
% Phi=rand(n,1);
% Theta=rand(n,1);
% X=r.*sin(Theta).*cos(Phi);
% Y=r.*sin(Theta).*sin(Phi);
% Z=r.*cos(Theta);  

% 
% U11=zeros(length(u3(1,:))/2,length(u3(1,:))/2);
% U21=zeros(length(u3(1,:))/2,length(u3(1,:))/2);
% U31=zeros(length(u3(1,:))/2,length(u3(1,:))/2);
% for i=1:length(u3(1,:))
%     U11(i)=u3(1,i);
%     U21(i)=u3(2,i);
%     U31(i)=u3(3,i);
% end
% save data11.txt u3 -ascii;
% save data21.txt U21 -ascii;
% save data31.txt U31 -ascii;
% figure,
% plot3(u3(1,:),u3(2,:),u3(3,:),'r.');hold on;
% C1=sqrt(U11.*U11+U21.*U21+U31.*U31);
% % figure,
% surf(U11,U21,U31,C1);
% colorbar