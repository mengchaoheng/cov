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
title('最大力矩顶点拟合成的面');   
xlabel('X轴（N*m）');        
ylabel('Y轴（N*m）');        
zlabel('Z轴（N*m）');       
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
cc=colorbar;
cc.Label.String = '最大力矩幅值（N*m）';

Color1=sqrt(X.*X+Y.*Y+Z.*Z);
figure,
surf(X,Y,Z,Color1);
% shading interp;
cc=colorbar;
cc.Label.String = '最大力矩幅值（N*m）';
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

A1=1/(0.3043*0.39088/(0.3043-0.1521));B1=1/(0.3043*0.39088/(0.3043-0.1521));C1=1/0.3043;D1=-1;% 法向量就是（A1，B1，C1）

A2=1;B2=0;C2=0;D2=-0.39088;
A3=0;B3=1;C3=0;D3=-0.39088;
a1=[0;0;0.3043];b1=[0.39088;0;0.1521];c1=[0.39088;0.39088;0];d1=[0;0.39088;0.1521];
a2=[0.39088;0;0];b2=[0.39088;0;0.1521];c2=[0.39088;0.39088;0];%
a3=[0;0.39088;0];b3=[0;0.39088;0.1521];c3=[0.39088;0.39088;0];
a4=a1;b4=[0.39088;0;0.3043];c4=[0.39088;0.39088;0.3043];d4=[0;0.39088;0.3043];
e4=[0.39088;0;0];f4=[0.39088;0.39088;0];g4=[0;0.39088;0];

V=[0.5;0.3;0.4];
V1(1)=Constrain(V(1),0,0.39088);
V1(2)=Constrain(V(2),0,0.39088);
V1(3)=Constrain(V(3),0,0.3043);

X1 = [a1(1); b1(1); c1(1); d1(1)];
Y1 = [a1(2); b1(2); c1(2); d1(2)];
Z1 = [a1(3); b1(3); c1(3); d1(3)];
Color1=sqrt(X1.*X1+Y1.*Y1+Z1.*Z1);
figure
h1=fill3(X1,Y1,Z1,Color1);
set(h1,'FaceAlpha',0.5);
hold on
X2 = [a2(1); b2(1); c2(1)];
Y2 = [a2(2); b2(2); c2(2)];
Z2 = [a2(3); b2(3); c2(3)];
Color2=sqrt(X2.*X2+Y2.*Y2+Z2.*Z2);
h2=fill3(X2,Y2,Z2,Color2);
set(h2,'FaceAlpha',0.5);
hold on
X3 = [a3(1); b3(1); c3(1)];
Y3 = [a3(2); b3(2); c3(2)];
Z3 = [a3(3); b3(3); c3(3)];
Color3=sqrt(X3.*X3+Y3.*Y3+Z3.*Z3);
h3=fill3(X3,Y3,Z3,Color3);
set(h3,'FaceAlpha',0.5);
cc=colorbar;
cc.Label.String = '最大力矩幅值（N*m）';
hold on
X4 = [a4(1); b4(1); c4(1);d4(1)];
Y4 = [a4(2); b4(2); c4(2);d4(2)];
Z4 = [a4(3); b4(3); c4(3);d4(3)];
h4=fill3(X4,Y4,Z4,[.5 .5 .5]);
set(h4,'FaceAlpha',0.1);
hold on
X5 = [ b4(1); e4(1);f4(1);c4(1)];
Y5 = [ b4(2); e4(2);f4(2);c4(2)];
Z5 = [ b4(3); e4(3);f4(3);c4(3)];
h5=fill3(X5,Y5,Z5,[.5 .5 .5]);
set(h5,'FaceAlpha',0.1);
hold on
X6 = [ c4(1);d4(1);g4(1);f4(1)];
Y6 = [c4(2);d4(2);g4(2);f4(2)];
Z6 = [c4(3);d4(3);g4(3);f4(3)];
h6=fill3(X6,Y6,Z6,[.5 .5 .5]);
set(h6,'FaceAlpha',0.1);
hold on
h7=quiver3(0,0,0,V(1),V(2),V(3),1,'-.g','LineWidth',1);
set(h7,'maxheadsize',0.5);hold on;
h8=quiver3(0,0,0,V1(1),V1(2),V1(3),1,'--m','LineWidth',1);
set(h8,'maxheadsize',0.5);hold on;
M = fun1(A1,B1,C1,D1,V1);
h9=quiver3(0,0,0,M(1),M(2),M(3),1,'r','LineWidth',1);
set(h9,'maxheadsize',0.7);hold on;
h10=plot3(M(1),M(2),M(3),'k.','MarkerSize',12);
title('第一卦限最大力矩顶点拟合成的面');   
axis([0 0.6 0 0.6 0 0.45]);grid on;
legend([h7 h8 h9 h10],'原始力矩','中间力矩','限幅后力矩','力矩与允许域边界交点');
xlabel('X轴（N*m）');        
ylabel('Y轴（N*m）');        
zlabel('Z轴（N*m）');  



% p_start=zeros(length(u3(:,)),3);%向量起点
% x=p_start(:,1);
% y=p_start(:,2);
% z=p_start(:,3);
% u=u3(:,1);
% v=u3(:,2);
% w=u3(:,3);
%  figure
% quiver3(x,y,z,u,v,w,1);

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