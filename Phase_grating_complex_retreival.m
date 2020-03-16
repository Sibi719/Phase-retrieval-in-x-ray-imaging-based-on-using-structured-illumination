
clc;
clear all;
close all;
M=512;
N=512;
d1=2*10^-6;
z=0.8;
Tp=(18*10^-6);
it=200;
fm=1/Tp;
lambda=0.10125*10^-9;
noise=randn(N).*sqrt(10^-4);
wn=1;
k=(2*pi)/lambda;
xm= (-M/2:M/2-1)*d1;
ym= (-M/2:M/2-1)*d1;
[Xm,Ym]=meshgrid(xm,ym);
xn= (-N/2:N/2-1)*d1;
yn= (-N/2:N/2-1)*d1;
[Xn,Yn]=meshgrid(xn,yn);

S=zeros(N);
R=(M/2)*d1;
A=abs(Xn)<=R&abs(Yn)<=R;
S(A)=1;

figure
imagesc(xn*10^3,yn*10^3,S);
colormap(gray);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\1.png");

T=imread("Lena.png");
T=rgb2gray(T);
T=double(T);
T=rescale(T,0.3,0.92);
Nr=size(T,1);
Nc=size(T,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
T = padarray(T,[Dr,Dc],0,'both');
figure
imagesc(xn*10^3,yn*10^3,T);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\2.png");

P=imread("dog.jpg");
P=rgb2gray(P);
P=double(P);
P=rescale(P,-0.4*pi,0.3*pi);
Nr=size(P,1);
Nc=size(P,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
P = padarray(P,[Dr,Dc],0,'both');
figure
imagesc(xn*10^3,yn*10^3,P);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\3.png");

Tr=flipdim(T,2);
Pr=flipdim(P,2);


Uo=ones(M);
a=0.65;
n=0.5;
theta=0;

Am=a.*Uo.*( 1+ (n.*sin(2.*pi.*(fm*cosd(theta).*Xm+fm*sind(theta).*Ym))));
Am=rescale(Am,0.85,0.92);
Am = padarray(Am,[Dr,Dc],0,'both');
A= (pi/2).*sin(2.*pi.*(fm*cosd(theta).*Xm+fm*sind(theta).*Ym));
A = padarray(A,[Dr,Dc],0,'both');

theta=0+90;

Amr=a.*Uo.*( 1+ (n.*sin(2.*pi.*(fm*cosd(theta).*Xm+fm*sind(theta).*Ym))));
Amr=rescale(Amr,0.85,0.92);
Amr = padarray(Amr,[Dr,Dc],0,'both');
B= (pi/2).*sin(2.*pi.*(fm*cosd(theta).*Xm+fm*sind(theta).*Ym));
B = padarray(B,[Dr,Dc],0,'both');

MA=exp(1j.*A);
MB=exp(1j.*B);

% figure
% imagesc(xn*10^3,yn*10^3,A);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colorbar
% colormap(gray)
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\4.png");

% figure
% imagesc(xn*10^3,yn*10^3,B);
% colorbar
% colormap(gray)
% xlabel('x(mm)');
% ylabel('y(mm)');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\5.png");



U1A=sqrt(T).*exp(1j.*P).*MA;

% figure
% imagesc(xn*10^3,yn*10^3,abs(U1A));
% colorbar
% colormap(gray)
% xlabel('x(mm)');
% ylabel('y(mm)');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\6.png");


% figure
% imagesc(xn*10^3,yn*10^3,(angle(U1A)));
% colorbar
% colormap(gray)
% xlabel('x(mm)');
% ylabel('y(mm)');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\7.png");


U1B=sqrt(T).*exp(1j.*P).*MB;

% figure
% imagesc(xn*10^3,yn*10^3,abs(U1B));
% colorbar
% colormap(gray)
% xlabel('x(mm)');
% ylabel('y(mm)');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\8.png");


% figure
% imagesc(xn*10^3,yn*10^3,(angle(U1B)));
% colorbar
% colormap(gray)
% xlabel('x(mm)');
% ylabel('y(mm)');
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\9.png");

[Uexp,Iexp]=FresProp(d1,z,lambda,N,sqrt(T).*exp(1j.*P));

[U2A,I2A]=FresProp(d1,z,lambda,N,U1A);%forward_Fresnel(k,z,lambda,Xn,Yn,U1A);

[U2B,I2B]=FresProp(d1,z,lambda,N,U1B);%forward_Fresnel(k,z,lambda,Xn,Yn,U1B);

I2A = I2A + (wn.*noise);
I2B = I2B + (wn.*noise);

% figure
% imagesc(xn*10^3,yn*10^3,I2A);
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\10.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,I2B);
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\11.png");
% 

AA=sqrt(I2A);
BB=sqrt(I2B);


%%GS algothm


U1_est_abs = rand(N).*S;
U1_est_phase= rand(N).*S;
% figure
% imagesc(xn*10^3,yn*10^3,U1_est_abs);
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\18.png");
% figure
% imagesc(xn*10^3,yn*10^3,U1_est_phase);
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\19.png");

U1_est =U1_est_abs.*exp(1j.*U1_est_phase);
beta=0.9;

for i=1:it

 U1_est_A=U1_est.*MA;
 [U2_est_A,~]=FresProp(d1,z,lambda,N,U1_est_A);
%  figure
% imagesc(xn*10^3,yn*10^3,abs(U2_est_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\20.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,angle(U2_est_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\21.png");

 U2_est_A= AA.*exp(1j.*angle(U2_est_A));
 [U1_esti_A,~] = FresProp(d1,-z,lambda,N,U2_est_A);
 
%  figure
% imagesc(xn*10^3,yn*10^3,abs(U1_esti_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% % saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\22.png");
% 
% 
%  figure
% imagesc(xn*10^3,yn*10^3,angle(U1_esti_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\23.png");

 U1_esti_A=U1_esti_A./MA;
%  
%  figure
% imagesc(xn*10^3,yn*10^3,abs(U1_esti_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\24.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,angle(U1_esti_A));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\25.png");
% 

 U1_esti_A(isnan(U1_esti_A))=0;
 U1_esti_A(isinf(U1_esti_A))=0;
 dd= U1_esti_A~=0;
 U1_est_abs = ((abs(U1_esti_A).*dd) + (imcomplement(dd).*(abs(U1_est)-(beta.*abs(U1_esti_A))))).*S;
  U1_est_abs=imgaussfilt(U1_est_abs,1);
%   figure
% imagesc(xn*10^3,yn*10^3,U1_est_abs);
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\26.png");
%   
 U1_est = U1_est_abs .*exp(1j.*(angle(U1_esti_A)+(0.00*rand(N)))).*S;
 
 U1_est_B=U1_est.*MB;
 [U2_est_B,~]=FresProp(d1,z,lambda,N,U1_est_B);
%  figure
% imagesc(xn*10^3,yn*10^3,abs(U2_est_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\27.png");
  
%  figure
% imagesc(xn*10^3,yn*10^3,angle(U2_est_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\28.png");
%   
%  
 
 
 U2_est_B= BB.*exp(1j.*angle(U2_est_B));
 
 
 [U1_esti_B,~] = FresProp(d1,-z,lambda,N,U2_est_B);
 

%  figure
% imagesc(xn*10^3,yn*10^3,abs(U1_esti_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\29.png");

% figure
% imagesc(xn*10^3,yn*10^3,angle(U1_esti_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\30.png");
%   
%  
 
 U1_esti_B=U1_esti_B./MB;
 
%  figure
% imagesc(xn*10^3,yn*10^3,abs(U1_esti_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\31.png");
%   
% figure
% imagesc(xn*10^3,yn*10^3,angle(U1_esti_B));
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\32.png");
  
 U1_esti_B(isnan(U1_esti_B))=0;
 U1_esti_B(isinf(U1_esti_B))=0;
 dd= U1_esti_B~=0;
 U1_est_abs = ((abs(U1_esti_B).*dd) + (imcomplement(dd).*(abs(U1_est)-(beta.*abs(U1_esti_B))))).*S;
 
%  figure
% imagesc(xn*10^3,yn*10^3,U1_est_abs);
% xlabel("x(mm)");
% ylabel("y(mm)");
% colormap(gray)
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\33.png");
%   

  U1_est_abs=imgaussfilt(U1_est_abs,1);
 U1_est = U1_est_abs.*exp(1j.*(angle(U1_esti_B)+(0.0*rand(N)))).*S;
 
  [Uee,Iee]=FresProp(d1,z,lambda,N,U1_est);
 
  err=((Iee - Iexp)./Iexp).^2;
  error(i)= (sum(sum(err)*d1*d1));
end

figure
imagesc(xn*10^3,yn*10^3,(abs(U1_est)).^2);
xlabel("x(mm)");
ylabel("y(mm)");
colormap(gray)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\12.png");


figure
imagesc(xn*10^3,yn*10^3,unwrap(angle(U1_est)));
xlabel("x(mm)");
ylabel("y(mm)");
colormap(gray)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\13.png");

i=1:it;
% figure
% plot(i,error);
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\14.png");

figure
plot(i,log(error));
xlabel("Iterations");
ylabel("Error");
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\15.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2A));
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\16.png");
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2B));
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"D:\Research_implementations\Paper_33\Phase_grating\"+it+"\17.png");

function [U,I]=forward_Fresnel(k,z,lambda,Xn,Yn,Uin)
 U=((exp(1j.*k.*z)/(1j*lambda*z)).*exp(((1j*k)/(2*z)).* (Xn.^2+Yn.^2)).*fftshift(fft2(Uin.*exp(((1j*k)/(2*z)).*(Xn.^2+Yn.^2)))));
 I=abs(U).^2;
end

function [U,I]=back_Fresnel(k,z,lambda,Xn,Yn,Uin)
 U=((1j*lambda*z)./exp(1j.*k.*z)).*exp(((-1j*k)/(2*z)).*(Xn.^2+Yn.^2)).*ifft2(ifftshift(Uin.*exp(((-1j*k)/(2*z)).*(Xn.^2+Yn.^2))));
 I=abs(U).^2;
end

function [Ri,R]=FresProp(dpix,d,lambda,Hsize,Hcrop)
 z=d;

%Spatial frequencies
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Nyquist cut-off for Sampling Hologram
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
 
%Evanescent cut-off 
uev = 1/lambda; %???????
 
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
 
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
 
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
 
%Truncate kernel
H = W.*H;
clear W;
 
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
 
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
Ri = ifft2(RR); 
R=abs(Ri).^2;

end
