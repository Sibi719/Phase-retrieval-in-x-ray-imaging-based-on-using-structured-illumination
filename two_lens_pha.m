clc;
clear all;
close all;
M=512;
N=512;
d1=(10*10^-3)/M;
z=0.02;
Tp=(260*10^-9);
it=600;
fm=1/Tp;
lambda=650*10^-9;
noise=randn(N).*sqrt(10^-4);
wn=0;
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
colormap(jet);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(jet);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"1.png");



R1=3*10^-3;
R2=1.8*10^-3;
alpha1=0.001%0.0000002;
beta1=0.00012;%0.00000001;
alpha2=alpha1/1.1;
beta2=beta1/2;
cx1=0;
cy1=0;
wn=0;
cx2=-R1/1.8;
cy2=0*10^-3;
L1= real(sqrt(R1^2-((Xm-cx1).^2 + (Ym-cy1).^2)));
L2= real(sqrt(R2^2-((Xm-cx2).^2 + (Ym-cy2).^2)));
phase1 = k*alpha1*L1/1.4;
phase2 = k*alpha2*L2;
P= phase1+phase2;
T1 = k*beta1*L1/8;
T2 = k*beta2*L2;
T= exp(-(T1+T2)).^2;
%T=rescale(T,0.3,0.92);
Nr=size(T,1);
Nc=size(T,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
T = padarray(T,[Dr,Dc],0,'both');

figure
imagesc(xn*10^3,yn*10^3,T);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(jet);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"2.png");

P=rescale(P,-0.4*pi,0.5*pi);
Nr=size(P,1);
Nc=size(P,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
P = padarray(P,[Dr,Dc],0,'both');
figure
imagesc(xn*10^3,yn*10^3,P);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(jet);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"3.png");

n=0.5;
theta1=0;
theta2=90;


A= pi*0.5*(square(2.*pi.*(fm*cosd(theta1).*Xm+fm*sind(theta1).*Ym))+1);
A = padarray(A,[Dr,Dc],0,'both');

B= pi*0.5*(square(2.*pi.*(fm*cosd(theta2).*Xm+fm*sind(theta2).*Ym))+1);
B = padarray(B,[Dr,Dc],0,'both');

MA=exp(1j.*A);
MB=exp(1j.*B);

figure
imagesc(xn*10^3,yn*10^3,A);
xlabel('x(mm)');
ylabel('y(mm)');
colorbar
colormap(jet)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"4.png");

figure
imagesc(xn*10^3,yn*10^3,B);
colorbar
colormap(jet)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"5.png");



U1A=sqrt(T).*exp(1j.*P).*MA;

% figure
% imagesc(xn*10^3,yn*10^3,abs(U1A));
% colorbar
% colormap(jet)
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
% saveas(gcf,"6.png");


% figure
% imagesc(xn*10^3,yn*10^3,(angle(U1A)));
% colorbar
% colormap(jet)
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
% saveas(gcf,"7.png");


U1B=sqrt(T).*exp(1j.*P).*MB;



[Uexp,Iexp]=FresProp(d1,z,lambda,N,sqrt(T).*exp(1j.*P));

[U2A,I2A]=FresProp(d1,z,lambda,N,U1A);%forward_Fresnel(k,z,lambda,Xn,Yn,U1A);

[U2B,I2B]=FresProp(d1,z,lambda,N,U1B);%forward_Fresnel(k,z,lambda,Xn,Yn,U1B);

I2A = I2A + (wn.*noise);
I2B = I2B + (wn.*noise);


AA=sqrt(I2A);
BB=sqrt(I2B);


%%GS algothm


U1_est_abs = rand(N).*S;
U1_est_phase= rand(N).*S;

U1_est =U1_est_abs.*exp(1j.*U1_est_phase);
beta=0.9;

for i=1:it

 U1_est_A=U1_est.*MA;
 [U2_est_A,~]=FresProp(d1,z,lambda,N,U1_est_A);


 U2_est_A= AA.*exp(1j.*angle(U2_est_A));
 [U1_esti_A,~] = FresProp(d1,-z,lambda,N,U2_est_A);
 

 U1_esti_A=U1_esti_A./MA;

 U1_esti_A(isnan(U1_esti_A))=0;
 U1_esti_A(isinf(U1_esti_A))=0;
 dd= U1_esti_A~=0;
 U1_est_abs = ((abs(U1_esti_A).*dd) + (imcomplement(dd).*(abs(U1_est)-(beta.*abs(U1_esti_A))))).*S;
 U1_est_abs=imgaussfilt(U1_est_abs,0.1);

 U1_est = U1_est_abs .*exp(1j.*(angle(U1_esti_A)+(0.00*rand(N)))).*S;
 
 U1_est_B=U1_est.*MB;
 [U2_est_B,~]=FresProp(d1,z,lambda,N,U1_est_B);

 
 
 U2_est_B= BB.*exp(1j.*angle(U2_est_B));
 
 
 [U1_esti_B,~] = FresProp(d1,-z,lambda,N,U2_est_B);
 

%  
 
 U1_esti_B=U1_esti_B./MB;

 U1_esti_B(isnan(U1_esti_B))=0;
 U1_esti_B(isinf(U1_esti_B))=0;
 dd= U1_esti_B~=0;
 U1_est_abs = ((abs(U1_esti_B).*dd) + (imcomplement(dd).*(abs(U1_est)-(beta.*abs(U1_esti_B))))).*S;
 

  U1_est_abs=imgaussfilt(U1_est_abs,0.1);
 U1_est = U1_est_abs.*exp(1j.*(angle(U1_esti_B)+(0.0*rand(N)))).*S;
 
  [Uee,Iee]=FresProp(d1,z,lambda,N,U1_est);
 
  err=((Iee - Iexp)./Iexp).^2;
  error(i)= (sum(sum(err)*d1*d1));
end

figure
imagesc(xn*10^3,yn*10^3,(abs(U1_est)).^2);
xlabel("x(mm)");
ylabel("y(mm)");
colormap(jet)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"12.png");


figure
imagesc(xn*10^3,yn*10^3,unwrap(angle(U1_est)));
xlabel("x(mm)");
ylabel("y(mm)");
colormap(jet)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"13.png");

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
% saveas(gcf,"14.png");

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
saveas(gcf,"15.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2A));
% colormap(jet);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(jet);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"16.png");
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2B));
% colormap(jet);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(jet);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"17.png");

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