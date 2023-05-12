% Part B
close all;
clc;
clear;

N = 8192;

L = 16; n = [-L/2:1:L/2-1];
x = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x(L/2+1) = x(L/2+1) + 80;

X1 = abs(fft(x,L));
X2 = abs(fft(x,N));

figure(1); stem(n,x);
figure(2); stem([0:1:L-1],X1);
figure(3); plot([0:1/(N/2):1],X2(1:N/2+1));
%figure(4); stem([0:1:N-1],X2);

pause

L = 32; n = [-L/2:1:L/2-1];
x = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x(L/2+1) = x(L/2+1) + 80;

X1 = abs(fft(x,L));
X2 = abs(fft(x,N));

figure(4); stem(n,x);
figure(5); stem([0:1:L-1],X1);
figure(6); plot([0:1/(N/2):1],X2(1:N/2+1));
%figure(4); stem([0:1:N-1],X2);

pause

L = 64; n = [-L/2:1:L/2-1];
x = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x(L/2+1) = x(L/2+1) + 80;

X1 = abs(fft(x,L));
X2 = abs(fft(x,N));

figure(7); stem(n,x);
figure(8); stem([0:1:L-1],X1);
figure(9); plot([0:1/(N/2):1],X2(1:N/2+1));
%figure(4); stem([0:1:N-1],X2);

pause

L = 128; n = [-L/2:1:L/2-1];
x = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x(L/2+1) = x(L/2+1) + 80;

X1 = abs(fft(x,L));
X2 = abs(fft(x,N));

figure(10); stem(n,x);
figure(11); stem([0:1:L-1],X1);
figure(12); plot([0:1/(N/2):1],X2(1:N/2+1));
%figure(4); stem([0:1:N-1],X2);

pause

L = 64; n = [-L/2:1:L/2-1];
x1 = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x1(L/2+1) = x1(L/2+1) + 80;

w = 0.54-0.46*cos(2*pi*(n+L/2)/(L-1));
x2 = x1.*w;

X1 = abs(fft(x1,N));
X2 = abs(fft(x2,N));

figure(13); stem(n,x1);
figure(14); stem(n,x2);
figure(15); plot([0:1/(N/2):1],X1(1:N/2+1));
figure(16); plot([0:1/(N/2):1],X2(1:N/2+1));

pause

L = 128; n = [-L/2:1:L/2-1];
x1 = -40*sinc(n/2) + cos(pi*n/16) + cos(pi*n/8) + cos(pi*n/4);
x1(L/2+1) = x1(L/2+1) + 80;

w = 0.54-0.46*cos(2*pi*(n+L/2)/(L-1));
x2 = x1.*w;

X1 = abs(fft(x1,N));
X2 = abs(fft(x2,N));

figure(17); plot([0:1/(N/2):1],X1(1:N/2+1));
figure(18); plot([0:1/(N/2):1],X2(1:N/2+1));


pause


% The specifications of the filter
point1 = 0.2*pi;
point2 = 0.4*pi;
point3 = 0.6*pi;
point4 = 0.8*pi;

dw1 = point2 - point1;
dw2 = point4 - point3;
dw = max(dw1,dw2);

w1 = (point1+point2)/2;
w2 = (point3+point4)/2;

Mhamming = ceil(8*pi/dw);

% 1st way - Matlab functions
Whamming = hamming(Mhamming+1)';


% Calculations
f1 = w1/pi;
f2 = w2/pi;

% Values for Hann and Hamming (since they have the same M)
n2 = 0:1:Mhamming;
a = Mhamming/2;
m=n2-a;


filterhamming = (-f1*sinc(f1*m)+f2*sinc(f2*m)).*Whamming;

fvtool(filterhamming);
%figure(19); %impz(filterhamming);


x3 = conv(x2,filterhamming);
figure(20); stem(n,x2);
figure(21); stem([-L/2:1:L/2+length(filterhamming)-2],x3);
X3 = abs(fft(x3,N));
figure(22); plot([0:1/(N/2):1],X2(1:N/2+1));
figure(23); plot([0:1/(N/2):1],X3(1:N/2+1));

