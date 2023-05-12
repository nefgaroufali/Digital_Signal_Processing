close all;
clc;
clear;

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

Mhann = ceil(8*pi/dw);
Mhamming = ceil(8*pi/dw);
Mblackman = ceil(12*pi/dw);

% 1st way - Matlab functions
Whann = hann(Mhann+1)';
Whamming = hamming(Mhamming+1)';
Wblackman = blackman(Mblackman+1)';

% step functions with width of each window
Uhann = 0:1:Mhann;
Uhamming = 0:1:Mhamming;
Ublackman = 0:1:Mblackman;

% 2nd way - Manually
W2hann = 0.5-0.5*cos(2*pi*Uhann/Mhann);
W2hamming = 0.54-0.46*cos(2*pi*Uhamming/Mhamming);
W2blackman = 0.42-0.5*cos(2*pi*Ublackman/Mblackman) + 0.08*cos(4*pi*Ublackman/Mblackman);

% Figures of the 1st Windows
figure(1); stem(Uhann,Whann);
figure(2); stem(Uhamming,Whamming);
figure(3); stem(Ublackman,Wblackman);

% Figures of the 2nd Windows
figure(4); stem(Uhann,W2hann);
figure(5); stem(Uhamming,W2hamming);
figure(6); stem(Ublackman,W2blackman);


% Calculations
f1 = w1/pi;
f2 = w2/pi;

% Values for Hann and Hamming (since they have the same M)
n = 0:1:Mhann;
a = Mhann/2;
m=n-a;

%Values for Blackman
n2 = 0:1:Mblackman;
a2 = Mblackman/2;
m2=n2-a2;

% 1st way Filter
synartisi = (-f1*sinc(f1*m)+f2*sinc(f2*m));
filterhann1 = (-f1*sinc(f1*m)+f2*sinc(f2*m)).*Whann;
filterhamming1 = (-f1*sinc(f1*m)+f2*sinc(f2*m)).*Whamming;
filterblackman1 = (-f1*sinc(f1*m2)+f2*sinc(f2*m2)).*Wblackman;

% 2nd way Filter
hb = ((-sin(w1*(n-a)) + sin(w2*(n-a)))./(pi*(n-a)));
hb2 = ((-sin(w1*(n2-a2)) + sin(w2*(n2-a2)))./(pi*(n2-a2)));
filterhann2 = hb.*Whann;
filterhamming2 = hb.*Whamming;
filterblackman2 = hb2.*Wblackman;

% 3rd way Filter
filterhann3 = fir1(Mhann,[f1 f2]);
filterhamming3 = fir1(Mhamming,[f1 f2]);
filterblackman3 = fir1(Mblackman,[f1 f2]);

% All filter graphs
figure(7); stem(n,filterhann1);
figure(8); stem(n,filterhann2);
figure(9); stem(n,filterhann3);

figure(10); stem(n,filterhamming1);
figure(11); stem(n,filterhamming2);
figure(12); stem(n,filterhamming3);

figure(13); stem(n2,filterblackman1);
figure(14); stem(n2,filterblackman2);
figure(15); stem(n2,filterblackman3);


% Characteristics of the filters
fvtool(filterhann1);
figure(16); impz(filterhann1);
figure(17); stepz(filterhann1);
figure(18); grpdelay(filterhann1);
figure(19); freqz(filterhann1);

% Characteristics of the filters
fvtool(filterhamming1);
figure(20); impz(filterhamming1);
figure(21); stepz(filterhamming1);
figure(22); grpdelay(filterhamming1);
figure(23); freqz(filterhamming1);


% Characteristics of the filters
fvtool(filterblackman1);
figure(24); impz(filterblackman1);
figure(25); stepz(filterblackman1);
figure(26); grpdelay(filterblackman1);
figure(27); freqz(filterblackman1);


% DFT of filter
X1old = abs(fft(filterhann1,8192));
X2old = abs(fft(filterhamming1,8192));
X3old = abs(fft(filterblackman1,8192));

X1(1:8192/2+1) = X1old(1:8192/2+1);
X2(1:8192/2+1) = X2old(1:8192/2+1);
X3(1:8192/2+1) = X3old(1:8192/2+1);

figure(28); plot([0:1/(8192/2):1],X1);
figure(29); plot([0:1/(8192/2):1],X2);
figure(30); plot([0:1/(8192/2):1],X3);


% Error function
for i = 1:1:4097
    branchhann1(i) = X1(i);
    branchhann2(i) = 1-X1(i);
    branchhann3(i) = X1(i);
    error_function1(i) = 0;
end

branchhann1 = branchhann1(1:818);
branchhann2 = branchhann2(1638:2456);
branchhann3 = branchhann3(3276:4097);

error_function1(1:818) = branchhann1(1:818);
error_function1(1638:2456) = branchhann2(1:819);
error_function1(3276:4097) = branchhann3(1:822);

% Error function Graph
figure(31); plot([0:1/(8192/2):1],error_function1);

% Error function
for i = 1:1:4097
    branchhamming1(i) = X2(i);
    branchhamming2(i) = 1-X2(i);
    branchhamming3(i) = X2(i);
    error_function2(i) = 0;
end

branchhamming1 = branchhamming1(1:818);
branchhamming2 = branchhamming2(1638:2456);
branchhamming3 = branchhamming3(3276:4097);

error_function2(1:818) = branchhamming1(1:818);
error_function2(1638:2456) = branchhamming2(1:819);
error_function2(3276:4097) = branchhamming3(1:822);

% Error function Graph
figure(32); plot([0:1/(8192/2):1],error_function2);

% Error function
for i = 1:1:4097
    branchblackman1(i) = X3(i);
    branchblackman2(i) = 1-X3(i);
    branchblackman3(i) = X3(i);
    error_function3(i) = 0;
end

branchblackman1 = branchblackman1(1:818);
branchblackman2 = branchblackman2(1638:2456);
branchblackman3 = branchblackman3(3276:4097);

error_function3(1:818) = branchblackman1(1:818);
error_function3(1638:2456) = branchblackman2(1:819);
error_function3(3276:4097) = branchblackman3(1:822);

% Error function Graph
figure(33); plot([0:1/(8192/2):1],error_function3);