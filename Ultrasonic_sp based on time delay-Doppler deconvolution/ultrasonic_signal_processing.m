%%%%%%%%
% Ultrasonic signal processing method based on time delay Doppler
% deconvolution
% Ultrasonic signal modeling:
% ultrasonic signal can be modeled by Gaussian modulated sine wave
%%%%%%%
clc
clear
close all
a = 0.92;
amp = [0.95,0.85,0.75];
% a1 = 0.95; a2 = 0.85; a3 = 0.75;
b = 1.4062e15;  % the bandwidth factor b represents the square of the bandwidth（37.5MHZ）^2
tao = 0.5e-6;   % delay time
de_t = [0.4e-6,0.425e-6,0.45e-6]; % target delay vector
I = 3;          % echo number
% tao1 = 0.25e-6; tao2 = 0.3e-6; tao3 = 0.35e-6;
fc = 25e6;      % center frequency
b1 = 2;
fs = 1e9;       % sampling frequency
T = 1e-6;
t = -T:1/fs:T-1/fs;
tt = 0:1/fs:T-1/fs;
k_it = [10 20 40 60];
% PRF = 5e5;
c = 1500;              % speed of sound
v = [0.1,0.06,0.08];   % movement speed
V_dop = 2*fc.*v./c;
% doppler1 = exp(1i*2*pi*V.*tt);
s1 = exp(-b*t.^2).*cos(2*pi*fc*t);
s1_dt = [zeros(1,round(tao*fs)),s1(1:(length(s1)-round(tao*fs)))];
r_s2 = zeros(I,round(max(de_t)*fs)+length(tt));
% x1 = complex(a*exp(-b*(t-tao).^2).*cos(2*pi*fc*(t-tao)),0).*dopp;

% N = length(t);
% f = fs*(0:N/2)/N;
% s1f = fft(real(s1),N);
% x1f = fft(real(x1),N);
% s1f = abs(s1f./N); s1F = s1f(1:N/2+1);
% x1f = abs(x1f./N); x1F = x1f(1:N/2+1);
% figure
% subplot(211);
% plot(f,s1F); xlabel('频率(HZ)'); ylabel('幅度'); title('发射信号频谱');
% subplot(212);
% plot(f,x1F); xlabel('频率(HZ)'); ylabel('幅度'); title('回波信号频谱');
% [sv,poss] = max(s1F); possf = f(poss);
% [xv,posv] = max(x1F); posvf = f(posv);
s2 = exp(-b*tt.^2).*cos(2*pi*fc*tt);
% triple echo signals
for i = 1:I
    r_s2(i,round(de_t(i)*fs)+1:round(de_t(i)*fs)+length(s2)) = amp(i).*s2.*exp(1j*2*pi*V_dop(i).*tt);
end
r_echo = sum(r_s2,1);
t_p = (0:1:length(r_echo)-1)./fs;
% x2=complex(a1*exp(-b*(tt-tao1).^2).*cos(2*pi*fc*(tt-tao1)+b1),0).*doppler1+...
%    complex(a2*exp(-b*(tt-tao2).^2).*cos(2*pi*fc*(tt-tao2)+b1),0).*doppler1+...
%    complex(a3*exp(-b*(tt-tao3).^2).*cos(2*pi*fc*(tt-tao3)+b1),0).*doppler1;
% add noise with SNR 10dB
r_echon = awgn(r_echo,0,'measured','db');
xzratio = snr(r_echo,r_echon-r_echo);

%%
maxDelay1 = 0.2e-6;
maxDoppler1 = 6700;
maxDelay2 = 0.7e-6;
tstart = 0.2e-6;
maxDoppler2 = 6700;
[a_fmag,delay_a,dopp_a] = computeAmbiguityFunction(s2,fs,maxDoppler1,maxDelay1);
% pay attention to the order of two signals in the cross-ambiguity function
% the first signal is the transmission signal
% the second signal is the echo signal
[c_fmag,delay_c,dopp_c] = computeCrossAF(s2, r_echon, fs, maxDoppler2, maxDelay2,tstart);
% inversion filling for drawing
a_fmag_d = [fliplr(a_fmag'),a_fmag'];
delay_a_d = [-fliplr(delay_a),delay_a];
a_nom_max = 0;
for i = 1:length(a_fmag_d(:,1))
    a_nom_max = max(a_nom_max,max(a_fmag_d(i,:)));
end
a_fmag_d_nom = a_fmag_d./a_nom_max;
% % the cross-ambiguity function does not need to be filled, but only normalized
c_nom_max = 0;
for i = 1:length(c_fmag(:,1))
    c_nom_max = max(c_nom_max,max(c_fmag(i,:)));
end
c_fmag_nom = c_fmag'./c_nom_max;
%%
% R-L deconvolution
row_c = size(c_fmag_nom,1);
col_c = size(c_fmag_nom,2);
fsmd = zeros(row_c,col_c,4);
fsmd_nom = fsmd;
new_psf = imresize(a_fmag_d_nom,[row_c col_c]);
for i = 1:4
    fsmd(:,:,i) = deconvlucy(c_fmag_nom,new_psf,k_it(i));
    fsmd_nom_max = 0;
    for j = 1:row_c
        fsmd_nom_max = max(fsmd_nom_max,max(fsmd(j,:,i)));
    end
    fsmd_nom(:,:,i) = fsmd(:,:,i)./fsmd_nom_max;
end
%%
% Auto-Cross-AF plot
figure;
imagesc(delay_a_d*c/2, dopp_a*c/(2*fc), abs(a_fmag_d_nom));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('自模糊度函数');
colormap jet
colorbar;
figure;
imagesc(delay_c*c/2, dopp_c*c/(2*fc), abs(c_fmag_nom));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('互模糊度函数');
colormap jet
colorbar;
%%
%R-L deconvolution
figure
imagesc(delay_c*c/2, dopp_c*c/(2*fc), abs(fsmd_nom(:,:,1)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 10 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*fc), abs(fsmd_nom(:,:,2)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 20 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*fc), abs(fsmd_nom(:,:,3)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 40 反射性密度函数');
colormap jet
colorbar;
figure
imagesc(delay_c*c/2, dopp_c*c/(2*fc), abs(fsmd_nom(:,:,4)));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title('k = 60 反射性密度函数');
colormap jet
colorbar;
%%
figure
subplot(211)
plot(t,s1); xlabel('时间'); ylabel('幅度'); title('发射信号');
subplot(212)
plot(t,s1_dt); xlabel('时间'); ylabel('幅度'); title('回波信号');
figure
plot(tt,s2); xlabel('时间'); ylabel('幅度'); title('发射信号');
figure
subplot(211)
plot(t_p,real(r_echo)); xlabel('时间'); ylabel('幅度'); title('回波信号（无噪声）');
subplot(212)
plot(t_p,real(r_echon)); xlabel('时间'); ylabel('幅度'); title('回波信号（带噪）');