function [Xh] = preproces(ecg,fs,graph)
%X:N*channel,N:time,
%   此处显示详细说明
% ----high-pass filter----
fp1=2;fs1=4;
Fs2=fs/2;
Wp=fp1/Fs2; Ws=fs1/Fs2;
Rp=2; Rs=20;
[n1,Wn1]=buttord(Wp,Ws,Rp,Rs);
[b_hp,a_hp]=butter(n1,Wn1,'high');

% ----low-pass filter----
Fs2=fs/2;
Fpass = 100;   % Passband Frequency
Fstop = 110;  % Stopband Frequency
Wp = Fpass/Fs2;Ws = Fstop/Fs2;
[n2,Wn2] = buttord(Wp,Ws,1,20);
[b_lp,a_lp] = butter(n2,Wn2);


% ----north filter----
d50 = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',48,'HalfPowerFrequency2',52, ...
    'DesignMethod','butter','SampleRate',fs);
d150 = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',98,'HalfPowerFrequency2',102, ...
    'DesignMethod','butter','SampleRate',fs);

d60 = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
    'DesignMethod','butter','SampleRate',fs);
d160 = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',118,'HalfPowerFrequency2',122, ...
    'DesignMethod','butter','SampleRate',fs);

% ----filtering----
NeedFiltering = check_notch(ecg(:,1),50);
if NeedFiltering
    ecgn = filtfilt(d50,ecg);
    ecgn = filtfilt(d150 ,ecgn);
else
    ecgn = ecg;
end

NeedFiltering = check_notch(ecg(:,1),60);
if NeedFiltering
    ecgn = filtfilt(d60,ecgn);
    ecgn = filtfilt(d160 ,ecgn);
else
    ecgn = ecgn;
end

Xb=filtfilt(b_hp,a_hp,ecgn);

Xh=filtfilt(b_lp,a_lp,Xb);

%将功率谱绘制为频率的函数。尽管噪声在基于时间的空间内伪装成信号的频率分量，但傅里叶变换将其显现为功率尖峰。
if graph
    for is=1:4,
        figure, set(gcf,'Color','white');
        pwelch(Xh(:,is),[],[],[],fs);
        title('Detrended ECG Welch spectrum');
        shg
    end
end

end
% fs=1000;
% fs_new=1000;
% 
% fp1=8;fs1=10;
% Fs2=fs_new/2;
% Wp=fp1/Fs2; Ws=fs1/Fs2;
% Rp=2; Rs=20;
% [n1,Wn1]=buttord(Wp,Ws,Rp,Rs);
% [b2,a2]=butter(n1,Wn1,'high');
% 
% % low-pass filter
% Fpass = 100;   % Passband Frequency
% Fstop = 110;  % Stopband Frequency
% Apass = 0.1;    % Passband Ripple (dB)
% Astop = 20;   % Stopband Attenuation (dB)
% h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, fs_new);
% Hlp = design(h, 'butter', ...
%     'MatchExactly', 'stopband', ...
%     'SOSScaleNorm', 'Linf');
% [b_lp,a_lp] = tf(Hlp);
% 
% d = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',48,'HalfPowerFrequency2',52, ...
%     'DesignMethod','butter','SampleRate',fs_new);
% d1 = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',98,'HalfPowerFrequency2',102, ...
%     'DesignMethod','butter','SampleRate',fs_new);
% 
% d = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
%     'DesignMethod','butter','SampleRate',fs_new);
% d1 = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',118,'HalfPowerFrequency2',122, ...
%     'DesignMethod','butter','SampleRate',fs_new);

