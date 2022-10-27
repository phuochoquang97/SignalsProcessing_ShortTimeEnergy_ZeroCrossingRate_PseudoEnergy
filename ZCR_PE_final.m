clear variables
close all
%% Doc tin hieu vao
[inputdata, Fs] = audioread('studio_female.wav');
%phat tin hieu vao
sound(inputdata, Fs);
%chuan hoa tin hieu vao
inputdata = inputdata /max(abs(inputdata));

%% Tinh thoi gian
t = [0 : 1 / Fs : length(inputdata) / Fs]; % time in sec 
t = t(1:end - 1);

%% Plot inputdata
subplot(3, 1, 1);
plot(t, inputdata);
xlabel('t(s)');
ylabel('Ampiltude (normalized)');
title('Input Data and Silence Segmentation Prediction');
grid on;

%% chia khung
% do dai moi khung la 20ms
f_duration = 0.02;
frames = framing(inputdata, Fs, f_duration);
%row = so hang ma tran frames = tong so khung
%col = so cot mat tran frames = kich thuoc moi khung
[row, col] = size(frames);

%% tinh so lan bang qua 0 moi khung
ZCRframe = calculateZCRframe(frames,row);

%% tinh gia nang luong cua tin hieu
En = calculateEn(frames,row);
En = En / max(En);
    L = 1;
    for j = 1 : length(En)
        Enwave(L : L + col - 1) = En(j);
        L = length(Enwave) + 1;
    end

%% tinh ti le bang qua 0 cua tin hieu
%tinh ti le bang qua 0 cua cac frame
ZCRrate = ZCRframe / col;
% chuan hoa ZCRrate
ZCRrate = ZCRrate / abs(max(ZCRrate));
%tinh ti le bang qua 0 cua tat ca cac mau tin hieu vao
ZCRwave = calculateZCRwave(ZCRrate, col);

%% calculate time
t1 = 0 : 1 / Fs : length(ZCRwave) / Fs;
t1 = t1(1 : end - 1);

%% plot ZCRwave and Pseudo-Energy
subplot(3, 1, 2);
plot(t1, ZCRwave, 'r');
hold on;
plot(t1,Enwave,'b');
title('Zero-Crossing Rate and Pseudo-Energy');
grid on;
xlabel('t(s)');
ylabel('Ampiltude (normalized)');
legend('ZCR','PE');

subplot(3, 1, 3);
plot(t, inputdata);
grid on;
xlabel('t(s)');
ylabel('Ampiltude (normalized)');
title('Signal Segmentation (Silence)');

%% Tim bien cua khoang lang
m = max(ZCRrate)-0.2;
for i=2 : length(ZCRrate)
    if En(i) < 0.09 || ZCRrate(i) > m
            % markb danh dau vi tri bat dau khoang lang
            markb(i) = i;
            if markb(i) ~= 0 && markb(i - 1) == 0
                x1 = [(markb(i) - 1)*col/Fs (markb(i) - 1)*col/Fs];
                y1 = [-1 1];
            end
            
            while En(i) < 0.1 && i < length(ZCRrate)
                i=i+1;
            end
            
            %marke danh dau vi tri ket thuc khoang lang
            marke(i)=i;
            x2=[(marke(i)-1)*col/Fs (marke(i)-1)*col/Fs];
            y2=[-1 1];
            
            %kiem tra neu khoang lang lon hon 200ms thi ke duong bien doc
            if (x2(1)-x1(1)) >= 0.2
                line(x1, y1, 'color', 'r', 'LineWidth',2);
                line(x2,y2, 'color', 'k', 'LineWidth', 2);
            end
    end
end
%chu thich
legend('Inputdata','Begin','End');

%% Chia tin hieu thanh cac khung
%inputdata: tin hieu dau vao
function [frames] = framing(inputdata,Fs,f_duration)
% inputdata: tin hieu vao
% Fs: tan so
% f_duration: do dai moi khung, tinh theo giay
% frames: ma tran ket qua, moi hang la mot khung

% tinh va lam tron kich thuoc moi khung
Framesize = round(f_duration * Fs);
% tinh so khung
numberFrame = floor(length(inputdata)/Framesize);

% chia khung
temp = 0;
for i = 1 : numberFrame
    frames(i,:) = inputdata(temp + 1 : temp + Framesize);
    temp = temp + Framesize;
end
end

%% Tinh so lan bang qua 0 moi khung
%ZCRframe: ket qua tra ve luu o ma tran ZCRframe
%frames: ma tran luu gia tri cua tung khung
%row: so hang cua ma tran frames
function [ZCRframe]=calculateZCRframe(frames, row)
for i = 1 : row
    currentFrame = frames(i, :);
    ZCRframe(i) = 0;
    for k = 1 : length(currentFrame)
        if currentFrame(k) > 0 
            index(k) = 1;
        else
            index(k) = -1;
        end
    end
    for k = 2 : length(currentFrame)
        ZCRframe(i) = ZCRframe(i) + abs(index(k) - index(k - 1));
    end
    ZCRframe(i) = ZCRframe(i) / 2;
end
end
%% Tinh En moi khung
%En: ket qua tra ve luu o ma tran En
%frames: ma tran luu gia tri cua tung khung
%row: so hang cua ma tran frames
function [En] = calculateEn(frames, row)
for i = 1 : row
    currentFrame = frames(i, :);
    En(i) = 0;
    for k = 1 : length(currentFrame)
         En(i) = En(i) + abs(currentFrame(k));
    end
end
end
%% Tinh ZCR cua tin hieu vao
%ZCRate: ma tran luu ZCRrate moi khung
%ZCRframe: ma tran luu so lan bang qua 0 cua moi khung
%col: so cot cua ma tran frames
%ZCRwave: ket qua tra ve la ma tran luu ZCRrate cua tat ca mau tin hieu ban dau
function [ZCRwave] = calculateZCRwave(ZCRrate, col)
    ZCRwave = 0;
    L = 1;
    for j = 1 : length(ZCRrate)
        ZCRwave(L : L + col - 1) = ZCRrate(j);
        L = length(ZCRwave) + 1;
    end
end


    





