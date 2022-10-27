clear on;
clear variables;

%% doc file am thanh
% input_data : du lieu dau vao
% Fs: Tan so lay mau
[input_data, Fs] = audioread('lab_female.wav');
input_data = input_data./max(abs(input_data)); 
 
%% chia khung cho tin hieu
f_duration = 0.02;  % Thoi gian moi khung
frames = framing(input_data, Fs, f_duration);  
[frame_count,frame_size] = size(frames);
data_size = length(input_data); %Tong so mau tin hieu
data_duration = data_size/Fs; % Thoi gian cua tin hieu

%% Tinh muc nang luong tai moi khung

ste = zeros(0,frame_count);
for i = 1 : frame_count
    ste(i) = sum(frames(i,:).^2);
end
 ste = ste./abs(max(ste));
 
% energy wave (ste wave)
ste_wave = 0;
for j = 1 : length(ste)
    l = length(ste_wave); 
    ste_wave(l : l + frame_size) = ste(j);
end


%% Tim vung co tieng noi
% Nguong xet nang luong ste_threshold 
ste_threshold = 0.01;

% begin_frames luu tru thu tu cac khung ma tai do bat dau tieng noi
% end_frames luu tru thu tu cac khung ma tai do ket thuc tieng noi
[begin_frames, end_frames] = findSpeech(ste, ste_threshold, frames, f_duration);
begin_count = length(begin_frames);
end_count = length(end_frames);

%% Tao du lieu dau ra 
% output_data luu tru thong tin tieng noi
% begin_data luu tru thu tu cac mau tai diem bat dau tieng noi
% end_data luu tru thu tu cac mau tai diem ket thuc tieng noi
output_data = zeros(size(input_data));
begin_data = zeros(size(begin_frames));
end_data = zeros(size(end_frames));
[output_data, begin_data, end_data] = findData(input_data, begin_frames, end_frames, frames);


%% Thiet lap cac duong bien 
begin_data = (begin_data/data_size)*data_duration;
begin_data(2,:) = begin_data(1,:);
end_data = (end_data/data_size)*data_duration;
end_data(2,:) = end_data(1,:);
data_yline = ones(size(begin_data));
data_yline(2,:) = -1;

%% Ve bieu do
% Ve tin hieu dau vao
t= 0 : 1/Fs : length(input_data)/Fs;
t = t(1:end-1);
subplot(3,1,1);
plot(t, input_data);
ylabel('Amplitude (normalized)');
xlabel('Time(s)');
title('Input Data');
legend('Input Data');
grid on;

% Ve bieu do nang nuong 
t_ste = 0: 1/Fs : length(ste_wave)/Fs;
t_ste = t_ste(1:end-1);
subplot(3,1,2)
plot(t_ste, ste_wave, 'r');
ylabel('Amplitude (normalized)');
xlabel('Time(s)');
title('Short-time Energy');
legend('STE');
grid on;

% Ve bieu do phan doan tieng noi
subplot(3,1,3);
plot(t,input_data,'b');
hold on;
plot(0,0,'Color', 'r');
hold on;
plot(1,1,'Color', 'k');
legend('Input Data','Begin', 'End ');
line(begin_data, data_yline, 'Color', 'r', 'HandleVisibility','off', 'LineWidth',2);
line(end_data, data_yline, 'Color', 'k', 'HandleVisibility','off','LineWidth',2);
ylabel('Amplitude (normalized)');
xlabel('Time(s)');
title('Signal Segmentation (Speech)');
grid on;

sound(output_data, Fs);

%% Tim cac khung co tieng noi
function [begin_frames, end_frames] = findSpeech(ste, ste_threshold, frames, f_duration)
[frame_count, frame_size] = size(frames);
begin_frames = []; % Luu tru cac khung bat dau tieng noi
end_frames = []; % Luu tru cac khung ket thuc tieng noi
m=1;
begin_count = 1; % Dem so khung bat dau
end_count= 1; % Dem so khung ket thuc

while m <= frame_count
    % Tim khung m bat dau tieng noi
    if ste(m) < ste_threshold
      m = m + 1;
      continue;
    end 
    
    % Tim khung i ket thuc tieng noi
    i = m;
    while (i <= frame_count) && (ste(i) >= ste_threshold) 
        i = i + 1;
    end
    i = i - 1;
    
    %Kiem tra dieu kien va luu tru cac khung
    if (i - m + 1)*f_duration > f_duration*2 
        begin_frames(begin_count ) = m;
        begin_count  = begin_count  + 1;
        end_frames(end_count) = i;
        end_count = end_count + 1;
    end
    
    % Tiep tuc xet cac khung tiep theo
    m = i + 1;
end 

% Kiem tra lai cac vung chua duoc xac dinh
for i = 1 : length(begin_frames) - 1
    % Neu duoi 0.2, khong duoc xem la khoang lang
    if (begin_frames(1,i+1) - end_frames(1,i))*f_duration < 0.2
        begin_frames(:,i+1) = 0 ;
        end_frames(:,i) = 0;
    end
end

% Xoa cac vi tri khong thoa man
begin_frames( :, ~any(begin_frames,1) ) = [];
end_frames( :, ~any(end_frames,1) ) = [];
end 

%% Tao du lieu dau ra
function [output_data, begin_data, end_data] = findData(input_data, begin_frames, end_frames, frames)

[frame_count, frame_size] = size(frames);
begin_count = length(begin_frames); % Dem so diem bat dau
end_count = length(end_frames); % Dem so diem ket thuc
   
% Xac dinh cac diem bat dau và ket thuc tieng noi
for i=1:begin_count
   begin_data(1,i) = (begin_frames(1,i)-1)*frame_size + 1;
   end_data(1,i) = end_frames(1,i)*frame_size;
end

% Luu tru cac vung co tieng noi
for i = 1 : begin_count
   output_data(begin_data(1,i):end_data(1,i)) = input_data(begin_data(1,i):end_data(1,i));
end 
end

%% Chia khung 
function [frames] = framing(input_data,Fs,f_duration)

frame_size = round(f_duration * Fs);  % so luong mau tai moi khung
data_size = length(input_data);   % so mau cua tin hieu
frame_count = floor(data_size/frame_size); % so khung duoc chia

% Luu tru du lieu vao cac khung
temp = 0;
for i = 1 : frame_count
    frames(i,:) = input_data(temp + 1 : temp + frame_size);
    temp = temp + frame_size;
end

end
