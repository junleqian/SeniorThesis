close all;


clear;
kstr = num2str(0.75);
filename = strcat('_relaxed_Peelnewr', kstr,'m200dur5000p0.8q0.4.mat'); 
load(filename, 'QsizeP');
load(filename, 'lambda');
%% load lambda vector for conveience
load(filename, 'MpP');
load(filename, 'MqP');

lambda_relaxed = lambda;
QsizeP_relaxed = QsizeP;
MpP_relaxed = MpP;
MqP_relaxed = MqP;

filename = strcat('Peelnewr', kstr,'m200dur5000p0.8q0.4.mat'); 
load(filename, 'QsizeP');
load(filename, 'lambda');
%% load lambda vector for conveience
load(filename, 'MpP');
load(filename, 'MqP');

MR = zeros(length(lambda), 2);
MM = zeros(length(lambda), 2);
for k =1:length(lambda)
   MM(k, 1) = lambda(k)/0.8; 
   MR(k, 1) = lambda_relaxed(k)/0.8;
   MM(k, 2)  = QsizeP(k);
   MR(k, 2) = QsizeP_relaxed(k);
   MM(k, 3) = MpP(k); 
   MR(k, 3) = MpP_relaxed(k);
   MM(k, 4) = MqP(k);
   MR(k, 4)= MqP_relaxed(k);
   MM(k, 5) = lambda(k)/(MpP(k)+MqP(k));
   MR(k, 5) = lambda_relaxed(k)/(MpP_relaxed(k)+MqP_relaxed(k));
end
MM = sortrows(MM, 1);
MR = sortrows(MR, 1);
for k=2:5
figure;
plot(MM(:, 1), MM(:, k), 'r','LineWidth',2);
hold on;

plot(MR(:, 1), MR(:, k), 'b','LineWidth', 1);

plot(MM(:, 1), MM(:, k),'cs');

plot(MR(:, 1), MR(:, k), 'm+');
end
