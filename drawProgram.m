
clear;
kstr = num2str(0.75);
filename = strcat('Peelnewr', kstr,'m200dur5000p0.8q0.4.mat'); 
load(filename, 'QsizeP');
load(filename, 'lambda');
%% load lambda vector for conveience
load(filename, 'MpP');
load(filename, 'MqP');
MM = zeros( length(lambda), 2);
for k =1:length(lambda)
   MM(k, 1) = lambda(k); 
   MM(k, 2)  = QsizeP(k);
   MM(k, 3) = MpP(k); 
   MM(k, 4) = MqP(k);
   MM(k, 5) = lambda(k)/(MpP(k)+MqP(k));
end
MM = sortrows(MM, 1);
for k=2:5
figure;
plot(MM(:, 1), MM(:, k),'b');
hold;
plot(MM(:, 1), MM(:, k),'r.');
end
