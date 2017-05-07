%% Problem2: The ICP algorithm for 2D LINE DATA WITH NOISE

clear all;
close all;
% Model data and source data assigning for 2d data noise
a = load('2D_line_noise.mat');
mod = a.model;
source = a.source;
figure(1)
hold on 
grid on 
plot(mod(1,:),mod(2,:),'LineStyle', ':','Color','r');
plot(source(1,:),source(2,:),'LineStyle', ':','Color','b');
legend('model','source');
hold off
%% k iteration for the algorithm - ICP 
for k = 1:1:20
[m,n] = size(mod);
[ms,ns] = size(source);
v = zeros(1,ns);
diff = zeros(1,ns);
% closest point algorithm 
for j = 1:1:ns
        mval = 9e99;
        val =sqrt(sum((mod - repmat(source(:,j),1,n)).^2));
        if val<=mval
            [minim,v(j)] = min(val);
        end
end
 modchanged = mod(:,v);
% application of Principal component analysis for finding the rotation
% matrix
 centroidmod = mean(modchanged,2);
centroidsource = mean(source,2);
%Cov(x) = E(xy) - 3*E(x)*E(y)
cov = source* modchanged' - 3*centroidsource*centroidmod';%Covariance Matrix Calculation
[U,~,V]=svd(cov);%Singular Value decomposition 
Ri=V*U';%Rotation Matrix
T = centroidmod - Ri*centroidsource;%translation vector
Changedpossource = Ri*source + repmat(T,1,ns);%position of the source changed
source = Changedpossource;
end
%% plotting the data 
figure(2)
hold on 
grid on 
plot(mod(1,:),mod(2,:),'LineStyle', ':','Color','r');
plot(source(1,:),source(2,:),'LineStyle', ':','Color','b');
legend('model','source');
hold off

