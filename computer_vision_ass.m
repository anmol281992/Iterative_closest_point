%% PROBLEM 4: The ICP algorithm for 3D DATA WITH NOISE
clear all;
close
%Model data and source data assigning for 3d Noise data
a = load('3D_Cat_Noise.mat');
mod = a.model;
source = a.source;
figure(1);
hold on 
grid on 
plot3(mod(1,:),mod(2,:),mod(3,:),'LineStyle', ':','Color','r');
plot3(source(1,:),source(2,:),source(3,:),'LineStyle', ':','Color','b');
legend('model','source');
hold off
%% 7 iteration for the algorithm - ICP 
for k = 1:1:7
[m,n] = size(mod);
[ms,ns] = size(source);
v = zeros(1,ns);
diff = zeros(1,ns);
% closest point algorithm 
for j = 1:1:n
        mval = 9e99;
        val =sqrt(sum((source - repmat(mod(:,j),1,ns)).^2));
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
cov = source* modchanged' - 3*centroidsource*centroidmod';
[U,~,V]=svd(cov);
Ri=V*U';% rotation matrix
T = centroidmod - Ri*centroidsource;%calculating translation matrix
Changedpossource = Ri*source + repmat(T,1,ns);%changing the source vaulue
source = Changedpossource;
end
% plotting the data
figure(2)
hold on 
grid on 
plot3(mod(1,:),mod(2,:),mod(3,:),'LineStyle', ':','Color','r');
plot3(source(1,:),source(2,:),source(3,:),'LineStyle', ':','Color','b');
legend('model','source');
hold off







