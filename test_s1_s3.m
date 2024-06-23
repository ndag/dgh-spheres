%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Facundo Memoli based on an idea due to Zane Smith  %
% circa 2015.                                                        %
% See https://github.com/ndag/dgh-spheres                            %
% Version of April/2023                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate curve
N = 1000


phi = linspace(0,2*pi,N)';


C = [cos(phi) sin(phi) cos(3*phi) sin(3*phi)]/sqrt(2);


% figure 
% hold on
% scatter3(C(:,1), C(:,2), C(:,3),50,'filled')

% associated dm
Sph1 = [cos(phi) sin(phi)];
d1 = L2_distance(Sph1',Sph1');
d1 = 2*asin(d1/2);

%% generate sphere
N3 = 3000
N3fps = 2000;

Sph3 = randn(N3,4);
dm = L2_distance(Sph3',Sph3');
I  = fps(dm,1,N3fps);
Sph3 = Sph3(I,:);
n3 = (sum(Sph3.^2,2)).^.5;
Sph3 = Sph3 ./ (n3 * ones(1,4));

%scatter3(Sph2(:,1), Sph2(:,2), Sph2(:,3))

d3 = L2_distance(Sph3',Sph3');
d3 = 2 * asin(d3/2);

%% generate map
dm = L2_distance(C',Sph3');
f = zeros(N3fps,1);
for k=1:N3fps
    ik = find(dm(:,k) == min(dm(:,k)));
    f(k) = ik(1);
end

%% evaluate distortion of the map
error = d3-d1(f,f);
max(abs(error(:)))
