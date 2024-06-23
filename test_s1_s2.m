%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Facundo Memoli based on an idea due to Zane Smith  %
% circa 2015.                                                        %
% See https://github.com/ndag/dgh-spheres                            %
% Version of April/2023                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000
N2 = 3000


%% generate curve

phi = linspace(0,2*pi,N)';


C = [cos(phi) sin(phi) cos(3*phi)];
nor = sum(C.^2,2).^.5;
C =  C./(nor*[1 1 1]);

% figure 
% hold on
% scatter3(C(:,1), C(:,2), C(:,3),50,'filled')

% associated dm
Sph1 = [cos(phi) sin(phi)];
d1 = L2_distance(Sph1',Sph1');
d1 = 2*asin(d1/2);

%% generate sphere
Sph2 = randn(N2,3);
n2 = (sum(Sph2.^2,2)).^.5;
Sph2 = Sph2 ./ (n2 * ones(1,3));

%scatter3(Sph2(:,1), Sph2(:,2), Sph2(:,3))

d2 = L2_distance(Sph2',Sph2');
d2 = 2 * asin(d2/2);

%% generate map
dm = L2_distance(C',Sph2');
f = zeros(N2,1);
for k=1:N2
    ik = find(dm(:,k) == min(dm(:,k)));
    f(k) = ik(1);
end

%% evaluate distortion of the map
error = d2-d1(f,f);
max(abs(error(:)))
