%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Facundo Memoli based on an idea due to Zane Smith  %
% circa 2015.                                                        %
% See https://github.com/ndag/dgh-spheres                            %
% Version of April/2024                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S^1 --->> S^n
N = 5000   % Number of points on curve
Nn = 5000  % Number of points on S^n
k = 2     % This code is for n =  2*k - 1 (S^1 --->> S^n). Ex: k=2 --> n=3; k=3 --> n=5


% generate curve
phi = linspace(0,2*pi,N)';

C = [];
for i=1:k
    C = [C [cos((2*i-1)*phi) sin((2*i-1)*phi)]];
end

C = C/sqrt(k);



% compute associated distance matrix
Sph1 = [cos(phi) sin(phi)];
d1 = L2_distance(Sph1',Sph1');
d1 = 2*asin(d1/2);

% generate sphere S^n
Sphn = randn(Nn,2*k);
norn = (sum(Sphn.^2,2)).^.5;
Sphn = Sphn ./ (norn * ones(1,2*k));


% generate map
f = zeros(N,1);
for j=1:Nn
    dj =  L2_distance(C',Sphn(j,:)');
    ik = find(dj == min(dj));
    f(j) = ik(1);
end

%% evaluate distortion of Zane's map
error = 0;
for i=1:Nn
    i 
    xi = Sphn(i,:);
    di = L2_distance(xi',Sphn');
    di = 2*asin(di/2);
    
    d1i = d1(f(i),f);
    
    ei = max(abs(di-d1i));
    error = max(error,ei);
end


disconj = 2*pi*(k-1)/(2*k-1) % this is the theoretical lower bound for distortion
error % this is the experimentally found distortion


