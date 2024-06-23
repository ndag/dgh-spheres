%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Facundo Memoli                                     %
%                                                                    %
% See https://github.com/ndag/dgh-spheres                            %
% Version of April/2024                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTH = 200;
nT = 500;

TH = linspace(-pi/3,pi/3,nTH);
Z = acot(3*(3-4*sin(TH).^2));


%% (**)  We'll test inequality \cos(t-2 \pi/3) \geq max_{th,th'} F_t(th,th') 
T = linspace(2*pi/3,pi,nT);

M = zeros(1,nT);
for k=1:nT
    t = T(k);
    Ft = zeros(nTH,nTH);
    for i=1:nTH
        thi = TH(i);
        zi = Z(i);
        for j=1:nTH
            thj = TH(j);
            zj = Z(j);
            Ft(i,j) = cos(thi-thj+t) * cos(zi) * cos(zj) + cos(3*(thi-thj+t)) * sin(zi) * sin(zj);
        end
    end
    M(k) = max(Ft(:));
end


plot(T,cos(T-2*pi/3),T,M)
legend('cos(t-2*pi/3)', 'max Ft')
title('Inequality (**)')

nv = sum(M > cos(T-2*pi/3));
disp(['Inequality (**) was violated ' num2str(nv) ' times out of ' num2str(nT) ' total values in the range of t'])


%% (***)  We'll test the inequality \cos(t+2 \pi/3) \geq max_{th,th'} F_t(th,th') 
figure

T = linspace(-pi,-2*pi/3,nT);

M = zeros(1,nT);
for k=1:nT
    t = T(k);
    Ft = zeros(nTH,nTH);
    for i=1:nTH
        thi = TH(i);
        zi = Z(i);
        for j=1:nTH
            thj = TH(j);
            zj = Z(j);
            Ft(i,j) = cos(thi-thj+t) * cos(zi) * cos(zj) + cos(3*(thi-thj+t)) * sin(zi) * sin(zj);
        end
    end
    M(k) = max(Ft(:));
end


plot(T,cos(T+2*pi/3),T,M)
legend('cos(t+2*pi/3)', 'max Ft')
title('Inequality (***)')

nv = sum(M > cos(T+2*pi/3));

disp(['Inequality (***) was violated ' num2str(nv) ' times out of ' num2str(nT) ' total values in the range of t'])


