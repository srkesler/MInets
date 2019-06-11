% Shelli Kesler 3/7/19
% Compare covariance matrices using eigenvectors
% See Garcia, C 2012. BMC Evol Biol, 12, pp. 222. 

load('MInetsImpR1NimpR2.mat') %R1 = impared, R2 = not impaired
nperm = 5000;
rwire = 5;

rng(nperm);

%compute eigenvectors V from matrix A
[V1,~] = eig(R1);
[V2,~] = eig(R2);

%calculate eigenvalues D where AV=VD for variance explained by eigenvectors
%of original sample
D1X1 = R1*V1;
D2X2 = R2*V2;
%and by eigenvectors of the compared sample
D1X2 = R1*V2;
D2X1 = R2*V1;

%compute metrics
for i = 1:length(R1)
    S1temp(i,:) = (norm(D1X1(:,i)-D2X1(:,i))^2) + (norm(D1X2(:,i)-D2X2(:,i))^2);
    S2temp(i,:) = norm((D1X1(:,i)+D2X2(:,i)) - (D1X2(:,i)+D2X1(:,i)))^2;
    S3temp(i,:) = norm((D1X1(:,i)+D1X2(:,i)) - (D2X1(:,i)+D2X2(:,i)))^2;
end

%differention of covariance matrices as the ability of the eigenvectors of
%each sample to explain the variation in the other
S1 = 2*sum(S1temp);

%difference in orientation between eigenvectors in the same ordinal
%position in the two matrices

S2 = sum(S2temp);

%difference in proportion of total variance explained by eigenvectors in
%the same ordinal position in the two matrices

S3 = sum(S3temp);

%create random networks
for i = 1:nperm; Rand1(:,:,i) = randmio_und(R1,rwire); end
for i = 1:nperm; Rand2(:,:,i) = randmio_und(R2,rwire); end

%compute random eigenvectors
for i = 1:nperm
    [V1R(:,:,i),~] = eig(Rand1(:,:,i));
    [V2R(:,:,i),~] = eig(Rand2(:,:,i));
end

%compute random eigenvalues
for i = 1:nperm
    D1X1R(:,:,i) = Rand1(:,:,i)*V1R(:,:,i);
    D2X2R(:,:,i) = Rand2(:,:,i)*V2R(:,:,i);
    D1X2R(:,:,i) = Rand1(:,:,i)*V2R(:,:,i);
    D2X1R(:,:,i) = Rand2(:,:,i)*V1R(:,:,i);
end

% Calculate metrics for random networks
for j = 1:nperm
    for i = 1:length(R1)
    S1tempR(i,j) = (norm(D1X1R(:,i,j)-D2X1R(:,i,j))^2) + (norm(D1X2R(:,i,j)-D2X2R(:,i,j))^2);
    S2tempR(i,j) = norm((D1X1R(:,i,j)+D2X2R(:,i,j)) - (D1X2R(:,i,j)+D2X1R(:,i,j)))^2;
    S3tempR(i,j) = norm((D1X1R(:,i,j)+D1X2R(:,i,j)) - (D2X1R(:,i,j)+D2X2R(:,i,j)))^2;
    end
end

for j = 1:nperm
    S1R(:,j) = 2*sum(S1tempR(:,j));
    S2R(:,j) = sum(S2tempR(:,j));
    S3R(:,j) = sum(S3tempR(:,j));
end

% calculate 95% CIs and p values: determine the proportion (i.e. the mean) of times
% the random networks produce a value greater than actual value
lolim = .025*nperm;
hilim = nperm-lolim;

ci = {'LL', 'UL', 'S', 'p'};
for i = 1:3
    y = sort(eval(['S' num2str(i) 'R']));
        ci{i+1,1} = y(round(lolim));
        ci{i+1,2} = y(round(hilim));
        ci{i+1,3} = eval(['S' num2str(i)]);
        ci{i+1,4} = mean(eval(['S' num2str(i) 'R']) >= eval(['S' num2str(i)]));
end

ci
save eigenStatsMI ci;