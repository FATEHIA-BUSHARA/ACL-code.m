function [Vf,Sigma,Beta, Results] = ACL_Strain(Diameter,Volume, Theta,lamda3,Force,edges)
% Analysis of ACL strained in microCT – 08/02/2021
% X, Z, (Lamda) for each angular interval
% F/A , Beata  for each deformation (and entropy, and energy)
% Ef for each sample
% X volume fraction each angular interval
% Z angular frequency of that interval(orientation)
% F/A represents the Sigma (lamda33)and
% Sigma is given for each diformation
% F is given for each diformation 
% A is calculated as sum(pi*(Thi /2)^2), (Thi is fascicles Diameter)
% The best equation fitting of  Sigma33 (Lamda3) on all the lamda3 gives Ef
% edges is Interval should be a vector for ex... edges = 0:10:90;

% % % % % % % % % % Radian %%%%%%%%%%%%%%%%%%%%
R = deg2rad(Theta);
% edges = 0:10:90;  % 0:0.1745:1.6; = deg2rad (BinEdges) ,BinWidth = 0.1745
BinEdges = deg2rad(edges);
Z= histcounts(R,BinEdges,'Normalization', 'probability');

% figure
% h = histogram(R,10 ,'Normalization','probability' );
% h.BinEdges = deg2rad(edges);
% Z = h.Values;  % %% Relative angular frequency in radian 
% xlabel('Orientation')
% ylabel('Relative Frequency')
% title('TK027-T2 Orientation Distribution ')
% title('TK024-T2 Orientation Distribution ')
% saveas(figure(1),'C:\Users\PhDs\OneDrive\ACL measurment\Figures\dif Orientation/P Orientation T2.jpg');


% % % % % % % % % Orientation on degree %%%%%%%%%%
% figure
% h1 = histogram(Theta,10 ,'Normalization','probability' );
% h1.BinEdges = 0:10:90;
% xlabel('Orientation')
% ylabel('Relative Frequency')
% Z1 = h1.Values; % Relative angular frequency in degree 

%% %%  X volume fraction each angular interval 
[X,Vf] = Volfraction(Volume,Theta,edges);
xy = Vf.* R; % Volume fraction using angle in radian
% xy = Vf.* B(:, 3);% Volume fraction using angle in degree
% c = [xy B(:, 3)];
En = entropy(xy);
Eg = energy(xy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy1 = Vf.* Theta; 
% figure
% h11 = histogram(xy1,10 ,'Normalization','probability' );
% h11.BinEdges = 0:10:90;
% % edges = 0:10:90;  % 0:0.1745:1.6; = deg2rad (BinEdges) ,BinWidth = 0.1745
% % h11.BinEdges = deg2rad(edges);
% xlabel('Weighted Orientation')
% ylabel('Relative Frequency')


%% Lamda calculation for each interval N
% Theta is the central frequency (phi)
% Y represent the third part of the equation model which is:
% Y = (1-1/Lamda).* cos(Theta).^2 
% cedges = 5:10:85;            % central angle in degree
cedges = edges+5;
cedges = cedges(1:(length(cedges)-1));
Lamda=zeros(1,numel(cedges)); % Template for the Lamda at each interval
Y=zeros(1,numel(cedges));     % Template for the Y at each interval

% lamda3 = 1.11111; % % at T0

for k = 1: numel(cedges)
Theta = deg2rad(cedges (k));     % Central angle in radian   
Lamda(k) = sqrt(((lamda3)^2 *(cos(Theta)).^2) + ((1/lamda3)*(sin(Theta)).^2));
Y (k) = (1-(1./Lamda(k))).* cos(Theta).^2 ;
end

%% Beta represents the left side of the equation
%%% Sigma represents the right side of the equation F/A 
% F is the resistance to elongation given from the experiment
% A = sum(pi*(Diameter/2).^2);

Beta = 2*pi*lamda3^2 *(sum(X .* Y .* Z));

A = sum(pi*(Diameter./2).^2);
% F = 8.4; % at T0
Sigma = Force/A;

%% Diameter measurements (mean , variance)and histogram distribution
% Diameter histogram, mean , variance, standard deviation
dmean = mean(Diameter);
dvar = var(Diameter);
% dstd = std(Diameter);
mes1 = [dmean, dvar];

% A is Area ,En entropy, Eg energy, mes1 diameter mean and variance

Results = [ A, En,Eg, mes1];
end

%% Compute energy

function eg=energy(p)
[n]=length(p);
eg=0;
for i=1:n
        eg=eg+p(i)^2;   
end
end

%% compute the volume fraction 

function [X,Vf] = Volfraction(a,Theta,n)

%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = Volume vector
% Theta = Orientation
% n= edges (vector)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = Volume fraction for each interval
% cm = total for each interval (vetor of 9 scalaar)
% Vf : volume fraction for each object respect to spesific interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % a1= [ 1 1 2 2 2  5 5  7 7 8 9]
% % % % Theta = a1'
% % % % a = [rand(length(Theta),1)]
% % % % n= 0:3:10
% % % % n1 = histcounts(Theta,n)

% a =(B(:, 2)); Volume
% Theta = B(:, 3); Orientation
% n =  0:10:90;
%  n1 = histcounts(Theta,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1 = histcounts(Theta,n); % Number of frequency
cm=zeros(1,length(n1));   % Template for the Total volume for each interval
Vf=zeros(length(a),1);    % Template for the volume fraction

%% Total volume for each interval

for h = 1:length(n1)
    s = 0;
    for i = 1:length(a)
        if (Theta(i) >= n(h) )&& (Theta(i)<= n(h+1))
        s= s+a(i);
        end
    cm(h)=s;
    end
end
X= cm/sum(a);

%% Vf volume fraction

for h = 1:length(n1)
    for i = 1:length(a)
        if (Theta(i) >= n(h) )&& (Theta(i)<= n(h+1))
        v(i)=  a(i)/cm(h);
        end
    end
    Vf=v';
end

end
