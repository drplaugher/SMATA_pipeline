% Code for some examples in thesis on 
% Assembling an Integrated Pipeline to Construct Patient-Specific Cancer Models
% By Daniel Plaugher 9/10/22

% home 
clear
clc
close all
addpath('/Users/drpla/Dropbox/BooleanNetworks');
addpath('/Users/drpla/Dropbox/BooleanNetworks/BNPBN/');
addpath('/Users/drpla/Dropbox/BooleanNetworks/BNPBN/BNPBN');
addpath('/Users/drpla/Dropbox/BooleanNetworks/BNPBN/BNPBN/BNPBN');

% office
% addpath('/Users/drpl222/Dropbox/BooleanNetworks/');
% addpath('/Users/drpl222/Dropbox/BooleanNetworks/BNPBN/');
% addpath('/Users/drpl222/Dropbox/BooleanNetworks/BNPBN/BNPBN/BNPBN/');
% addpath('/Users/drpl222/Dropbox/BooleanNetworks/BNPBN/BNPBN/BNPBN/');



%%%%     The following sections are for defining example models
%%  Simple 3-cycle
n = 3; 
p = 2;

%%% for SDDS_Build
syms x1 x2 x3

f= [x3
    x1;
    x2]; 
[varF,nv,F]=SDDS_Build(syms,f,p);

%% CK example
n = 6; 
p = 2;

%%% for SDDS_Build
syms x1 x2 x3 x4 x5 x6

f= [x4+1;
    x1;
    x1;
    x5+1;
    x3*x6;
    x2]; 
[varF,nv,F]=SDDS_Build(syms,f,p);

%% Modularity Example
n = 6; 
p = 2;

% for SDDS_Build
syms x1 x2 x3 x4 x5 x6 

f= [x3;
    x1;
    x2;
    x1*x6+x1+x6;
    x4;
    x5]; 
[varF,nv,F]=SDDS_Build(syms,f,p);


% Set up for di-graph
[row,col]=find(varF~=-1);
t=NaN(1,length(row)); % tail
h=NaN(1,length(row)); % head
for i=1:length(row)
    t(1,i)= varF(row(i),col(i)); %builds vector of tails
    h(1,i) = col(i); % builds vector of heads
end

G = digraph(t,h,[],n);
plot(G,'Layout','layered')

% Finding modules
[bin,binsize] = conncomp(G);
mods = find(binsize~=1); % modules found
mod_nodes=-1*ones(max(binsize),length(mods)); % nodes in each module in columns
for i=1:length(mods)
    mod_nodes(1:binsize(mods(i)),i)=find(bin==mods(i))';   
end
modTable= [mods; NaN(1,length(mods)); mod_nodes]; % table with modules and their nodes

% to find module rankings
pl = plot(G);
pl.MarkerSize = 7;
pl.NodeCData = bin;
colormap(hsv(length(binsize)))

C = condensation(G);
pl2 = plot(C);
pl2.MarkerSize = 7;
pl2.NodeCData = 1:length(binsize);
colormap(hsv(length(binsize)))

%% CA Example
n = 5; 
p = 2;

%%% for SDDS_Build
syms x1 x2 x3 x4 x5 

f= [1+x3+x5+x3*x5;
    1+x1+x1*x4;
    1+x2;
    x3;
    1+x4]; 
[varF,nv,F]=SDDS_Build(syms,f,p);


%% 

%%%%% The sections below are for analysing the defined models above

%% Set up for Markov Chains
c = 0.9*ones(2,n); % normal propensities for SDDS
 
[~,Avec] = bnAsparse(F,varF,nv);
% shows the atractors and their respective basins
[ab,dd] = bnAttractor(Avec); % dd~ steps from attr.
attrs =  unique(ab(ab<0)); % finds attractors

disp(attrs)
a1= find(ab==-1); % finds attractor 
a2= find(ab==-2); % finds attractor 
a3= find(ab==-3); % finds attractor 
a4= find(ab==-4); % finds attractor 
 
% to find binary representation of attractors
att=a2;
for i=1:length(att)
    x = dec2multistate(att(i)-1,p,n) % binary represention of attractor
end


Trans=multistateA(F,varF,nv,c,p); % transistion matrix--probability of moving from one node to another 

%%%%% Google Matrix %%%%%
K=(1/2^n)*ones(2^n,2^n); % K matrix for noise
g=0.9;
G=g*Trans+(1-g)*K; % google matrix, adds noise to make regular for Perron-Frobenius Thm
G_mc=dtmc(G); % markov chain for google matrix
G_pwr=G^10000; % high powers make columns converge to stationary distribution
G_dist=asymptotics(G_mc); % stationary distribution of google matrix--time spent at node
[B, I]=sort(G_dist,'descend'); % sorts in descending order and oututs associated index
Ranks=[I' B']; % table with index and dist.

%% Mutations and Control

% for nodes
F1 = TruthTable_del_n_temp(F,nv,varF,p, 1,0); % regulates 
F2 = TruthTable_del_n_temp(F1,nv,varF,p, 6,0); % regulates 

% for edges
%  F1 = TruthTable_del_a_temp(F,nv,varF,p,2,3,0); %edge control
%  F2 = TruthTable_del_a_temp(F1,nv,varF,p,3,1,0); %edge control

%% Simulations
nins = 1000; % number of initializations
nsteps=100; % number of steps for SDDS
g=0.01; % noise

[Y,My]=SDDS_simNoise(g,F,varF,nv,p,c,n, nsteps,nins); 
%[Y,My]=SDDS_sim(F1,varF,nv,p,c,n, nsteps,nins); % simulation w/o noise
Ylast=Y(:,end); % long-term trajectories
 
X = 0:1:nsteps; % time steps

figure('Name', 'Simulation')
plot(X,Y)
%plot(X,Y(1,:),'b',X,Y(2,:),'m',X,Y(3,:),'k',Y(4,:),Y(5,:),'LineWidth',1.5,'MarkerSize',10)%,,X,Y(4,:),'y'
legend('x1','x2', 'x3','x4','x5')
xlabel('Time Steps')
ylabel('Average Frequencies')
