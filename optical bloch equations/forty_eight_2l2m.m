clear;
warning('off','all');
format short;

import bloch.*
import TlF.*

time_of_flight=50000; %[ns]
pulse_time=300;%[ns]
switching_freq=1/(2*pulse_time);%[GHz]
no_switches=floor(time_of_flight/pulse_time);
t_step=time_of_flight/no_switches; %[ns]
% t_step=1*1000;

% Physical constants

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=1/99; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=0.32*a0; % [m]
B_0=6.686667; % [GHz]
T=7; %[K]
r_expval_m=5.22*r_expval;

%/////////////////////////////////////////////////////////////////////////%

% States
J_g=[0,1,2];
J_e=[1];

St=generateStates(J_g,J_e,0);
% n=size(St,1);

%/////////////////////////////////////////////////////////////////////////%
% Transition matrix

rabi_matrix=load('TransitionMatrix_EtoE_mixing_full.mat');
rabi_matrix=rabi_matrix.rabi_matrix;
n=size(rabi_matrix,1);


%/////////////////////////////////////////////////////////////////////////%

% 
% rabi_matrix=zeros(n,n,3);
% disp('Transition matrix')
% ind=0;
% for p=[-1,0,1]
%     ind=ind+1;
%     for f=1:n     %f - final states
%         for i=1:n %i - initial states
%             rabi_matrix(i,f,ind)=dipoleTransitionMatrixElement(St(f,:),St(i,:),p);
%         end
%     end
% end
% return

%Elimination of states that don't participate in the process
for i=104:-1:101
    rabi_matrix(i,:,:)=[];
    rabi_matrix(:,i,:)=[];
end

% rabi_matrix(101,:,:)=[];
% rabi_matrix(:,101,:)=[];

for i=100:-1:37
    rabi_matrix(i,:,:)=[];
    rabi_matrix(:,i,:)=[];
end

n=36+3+5;


%/////////////////////////////////////////////////////////////////////////%
% Initial distribution
disp('Initial distribution')

Ini=boltzmann_distribution(B_0*10^9,T,0:20);
IC=zeros(n);
ind=0;
for j=J_g    
    for i=1:4*(2*j+1)
        k=ind+i;
        IC(k,k)=Ini(j+1)/(4*(2*j+1));
    end
    ind=ind+4*(2*j+1);
end
IC=IC./sum(Ini(1:max(J_g)+1));


%/////////////////////////////////////////////////////////////////////////%

% Variables

syms w_0 w_1 w_2 w_e w_e2 real; %energies
syms w_m0 w_m1 w_L w_L2 real; %light frequencies
syms W_m0 W_m1 W_L W_L2 real; %Rabi frequencies
syms d_m0 d_m1 d_L d_L2 real; %detunings
syms D_0 D_1 D_10 D_11 D_2 D_21 D_22 real; %splittings
syms a_m0 p_m0 a_m1 p_m1 real; % microwave angles: a_m = polarization angle, p_m = beam direction angle with respect to quant axis z 
syms pol_L real; % laser polarization indicator variable: 1 - sigma+, 1 - sigma-

% Their values

Del_10=22.24*10^(-6)*2*pi;
Del_1=0.17595*10^(-3)*2*pi;
Del_11=14.54*10^(-6)*2*pi;
Del_21=44.52*10^(-6)*2*pi;
Del_2=0.27881*10^(-3)*2*pi;
Del_22=33.22*10^(-6)*2*pi;
Del_32=63.95*10^(-6)*2*pi;
Del_3=0.38458*10^(-3)*2*pi;
Del_33=54.38*10^(-6)*2*pi;
Del_0=13.3*10^(-6)*2*pi;
det_L=(Del_2/2+Del_21);
det_m0=0;
det_m1=0;
det_m3=0;
polangle_m0=pi/4;
dirangle_m0=pi/4;
polangle_m1=pi/4;
dirangle_m1=pi/4;
polangle_m3=pi/4;
dirangle_m3=pi/2;


% Light Properties

P_L=0.05; %[W]
diam_L=0.003; %[m]
P_L2=0.00; %[W]
diam_L2=0.003; %[m]
P_m0=0.000; %[W]
diam_m0=0.05; %[m]
P_m1=0.5; %[W]
diam_m1=0.05; %[m]

Rabi_m0=10^(-9)*q_e*r_expval_m*sqrt(8*P_m0/(pi*c*eps_0*diam_m0^2))/hbar; %[GHz]
Rabi_m1=10^(-9)*q_e*r_expval_m*sqrt(8*P_m1/(pi*c*eps_0*diam_m1^2))/hbar; %[GHz]
Rabi_L=10^(-9)*q_e*r_expval*sqrt(8*P_L/(pi*c*eps_0*diam_L^2))/hbar;
Rabi_L2=10^(-9)*q_e*r_expval*sqrt(8*P_L2/(pi*c*eps_0*diam_L2^2))/hbar;

disp('Rabi rate')
disp('R_L [\Gamma]')
disp(vpa(Rabi_L/Gamma,3))
disp('R_L2 [\Gamma]')
disp(vpa(Rabi_L2/Gamma,3))
disp('R_m0 [\Gamma]')
disp(vpa(Rabi_m0/Gamma,3))
disp('R_m0 [\Delta_1]')
disp(vpa(Rabi_m0/(2*pi*Del_1),3))
disp('R_m1 [\Gamma]')
disp(vpa(Rabi_m1/Gamma,3))

%/////////////////////////////////////////////////////////////////////////%
%%
% Branching ratios

disp('Branching ratios')

BR=load('BranchingRatios_EtoE_mixing_full.mat');
BR=BR.BR;


% BR=zeros(n);
% 
% transition_strengths=zeros(n);
% for i=1:n
%     for f=1:n
%         for p=1:3
%             transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(f,i,p)^2;
%         end
%     end
% end
% for i=1:n
%     sums=0;
%     for f=1:n
%         sums=sums+transition_strengths(i,f);
%     end
%     for f=1:n
%         BR(i,f)=transition_strengths(i,f)/sums;
%     end
% end

for i=104:-1:101
    BR(i,:,:)=[];
    BR(:,i,:)=[];
end

% BR(101,:,:)=[];
% BR(:,101,:)=[];

for i=100:-1:37
    BR(i,:,:)=[];
    BR(:,i,:)=[];
end


for i=1:n-8
    BR(i,:)=0;
end

%/////////////////////////////////////////////////////////////////////////%

% Dissipator

disp('Dissipator')

L=Dissipator(n);

syms G real;
assume(G,'positive') 

%Current implementation includes both J'=1, F'1=3/2, F'=1,2
DR=zeros(1,n-8);
DR=[DR,G,G,G,G,G,G,G,G];


L.fromBranching(BR,DR);

%/////////////////////////////////////////////////////////////////////////%
%%
syms t pol_mod pol_mod_m real;

%Laser 1
WLp=W_L/sqrt(2)*sin(pi/2*(1+sign(sin(2*pi*switching_freq*t)))/2)*pol_mod;
WLm=W_L/sqrt(2)*sin(pi/2*(1+sign(sin(2*pi*switching_freq*t)))/2)*pol_mod;
WLz=W_L*cos(pi/2*(1+sign(sin(2*pi*switching_freq*t)))/2)*pol_mod+(1-pol_mod)*W_L;


WL_pol=[WLm,WLz,-WLp];

%Laser 2
WLp2=0;
WLm2=0;
WLz2=W_L2;

WL_pol2=[WLm2,WLz2,-WLp2];

%Microwaves
Wm1p=W_m1/sqrt(2)*(sin(p_m1)+1i*cos(a_m1)*cos(p_m1));
Wm1m=W_m1/sqrt(2)*(sin(p_m1)-1i*cos(a_m1)*cos(p_m1));
Wm1z=-W_m1*cos(p_m1)*sin(a_m1);

Wm1_pol=[Wm1m,Wm1z,-Wm1p];

Wm1_pol=subs(Wm1_pol,p_m1,pi/2*(1-pol_mod_m)+pol_mod_m*pi/2*(1+sign(sin(2*pi*switching_freq*t+pi/2)))/2);

%/////////////////////////////////////////////////////////////////////////%
%%
% Hamiltonian

import bloch.*
import TlF.*

disp('Hamiltonian')

H=Hamiltonian(n);

H.addEnergies([w_0,w_0-D_0,w_0-D_0,w_0-D_0,...
    w_1,w_1+D_10,w_1+D_10,w_1+D_10,...
    w_1+D_1+D_10,w_1+D_1+D_10,w_1+D_1+D_10,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11,w_1+D_1+D_10+D_11...
    w_2,w_2,w_2,w_2+D_21,w_2+D_21,w_2+D_21,w_2+D_21,w_2+D_21,...
    w_2+D_2+D_21,w_2+D_2+D_21,w_2+D_2+D_21,w_2+D_2+D_21,w_2+D_2+D_21,...
    w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,w_2+D_2+D_21+D_22,...
    w_e,w_e,w_e,w_e2,w_e2,w_e2,w_e2,w_e2]);


for i=1:n
    for f=i:n
        for j=1:3
            %Laser (J=2 to Je=1 F1=3/2 F=1)
            if rabi_matrix(i,f,j)~=0 && f<=n-5 && f>n-8  && i<=n-8 && i>4
                Wr=(-1)^(j)*WL_pol(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,conj(Wr),w_L);
%             Laser (J=0 to Je=1 F1=3/2 F=2)
            elseif rabi_matrix(i,f,j)~=0 && f>n-5 && i<=4
                Wr=(-1)^(j)*WL_pol2(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,conj(Wr),w_L2);
            %Microwaves (J=1 to J=2)    
            elseif rabi_matrix(i,f,j)~=0 && f<=n-8 && f>16 && i<=n-8-20 && i>4
                Wr=(-1)^(j)*Wm1_pol(4-j)*rabi_matrix(i,f,j);
                H.addCoupling(i,f,conj(Wr),w_m1);
            %Microwaves (J=0 to J=1)
%             elseif rabi_matrix(i,f,j)~=0 && f<=16 && i<=4
%                 Wr=(-1)^(j)*Wm0_pol(j)*rabi_matrix(i,f,j);
%                 H.addCoupling(i,f,Wr,w_m0);
            end
        end
    end
end



H.defineEnergyDetuning(w_2,w_e,d_L,w_L);
H.defineEnergyDetuning(w_1,w_2,d_m1,w_m1);
% H.defineEnergyDetuning(w_0,w_e,d_L2,w_L2); 
% H.defineStateDetuning(1,5,d_m0);


H.defineZero(w_2);
H.unitaryTransformation();


H.defineZero(w_0);
H.unitaryTransformation();
H.defineStateDetuning(2,42,d_L2);

% H.addSidebands(w_L,9,8.5,1*Gamma);

disp(H.transformed)


%/////////////////////////////////////////////////////////////////////////%
%%
import bloch.*
% Bloch Equations

disp('Optical Bloch Equations')


Eq=BlochEqns(H,L);

Eq.eqnsRHS=subs(Eq.eqnsRHS,[D_1,D_2,D_0,D_10,D_11,D_21,D_22,G,a_m1,w_e2,w_0],[Del_1,Del_2,Del_0,Del_10,Del_11,Del_21,Del_22,Gamma,dirangle_m1,0,0]);

Eq.necessaryVariables();

return
%/////////////////////////////////////////////////////////////////////////%
%%

% for v=0
v=0;
R_m0=0;
R_m1=Rabi_m1;
R_L2=0;
R_L=Rabi_L;

fprintf('v = %.2f \n',v);
de_L=det_L+v;

disp('Solving')
Eq.evolve(0,time_of_flight,IC,[R_L,R_L2,R_m1,de_L,0,0,1,1]);
disp('Solved')


disp('J0')
disp(real(Eq.evolution(1,1,end)))
disp(real(Eq.evolution(2,2,end)))
disp(real(Eq.evolution(3,3,end)))
disp(real(Eq.evolution(4,4,end)))
disp('Js')
disp(real(Eq.evolution(1,1,end)+Eq.evolution(2,2,end)+Eq.evolution(3,3,end)+Eq.evolution(4,4,end)))
disp(real(Eq.evolution(5,5,end)+Eq.evolution(6,6,end)+Eq.evolution(7,7,end)+Eq.evolution(8,8,end)+Eq.evolution(9,9,end)+Eq.evolution(10,10,end)+Eq.evolution(11,11,end)+Eq.evolution(12,12,end)+Eq.evolution(13,13,end)+Eq.evolution(14,14,end)+Eq.evolution(15,15,end)+Eq.evolution(16,16,end)))
disp(real(Eq.evolution(17,17,end)+Eq.evolution(18,18,end)+Eq.evolution(19,19,end)+Eq.evolution(20,20,end)+Eq.evolution(21,21,end)+Eq.evolution(22,22,end)+Eq.evolution(23,23,end)+Eq.evolution(24,24,end)+Eq.evolution(25,25,end)+Eq.evolution(26,26,end)+Eq.evolution(27,27,end)+Eq.evolution(28,28,end)+Eq.evolution(29,29,end)+Eq.evolution(30,30,end)+Eq.evolution(31,31,end)+Eq.evolution(32,32,end)+Eq.evolution(33,33,end)+Eq.evolution(34,34,end)+Eq.evolution(35,35,end)+Eq.evolution(36,36,end)))
disp(real(Eq.evolution(44,44,end)+Eq.evolution(43,43,end)+Eq.evolution(42,42,end)+Eq.evolution(41,41,end)+Eq.evolution(40,40,end)+Eq.evolution(38,38,end)+Eq.evolution(39,39,end)+Eq.evolution(37,37,end)))


% end

%/////////////////////////////////////////////////////////////////////////%
%Plots
%%
figure
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(1,1,:)+Eq.evolution(2,2,:)+Eq.evolution(3,3,:)+Eq.evolution(4,4,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(5,5,:)+Eq.evolution(6,6,:)+Eq.evolution(7,7,:)+Eq.evolution(8,8,:)+Eq.evolution(9,9,:)+Eq.evolution(10,10,:)+Eq.evolution(11,11,:)+Eq.evolution(12,12,:)+Eq.evolution(13,13,:)+Eq.evolution(14,14,:)+Eq.evolution(15,15,:)+Eq.evolution(16,16,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(17,17,:)+Eq.evolution(18,18,:)+Eq.evolution(19,19,:)+Eq.evolution(20,20,:)+Eq.evolution(21,21,:)+Eq.evolution(22,22,:)+Eq.evolution(23,23,:)+Eq.evolution(24,24,:)+Eq.evolution(25,25,:)+Eq.evolution(26,26,:)+Eq.evolution(27,27,:)+Eq.evolution(28,28,:)+Eq.evolution(29,29,:)+Eq.evolution(30,30,:)+Eq.evolution(31,31,:)+Eq.evolution(32,32,:)+Eq.evolution(33,33,:)+Eq.evolution(34,34,:)+Eq.evolution(35,35,:)+Eq.evolution(36,36,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(37,37,:)+Eq.evolution(38,38,:)+Eq.evolution(39,39,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(40,40,:)+Eq.evolution(41,41,:)+Eq.evolution(42,42,:)+Eq.evolution(43,43,:)+Eq.evolution(44,44,:)))
legend('J=0','J=1','J=2','Je=1, F=1','Je=1, F=2')
xlabel('Time [us]')
drawnow
%%
figure
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(1,1,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(2,2,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(3,3,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(4,4,:)))
legend('F=0, M_F=0','F=1 M_F=-1','F=1 M_F=0','F=1 M_F=1')
xlabel('Time [us]')
drawnow
% %%
% Time=Eq.evTime;
% Populations=[];
% for i=1:44
%     Populations=[Populations;real(squeeze(Eq.evolution(i,i,:))')];
% end
% Data=[Time;Populations];
% save('RC_poldep_P2F1_and_R0F2.mat','Data')
% return




