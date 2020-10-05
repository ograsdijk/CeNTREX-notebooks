% clear all;
import bloch.*
import TlF.*

warning('off','all');

time_of_flight=100; %[us]
switching_freq=15/100;%[MHz]
no_switches=floor(time_of_flight*switching_freq);
t_step=time_of_flight*1000/no_switches; %[ns]


%Physical constants 

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=1/99; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=0.32*a0; % [m]
mu_r_expval=5.22*r_expval;

% States
J_g=[0,1,3];
J_e=[1];

St=generateStates(J_g,J_e,1);
n=size(St,1);


% Transition matrix
% rabi_matrix=load('TransitionMatrix_EtoF_nomixing.mat');
% rabi_matrix=rabi_matrix.rabi_matrix;
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


% rabi_matrix(end-11,:,:)=[];
% rabi_matrix(:,end-11,:)=[];
% rabi_matrix(end-7:end,:,:)=[];
% rabi_matrix(:,end-7:end,:)=[];

% rabi_matrix(end-11:end-8,:,:)=[];
% rabi_matrix(:,end-11:end-8,:)=[];
% rabi_matrix(end-4:end,:,:)=[];
% rabi_matrix(:,end-4:end,:)=[];
%%
rabi_matrix=load('TransitionMatrix_EtoE_mixing.mat');
rabi_matrix=rabi_matrix.rabi_matrix;

X1=rabi_matrix(:,:,1)+rabi_matrix(:,:,2)+rabi_matrix(:,:,3);
Z1=X1;
for i=2:76
    for j=1:i-1
        Z1(i,j)=0;
    end
end
Z1=Z1+Z1';
Y1=X1-X1';
Y1(Y1<1e-10)=0;

%%
TransMatrix=Z;
filename="TransitionDipoles_EtoF.mat";
% % 
save(filename,'TransMatrix');
%%
% disp(abs(rabi_matrix(:,:,1))-abs(rabi_matrix(:,:,1)'))
% disp(rabi_matrix(:,:,2)+rabi_matrix(:,:,2)')
% disp(rabi_matrix(:,:,3)+rabi_matrix(:,:,3)')

% rabi_matrix(end-11,:,:)=[];
% rabi_matrix(:,end-11,:)=[];
% rabi_matrix(end-7:end,:,:)=[];
% rabi_matrix(:,end-7:end,:)=[];
% rabi_matrix(1:4,:,:)=[];
% rabi_matrix(:,1:4,:)=[];
rabi_matrix(17:36,:,:)=[];
rabi_matrix(:,17:36,:)=[];
rabi_matrix(end-11,:,:)=[];
rabi_matrix(:,end-11,:)=[];
rabi_matrix(end-7:end,:,:)=[];
rabi_matrix(:,end-7:end,:)=[];
rabi_matrix(end-25:end-3,:,:)=[]; %elimination of most J=3 levels
rabi_matrix(:,end-25:end-3,:)=[];

n=n-9-23;

% n=n-9;
%%

% Branching ratios
disp('Branching ratios')
% BR=load('BranchingRatios_EtoE_nomixing.mat');
% BR=BR.BR;


BR=zeros(n);

transition_strengths=zeros(n);
for i=1:n
    for f=1:n
        for p=1:3
            transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(f,i,p)^2;
        end
    end
end
for i=1:n
    sums=0;
    for f=1:n
        sums=sums+transition_strengths(i,f);
    end
    for f=1:n
        BR(i,f)=transition_strengths(i,f)/sums;
    end
end

for i=1:n-3
    BR(i,:)=0;
end

% format rat
% disp(BR)
% format short

% Dissipator


disp('Dissipator')

L=Dissipator(n);

syms G real;
assume(G,'positive') 



DR=zeros(1,n-3);
DR=[DR,G,G,G];


L.fromBranching(BR,DR);



%%


% Variables

import bloch.*

syms w_0 w_1 w_e w_e2 w_3 real; %energies
syms w_m w_L w_L2 real; %light frequencies
syms W_m W_L W_L2 real; %Rabi frequencies
syms d_m d_L d_L2 real; %detunings
syms D_0 D_11 D_1 D_12 real; %splittings
syms a_m a_L real; %polarization angles
syms p_m real; %microwave angle with respect to y (propagation direction of the laser)
syms pol_L pol_L2 pol_m real;

assume(W_L,'positive');
% assume(W_L2,'positive');
assume(W_m,'positive');
%%
%Hamiltonian
import bloch.*
disp('Hamiltonian')
H=Hamiltonian(n);
% H.addEnergies([w_1,w_1+D_11,w_1+D_11,w_1+D_11,...
%     w_1+D_11+D_1,w_1+D_11+D_1,w_1+D_11+D_1,...
%     w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12...
%     w_e,w_e,w_e,w_e2,w_e2,w_e2,w_e2,w_e2]);
H.addEnergies([w_0+D_0,w_0,w_0,w_0,...
    w_1,w_1+D_11,w_1+D_11,w_1+D_11,...
    w_1+D_11+D_1,w_1+D_11+D_1,w_1+D_11+D_1,...
    w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12,w_1+D_11+D_1+D_12...
    w_3,w_3,w_3,w_3,w_3,...
    w_e,w_e,w_e]);

%Laser
WLp=W_L*pol_L/sqrt(2);
WLm=W_L*pol_L/sqrt(2);
WLz=W_L*(1-pol_L);


WL_pol=[WLm,WLz,WLp];



Wmp=W_m*exp(1i*pi/4)*pol_m;
Wmm=W_m*exp(-1i*pi/4)*pol_m;
Wmz=W_m*(1-pol_m);

Wm0_pol=[Wmm,Wmz,Wmp];

% Wmp=W_m*(sin(p_m)+1i*cos(a_m)*cos(p_m))/sqrt(2);
% Wmm=W_m*(sin(p_m)-1i*cos(a_m)*cos(p_m))/sqrt(2);
% Wmz=W_m*cos(p_m)*sin(a_m);

% Wm0_pol=[Wmm,Wmz,Wmp];


%Couplings
for i=1:n
    for f=i:n
        for j=1:3
            if rabi_matrix(i,f,j)~=0 && f>n-3 && i<=n-8
                Wr=(-1)^(j)*WL_pol(j)*abs(rabi_matrix(i,f,j));
                H.addCoupling(i,f,Wr,w_L);
%             elseif rabi_matrix(i,f,j)~=0 && f<=n-5 && f>n-8 && i<=n-8
%                 Wr=(-1)^(j)*WL_pol2(j)*rabi_matrix(i,f,j);
%                 H.addCoupling(i,f,Wr,w_L2);
            elseif rabi_matrix(i,f,j)~=0 && f<=n-8 && i<=4
                Wr=(-1)^(j)*Wm0_pol(j)*abs(rabi_matrix(i,f,j));
                H.addCoupling(i,f,Wr,w_m);
            end
        end
    end
end


H.defineEnergyDetuning(w_1,w_e,d_L,w_L);
H.defineEnergyDetuning(w_0,w_1,d_m,w_m);
% H.defineEnergyDetuning(w_1,w_e2,d_L2,w_L2);


H.defineZero(w_0);

H.unitaryTransformation();

disp(H.transformed)


%%
Del_0=2*pi*13.3*10^(-6);
Del_1=2*pi*176*10^(-6);
Del_11=2*pi*22.24*10^(-6);
Del_12=2*pi*14.54*10^(-6);
det_L=0;
det_m=0;

% Light Properties
%%


P_L=0.01; %[W]
diam_L=0.0014; %[m]
P_m=0.4; %[W]
diam_m=0.03; %[m]


Rabi_m=10^(-9)*q_e*mu_r_expval*sqrt(8*P_m/(pi*c*eps_0*diam_m^2))/hbar; %[GHz]
Rabi_L=10^(-9)*q_e*r_expval*sqrt(8*P_L/(pi*c*eps_0*diam_L^2))/hbar;

disp(vpa(Rabi_L/Gamma,5))
disp(vpa(Rabi_m/Del_1,5))
disp(vpa(Rabi_m/Gamma,5))


%%
import bloch.*
Eq=BlochEqns(H,L);

%%
IC=zeros(n);

for i=5:16
IC(i,i)=1/12;
end

Eq.eqnsRHS=subs(Eq.eqnsRHS,[G,D_0,D_11,D_1,D_12,w_3],[Gamma,Del_0,Del_11,Del_1,Del_12,0]);
%%
Eq.necessaryVariables();
return
%%
import bloch.*
no_switches=100;
det_m=0*Gamma;
% % Rabi_L=2.5*Gamma;
% Rabi_m=19.4*Gamma;
% angle_m=pi/4;
% pol_angle=pi/6;

Time_intervals=zeros(1,no_switches+1);
% Fields_time=zeros(1,2*no_switches);
% Laser_fields_p=zeros(1,2*no_switches);
% Laser_fields_m=zeros(1,2*no_switches);
% Micro_fields1=zeros(1,2*no_switches);
disp('Solving')
for sw=1:no_switches
    
    if mod(sw,4)==1
        polarization_L=1;
        polarization_m=0;
        R_L=Rabi_L;
        R_m=Rabi_m;
        t_step=10/Gamma;
    elseif mod(sw,4)==2
        polarization_L=0;
        polarization_m=0;
        R_L=Rabi_L;
        R_m=Rabi_m;
        t_step=10/Gamma;
    elseif mod(sw,4)==3
        polarization_L=1;
        polarization_m=1;
        R_L=Rabi_L;
        R_m=Rabi_m;
        t_step=10/Gamma;
    else
        polarization_L=0;
        polarization_m=1;
        R_L=Rabi_L;
        R_m=Rabi_m;
        t_step=10/Gamma;
 
    end

    Time_intervals(sw+1)=Time_intervals(sw)+t_step;
%     Fields_time(2*sw-1)=Time_intervals(sw);
%     Fields_time(2*sw)=Time_intervals(sw+1)-0.001;
%     Laser_fields_p(2*sw-1)=R_L*(polarization_L+1)/2;
%     Laser_fields_p(2*sw)=R_L*(polarization_L+1)/2;
%     Laser_fields_m(2*sw-1)=R_L*(1-polarization_L)/2;
%     Laser_fields_m(2*sw)=R_L*(1-polarization_L)/2;
%     Micro_fields1(2*sw-1)=R_m;
%     Micro_fields1(2*sw)=R_m;

if sw==1
    Eq.evolve(0,t_step,IC,[R_L,R_m,det_L,det_m,polarization_L,polarization_m]);
else
    Eq.intTime=[Time_intervals(sw-1),Time_intervals(sw)];
    Eq.lastSol=prevSol;
    Eq.extendEvolution(Time_intervals(sw+1),[R_L,R_m,det_L,det_m,polarization_L,polarization_m]);
    if sw==floor(no_switches/7)||sw==floor(no_switches/4)||sw==floor(no_switches/3)||sw==floor(no_switches/2)||sw==floor(2*no_switches/3)||sw==floor(4*no_switches/5)||sw==floor(6*no_switches/7)
        Eq.evolution(:,:,2:2:end)=[];
    end
end

prevSol=Eq.lastSol;

fprintf('Done %.2f%% \n',sw/no_switches*100);

end

%%
u1=zeros(1,19);
u1(1)=1;
u1=u1./norm(u1);
u2=zeros(1,19);
u2(2)=1;
u2=u2./norm(u2);
u3=zeros(1,19);
u3(3)=1;
u3=u3./norm(u3);
u4=zeros(1,19);
u4(4)=1;
u4=u4./norm(u4);

u5=zeros(1,19);
u5(5)=-2/3;
u5(7)=2/3;
u5(10)=-sqrt(2)/6;
u5(14)=sqrt(2)/6;
u5=u5./norm(u5);
u6=zeros(1,19);
u6(5)=-2/3;
u6(7)=-2/3;
u6(10)=sqrt(2)/6;
u6(14)=sqrt(2)/6;
u6=u6./norm(u6);
u7=zeros(1,19);
u7(6)=sqrt(6)/3;
u7(9)=-sqrt(3)/6;
u7(13)=1/2;
u7=u7./norm(u7);
u8=zeros(1,19);
u8(8)=-sqrt(6)/3;
u8(11)=sqrt(3)/6;
u8(15)=1/2;
u8=u8./norm(u8);
u9=zeros(1,19);
u9(12)=1;
u10=zeros(1,19);
u10(16)=1;


u11=zeros(1,19);
u11(6)=-1/12;
u11(8)=0.2357;
u11(9)=-0.2754;
u11(11)=11/12;
u11(13)=-0.0229;
u11(15)=-0.1443;
u11=u11./norm(u11);
u12=zeros(1,19);
u12(5)=-1/9;
u12(6)=-0.22;
u12(7)=-0.044;
u12(8)=-0.3111;
u12(9)=-0.7271;
u12(10)=-0.1257;
u12(11)=-0.2200;
u12(13)=-0.0606;
u12(14)=-0.3143;
u12(15)=-0.3810;
u12=u12./norm(u12);
u13=zeros(1,19);
u13(6)=-1/2;
u13(9)=0.0795;
u13(13)=0.8624;
u13=u13./norm(u13);
u14=zeros(1,19);
u14(5)=0.3143;
u14(6)=-0.0778;
u14(7)=-0.0157;
u14(8)=-0.1100;
u14(9)=-0.2571;
u14(10)=-0.0444;
u14(11)=-0.0778;
u14(13)=-0.0214;
u14(14)=8/9;
u14(15)=-0.1347;
u14=u14./norm(u14);
u15=zeros(1,19);
u15(6)=-0.1443;
u15(8)=0.4082;
u15(9)=-0.4771;
u15(11)=-0.1443;
u15(13)=-0.0397;
u15(15)=3/4;
u15=u15./norm(u15);
u16=zeros(1,19);
u16(6)=1/30;
u16(7)=-0.3300;
u16(8)=0.0471;
u16(9)=0.1102;
u16(10)=-14/15;
u16(11)=1/30;
u16(13)=0.0092;
u16(15)=0.0577;
u16=u16./norm(u16);

u17=zeros(1,19);
u17(8)=1/sqrt(3);
u17(11)=1/sqrt(6);
u17(15)=1/sqrt(2);
u17=u17./norm(u17);
u18=zeros(1,19);
u18(6)=-1/sqrt(3);
u18(9)=-1/sqrt(6);
u18(13)=1/sqrt(2);
u18=u18./norm(u18);
u19=zeros(1,19);
u19(19)=1;

Eq.changeBasis(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19);
%%
figure

plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolutionTr(5,5,:)+Eq.evolutionTr(6,6,:)+Eq.evolutionTr(7,7,:)+Eq.evolutionTr(8,8,:)+Eq.evolutionTr(9,9,:)+Eq.evolutionTr(10,10,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolutionTr(11,11,:)+Eq.evolutionTr(12,12,:)+Eq.evolutionTr(13,13,:)+Eq.evolutionTr(14,14,:)+Eq.evolutionTr(15,15,:)+Eq.evolutionTr(16,16,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolutionTr(17,17,:)+Eq.evolutionTr(18,18,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(1,1,:)+Eq.evolution(2,2,:)+Eq.evolution(3,3,:)+Eq.evolution(4,4,:)))
hold on
plot(Eq.evTime(1,:)./1000,squeeze(Eq.evolution(17,17,:)+Eq.evolution(18,18,:)+Eq.evolution(19,19,:)))
legend('Bright States','All Dark States','Disjoint Dark States','J=0','Je=1')
xlim([0 35])
xlabel('Time [\mu s]')
ylabel('Populations')
drawnow

%%
pts=squeeze(single(real(Eq.evolution(end-2,end-2,:)+Eq.evolution(end-1,end-1,:)+Eq.evolution(end,end,:))));
%%
disp(Eq.evTime(end))
tr=trapz(Eq.evTime(:),pts(:));% [ns]
rate=tr*Gamma/(Eq.evTime(end)); %
invrate_g=Gamma/rate;
ph_num=tr*Gamma;

disp(tr)
disp(rate)
disp(invrate_g)
disp(ph_num)

%%
figure
plot(Eq.evTime(1,:)./1000,pts)
xlabel('Time [\mu s]')
ylabel('Excited State Population')
drawnow