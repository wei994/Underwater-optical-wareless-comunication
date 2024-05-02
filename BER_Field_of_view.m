 tic %
clear all
%close all 

%SiPM Parameters-30020 series
PDE=0.24; %photo detection efficiency%3mm Sensor 20um - 30020, overvoltage @ 2.5V
lamda=450e-9;%Wavelength of Blue light 
A=0.0044;   %SiPM surface area
G=1e6;% Gain
pap=0.002; %probability of afterpulsing
pct=0.0003; %probability of crosstalk
dCR=6.6*10^6;  %dark count rate6.6*10^6
c=2.25e8;
hp=6.626e-34;% Planck's constant h,
e=1.6e-19;% Charge of Electron
Nspad=10998; %number of SPAD (microcells)
td= 100e-9; %tau_d deadtime

Pac=1+pap+pct;% probability of current
Res=PDE*lamda*G*e*Pac/(hp*c);% what? % SiPM responsitivity
Fe=1.1; %funtion of the doping profile



%Laser Diode Parameters  
    
ext=0.4;          %extinction ratio
Ptx1=0.1; %power of bit 1 in watt
Ptx0=ext*Ptx1; %0.0008  %power of bit 0 in watt
% Ptx = [0.0001,0.01,1]
ce=0.151; %diffuse coefficient of attenuation for clear waters 

%Rb1=[1:1:10].*1e6;
%Zb=[50 60 70 80 90 100 110 120 130 140 150];%Link range
%Z = 60;
%for b=1:length(Zb)%calculate the BER of each range
Z = 60;
figure
Esun_lamda = [0.02 0.2];
theta_FOVr0 = [4:2:28];
for i = 1:length(Esun_lamda)
for b = 1:length(theta_FOVr0)
Rb=1.*1e6;
% Rb = Rb1(b);% choose the bit rate

%Z=Zb(b);%%choose the range
%% path loss
Lch=exp(-ce.*Z);%channel loss%diffuse coefficient ce
h_L=Lch; %channel loss factor

%% Turbulence
lenth=1e7;%%what is the unit of measure

muuT=-0.3516;%% mean of log-amplitude coefﬁcient of turbulence
sigmaaT=0.7032;%%variance of log-amplitude coefﬁcient of turbulence
% SigmaaTT = [0.7032,0.9872,1.21]
% sigmaaT = SigmaaTT(i);
Norm_distT=makedist("Normal","mu",muuT,"sigma",sigmaaT);
gen_Norm_numbers_YT=random(Norm_distT,[1,lenth]);

h_T1=exp(gen_Norm_numbers_YT);
h_T=h_T1; %Turbulence factor,Its a distribution

%% pointing errors

sigma_rp = 1;%linspace(0.2,2,10);%(1:1:10).*1e-1;
 %%variance of radial displacement of the transmitted beam spot at L for the
 %center fo the Rx lens

Ray_dist = makedist('Rayleigh','B',sigma_rp);

r_p=random(Ray_dist,[1,lenth]); %


A_apr = A;%A; %pi.*r_apr.^2;Pi.*r^2 %Area of lens %SiPM surface area
r_apr = sqrt(A_apr./pi);%% radius of lens
D_apr = 2.*r_apr;%% diameters of lens

div_angle = 8.726646e-5;% divergence angle
w_z=div_angle.*Z;% approximation of beam width
%w_z = 3;

v_sm=(sqrt(pi).*r_apr)./(sqrt(2).*w_z);% variables

A_0 = (erf(v_sm)).^2; %factor at no pointing error

w_zeq2 = (w_z.^2.*sqrt(pi).*erf(v_sm))./(2.*v_sm.*exp(-v_sm.^2)); %equivalent beam width

h_P = A_0.*exp((-2.*r_p.^2)./w_zeq2); %pointing error facor

%% total channel gain

h_gain = h_L.*h_T; %.*h_P - no pointing error 
%

%% noise
    %%bandwith
Br=Rb/2;%1e7;% receiver low-pass filter Bandwidth / bandwidth of the receiver low-pass filter Be=1e9; %PD bandwidth (in Hz)(Electrical bandwidth)/ bandwidth of the receiver low-pass filter
RL=1000; %load resistor of TIA ohms
K=1.38064852e-23; % boltzman constant (m2 kg s-2 K-1)
Te=300;% equivqlent noise temperature [K]
sgmath2=4*K*Te*Br/RL; %thermal noise variance


Bo=2e-9; %Optical 'filter' bandwidth  at the Rx
%Esun_lamda = [0.02 0.2]./(1e-9);
%Esun_lamda_0 = Esun_lamda(i);
Esun_lamda_0=Esun_lamda(i)./(1e-9); %Spectral radiance of the background radiations  
tw=0.97; %water transmittance h_prime 0.95

Kd=0.08; %diffuse coefficient of attenuation for clear waters  
h_Lb=exp(-Kd.*Z);%channel loss

% muu = 0;
% theta_FOVr = 0.07;%4 deg. = 0.07rad,  1 rad = 60 deg
theta_FOVr = 0.07.*theta_FOVr0(b)./4;
Pb=Esun_lamda_0.*h_Lb.*tw.*Bo.*A_apr.*pi.*theta_FOVr.^2;% background power
% I_b = 0; %no solar noise %Res.*Pb.*RL;
I_b = Res.*Pb.*RL
% sigma2_b = 0; %no solar noise %2.*e.*G.*Fe.*Br.*I_b.*RL;
sigma2_b = 2.*e.*G.*Fe.*Br.*I_b.*RL;

%%%%%%%%%%%%%%Considering Optimal Threshold for an SiPM Output Photocurrent%%%%%%%%%%%%%%%
Is_1msm1=Res.*Ptx1.*RL;

I_d=dCR.*Pac.*G.*e.*RL;

Is_0msm0=ext.*Is_1msm1;


coeff=2*e*G*Fe*Br;

sgma2_d=coeff.*I_d;

sgma2_dth=sgma2_d.*RL  + sgmath2;


sgma2_1sm1=coeff.*Is_1msm1.*RL;

sgma2_0sm0=coeff.*Is_0msm0.*RL;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Threshold_opt = ((((Is_0msm0.*h_gain)+ I_d + I_b).*...
                (sgma2_1sm1.*h_gain + sgma2_dth + sigma2_b)-((Is_1msm1.*...
                h_gain)+ I_d + I_b).*(sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b))./...
                ((sgma2_1sm1.*h_gain + sgma2_dth + sigma2_b)-(sgma2_0sm0.*...
                h_gain + sgma2_dth + sigma2_b))+sqrt((((Is_0msm0.*h_gain+ I_d + I_b).*...
                (sgma2_1sm1.*h_gain+ sgma2_dth + sigma2_b)-(Is_1msm1.*h_gain+ I_d + I_b).*...
                (sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b))./((sgma2_1sm1.*...
                h_gain + sgma2_dth + sigma2_b)-(sgma2_0sm0.*h_gain +...
                sgma2_dth + sigma2_b))).^2+(((Is_1msm1.*h_gain)+ I_d + I_b).^2.*...
                (sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b)-...
                ((Is_0msm0.*h_gain)+ I_d + I_b).^2.*(sgma2_1sm1.*h_gain +...
                sgma2_dth + sigma2_b))./((sgma2_1sm1.*h_gain + sgma2_dth + sigma2_b)-...
                (sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b))-(((sgma2_0sm0.*...
                h_gain + sgma2_dth + sigma2_b).*(sgma2_1sm1.*h_gain +...
                sgma2_dth + sigma2_b))./((sgma2_1sm1.*h_gain + sgma2_dth + sigma2_b)-...
                (sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b))).*...
                (log((sgma2_0sm0.*h_gain + sgma2_dth + sigma2_b)./...
                (sgma2_1sm1.*h_gain + sgma2_dth + sigma2_b)))));
                       
Numr_trm1 = (Is_1msm1.*h_gain)+ I_d + I_b;%nurmical simulated current of bit 1 

Denm_trm1 =  sqrt(2*((sgma2_1sm1.*h_gain) + sgma2_dth + sigma2_b)); 

Numr_trm0 = (Is_0msm0.*h_gain)+ I_d + I_b;

Denm_trm0 =  sqrt(2*((sgma2_0sm0.*h_gain) + sgma2_dth + sigma2_b)); 

Avg_BER_Inst =  0.25.*erfc(((Threshold_opt - Numr_trm0)./Denm_trm0)) +...
                0.25.*erfc(((Numr_trm1 - Threshold_opt)./Denm_trm1));


Avg_BER1(b,:)= mean(Avg_BER_Inst);


end


%figure;
%semilogy(Rb1,Avg_BER1,'-ob')
semilogy(theta_FOVr0,Avg_BER1)
grid on
ylabel('{\itP}_{e}','Rotation',0');
xlabel('theta_F_O_V_r');
hold on
%legend('{\itR}_{b} = 1 Mbps','{\itR}_{b} = 10 Mbps','{\itR}_{b} = 100 Mbps', 'Location','east');
%axis([50 150 1e-10 1])
end

toc

%%
% legend("E_sun(0)=0.02","E_sun(0)=0.2")
% legend("L=60m ","L=90m","L=120m")