function [T_w_bed_out,UA]= BedTwo(T_bed,T_w_bed_in,m_dot_water_bed,fin_pitch_bed,Metal,Pge) 

Ref='water'; 
SecFld='water';
%==========================================================================
% Input geometrical characteristics  L_fin_bed=340E-3; H_fin_bed=0.0285;   %fin_pitch_bed=1.5E-3;
Dp=0.31; L_module=3400E-3; 
N_t_module=12; N_pass_bed=2;
N_module=7*2*4; D_bed_o=(5/8)*0.0254; t_bed=0.8E-3; D_bed_i=D_bed_o-(2*t_bed);
W_fin_bed=0.105E-3; epslon_bed=0.0015E-3; 
N_fin_module=round(L_module/fin_pitch_bed);
%==========================================================================
% Array booking
N_inc=N_fin_module*2; 
site=zeros(1,N_inc);Segma_R=zeros(1,N_inc);UA_bed=zeros(1,N_inc); T_w_out=zeros(1,N_inc);
DT_b=zeros(1,N_inc);DT_s=zeros(1,N_inc);LMTD=zeros(1,N_inc);dq=zeros(1,N_inc);
T_w_in=zeros(1,N_inc);T_w_out=zeros(1,N_inc);htc_w_bedi=zeros(1,N_inc);
%==========================================================================
% Geometrical calculation
a=L_fin_bed/N_t_module; 
b=H_fin_bed; 
r_fin_bed=sqrt(a*b/pi); 
D_fin_bed=r_fin_bed*2;
fin_space_bed=fin_pitch_bed-W_fin_bed; A_w_bed=pi*D_bed_i^2/4;
N_t_bed=N_module*N_t_module; N_t_bed_pass=N_t_bed/N_pass_bed;
T_bed_C=T_bed-273; R_cont_TSG=Rc(T_bed_C,Dp); R_cont_FSG=Rc(T_bed_C,Dp);
%==========================================================================
% Area calculation
L_FT=W_fin_bed*N_fin_module; L_UFT=L_module-L_FT; A_UFT=pi*D_bed_o*L_UFT; 
A_UFTM=A_UFT*N_t_module;
A_FS=(L_fin_bed*H_fin_bed)-(0.25*pi*D_bed_o^2*N_t_module); 
A_FTP=(2*L_fin_bed*W_fin_bed)+(2*H_fin_bed*W_fin_bed);
A_F=A_FS+A_FTP; A_FM=A_F*N_fin_module;
A_M=A_FM+A_UFTM; A_bed=A_M*N_module;
%==========================================================================
% Mass calculation
Roh_F=3661.85; C_F=0.896; Roh_T=8954; C_T=0.3831; Roh_S=690.987; C_S=0.921; 
%Roh_S=708.299;
M_TM=0.25*pi*(D_bed_o^2-D_bed_i^2)*L_module*N_t_module*Roh_T;
M_FM=A_FS*W_fin_bed*Roh_F*N_fin_module;
MC_M=M_TM*C_T+M_FM*C_F;
%==========================================================================
V_TM=0.25*pi*L_module*N_t_module*D_bed_o^2;
V_FM=A_FS*W_fin_bed*N_fin_module;
V_M=L_module*L_fin_bed*H_fin_bed;
V_SM=V_M-V_FM-V_TM; M_SM=V_SM*Roh_S; M_SBed=(M_SM*N_module);
%==========================================================================
[K_SG, C_Ad]=KC_mix(Metal, Pge);
M_bed_ads=M_SBed*(1-Pge/100); MC_SBed=M_SBed*C_S;
MC_Ad=(M_SBed*Pge/100)*C_Ad; MC_bed=MC_M*N_module;
%==========================================================================
% Outside surface heat transfer resistance
A_s_fin_bed=(a*b)-(pi*D_bed_o^2/4);
d_SG=fin_space_bed/2;
RA1=R_cont_FSG/(2*A_s_fin_bed);
% Axial contact thermal resistance(Per increament)
RA2=d_SG/(2*A_s_fin_bed*K_SG);% Axial thermal resistance (Per increament)
RA=RA1+RA2; % Total axial thermal resistance (Per increament)
RB1=R_cont_TSG/(2*pi*D_bed_o*d_SG); % Radial Contact thermal resistance (Per increment)
RB2=(log(D_fin_bed/D_bed_o))/(4*pi*K_SG*d_SG); % Radial thermal resistance (Per increament)
RB=RB1+RB2; % Total radial thermal resistance (Per increment)
R_bed=RA*RB/(RA+RB); % Bed side thermal resistance for one increment
%==========================================================================
% Tube wall heat heat transfer resistance
K_t_bed=310*1E-3; %kW/m.K
R_t_bed=(log(D_bed_o/D_bed_i))/(2*pi*K_t_bed*fin_pitch_bed);
%==========================================================================
% Start of calculation loop
T_w_in(1)=T_w_bed_in;

for i=1:N_inc
    site(i)=i;
    %==========================================================================
    % Water side heat transfer resistance
    % Thermo-physical properties
    T_w_bed=T_w_in(i);
    P_atm=refpropm('P','T',373.15,'Q',0,SecFld);
    Roh_w_bed=refpropm('D','T',T_w_bed,'P',P_atm,SecFld); %Kg/m3
    Meu_w_bed=refpropm('V','T',T_w_bed,'P',P_atm,SecFld); %Pa s
    Cp_w_bed=refpropm('C','T',T_w_bed,'P',P_atm,SecFld)*1E-3; %kJ/kg.K
    K_w_bed=refpropm('L','T',T_w_bed,'P',P_atm,SecFld)*1E-3; %kW/m.K
    Pr_w_bed=Cp_w_bed*Meu_w_bed/K_w_bed;
    % Calculations
    VeL_t_bed=m_dot_water_bed/(Roh_w_bed*A_w_bed*N_t_bed_pass);
    Re_w_bed=Roh_w_bed*VeL_t_bed*D_bed_i/Meu_w_bed;
    if(Re_w_bed<3000)
        f_bed_i=64/Re_w_bed;
        Nus_w_bed=4.36;
    elseif(Re_w_bed>=3000)
        f_bed_i=(-1.8*log10((6.9/Re_w_bed)+(epslon_bed/(3.7*D_bed_i))^1.11))^-2;
        Num5=(Pr_w_bed*f_bed_i/8)*(Re_w_bed-1000);
        Den5=1+(12.7*((f_bed_i/8)^0.5)*((Pr_w_bed^(2/3))-1));
        Nus_w_bed=Num5/Den5;
    end
    htc_w_bed=Nus_w_bed*K_w_bed/D_bed_i;
    R_w_bed=1/(htc_w_bed*pi*D_bed_i*fin_pitch_bed);
    %==========================================================================
    % Thermal balance and outlet temperature calculation
    Segma_R(i)=(R_w_bed+R_t_bed+R_bed);
    UA_bed(i)=1/Segma_R(i);
    T_w_out(i)=T_bed+(T_w_in(i)-T_bed)*(exp(-UA_bed(i)/(m_dot_water_bed*4.18/N_t_bed_pass)));
    if(T_w_bed_in>T_bed)
        DT_b(i)=T_w_in(i)-T_bed; DT_s(i)=T_w_out(i)-T_bed;
    else
        DT_b(i)=T_bed-T_w_in(i); DT_s(i)=T_bed-T_w_out(i);
    end
    LMTD(i)=(DT_b(i)-DT_s(i))/(log(DT_b(i)/DT_s(i)));
    dq(i)=UA_bed(i)*LMTD(i)*(12*7*4);
    %==========================================================================
    htc_w_bedi(i)=htc_w_bed;
    %==========================================================================
    if(i==N_inc)
    
    else
        T_w_in(i+1)=T_w_out(i);
    end
end
T_w_bed_out=T_w_out(N_inc);
%==========================================================================
% Chick the tolerance in calculations
Q_bed1=sum(dq);
Q_bed2=m_dot_water_bed*4.18*(abs(T_w_bed_in-T_w_bed_out));
Tolerance=abs(Q_bed1-Q_bed2);
%==========================================================================
if(T_w_bed_in>T_bed)
    DT_b_bed=T_w_bed_in-T_bed; DT_s_bed=T_w_bed_out-T_bed;
else
    DT_b_bed=T_bed-T_w_bed_in; DT_s_bed=T_bed-T_w_bed_out;
end
LMTD_bed=(DT_b_bed-DT_s_bed)/(log(DT_b_bed/DT_s_bed));
UA=Q_bed1/(LMTD_bed);
U_bed=Q_bed1/(A_bed*LMTD_bed);
end
