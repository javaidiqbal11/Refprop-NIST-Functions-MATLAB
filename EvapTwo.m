function [T_w_evap_out]=EvapTwo(T_evap, T_w_evap_in, m_dot_w_evap)

Ref='water'; 
SecFld='water';
% Input geometrical characteristics
D_evap_i=16.33E-3;
D_evap_f=18.85E-3;
D_evap_r=17.75E-3;
t_evap_f=0.20E-3;
Pt_evap=0.45E-3;
L_t_evap=3894e-3;
N_t_evap=224;
N_pass_evap=3;
Eta_evap_f=0.98;
epslon_evap=0.0015E-3;
%==========================================================================
% Geometrical characteristics calculations
r1_evap=D_evap_r/2;
r2_evap=D_evap_f/2;
r2_c_evap=r2_evap+(t_evap_f/2);
A_f_evap=2*pi*(r2_c_evap^2-r1_evap^2);
A_uf_evap=pi*D_evap_r*(Pt_evap-t_evap_f);
A_fs_evap=2*(pi/4)*(D_evap_f^2-D_evap_r^2);
A_ft_evap=pi*D_evap_f*t_evap_f;
A_evap_seg=A_f_evap+A_uf_evap;
Term1_evap=1.3*Eta_evap_f;
Term2_evap=A_fs_evap/(A_evap_seg*(0.25*pi*(D_evap_f^2-D_evap_r^2)/D_evap_f)^0.25);
Term3_evap=A_ft_evap/(A_evap_seg*D_evap_f^0.25);
Term4_evap=A_uf_evap/(A_evap_seg*D_evap_r^0.25);
D_evap=(Term1_evap*(Term2_evap+Term3_evap+Term4_evap))^-4;
N_evap_seg=L_t_evap/Pt_evap;
A_evap=A_evap_seg*N_evap_seg*N_t_evap;
%==========================================================================
itrE=0;
T_w_evap_out=T_w_evap_in;
for n=1:10
    itrE=itrE+1;
    % A- Refrigerant side convective Resistance
    % A-1 Thermo-physical properties calculations
    T_w_evap=(T_w_evap_out+T_w_evap_in)/2;
    T_ref_evap_surf=T_w_evap;
    T_ref_evap_sat=T_evap;
    T_ref_evap_v=(T_ref_evap_surf+T_ref_evap_sat)/2;
    Roh_ref_evap_L=refpropm('D','T',T_ref_evap_sat,'Q',0,Ref); %Kg/m3
    Roh_ref_evap_v=refpropm('D','T',T_ref_evap_v,'Q',1,Ref); %Kg/m3
    h_ref_evap_L=refpropm('H','T',T_ref_evap_sat,'Q',0,Ref)*1E-3; %kJ/kg
    h_ref_evap_v=refpropm('H','T',T_ref_evap_v,'Q',1,Ref)*1E-3; %kJ/kg
    h_ref_evap_fg=abs(h_ref_evap_v-h_ref_evap_L); %kJ/kg
    Cp_ref_evap_v=refpropm('C','T',T_ref_evap_v,'Q',1,Ref)*1E-3; %kJ/kg.K
    Cp_ref_evap_L=refpropm('C','T',T_ref_evap_sat,'Q',0,Ref)*1E-3; %kJ/kg.K
    Meu_ref_evap_v=refpropm('V','T',T_ref_evap_v,'Q',1,Ref); %Pa.s
    Meu_ref_evap_L=refpropm('V','T',T_ref_evap_sat,'Q',0,Ref); %Pa.s
    Neu_ref_evap_v=Meu_ref_evap_v/Roh_ref_evap_v; 
    K_ref_evap_v=refpropm('L','T',T_ref_evap_v,'Q',1,Ref)*1E-3; %kW/m.K
    K_ref_evap_L=refpropm('L','T',T_ref_evap_sat,'Q',0,Ref)*1E-3; %kW/m.K
    %--------------------------------------------------------------------------
    % A-2-1 Calculations 'Film Boiling'
    gr=9.81; C1_ref_evap=0.62;
    DeltaT_surf_sat=abs(T_ref_evap_surf-T_ref_evap_sat);
    h_ref_evap_fg_bar=h_ref_evap_fg+0.8*Cp_ref_evap_v*(T_ref_evap_surfT_ref_evap_sat);
    Num1=gr*(Roh_ref_evap_L-Roh_ref_evap_v)*h_ref_evap_fg_bar*D_evap^3;
    Den1=Neu_ref_evap_v*K_ref_evap_v*DeltaT_surf_sat;
    Nus_ref_evap=C1_ref_evap*((Num1/Den1)^0.25);
    htc_ref_evap1=Nus_ref_evap*K_ref_evap_v/D_evap;
    %--------------------------------------------------------------------------
    % A-2-2 Calculations 'Nuclate Boiling'
    C_s_f_evap=0.0068; 
    n=1;
    Ja_ref_evap=Cp_ref_evap_L*DeltaT_surf_sat/h_ref_evap_fg;
    Pr_ref_evap_L=Cp_ref_evap_L*Meu_ref_evap_L/Meu_ref_evap_L;
    TERM_A=Meu_ref_evap_L*h_ref_evap_fg/DeltaT_surf_sat;
    TERM_B=gr*(Roh_ref_evap_L-Roh_ref_evap_v)/72.81E-3;
    TERM_C=Ja_ref_evap/(C_s_f_evap*Pr_ref_evap_L^n);
    htc_ref_evap=TERM_A*(TERM_B^0.5)*(TERM_C^3);
    R_ref_evap=1/(htc_ref_evap*A_evap);
    %==========================================================================
    % B- Water side convective resistance
    % B-1 Thermo-physical properties
    P_atm=refpropm('P','T',373.15,'Q',0,'water');
    Roh_w_evap=refpropm('D','T',T_w_evap,'P',P_atm,SecFld); %Kg/m3
    Meu_w_evap=refpropm('V','T',T_w_evap,'P',P_atm,SecFld); %Pa s
    Cp_w_evap=refpropm('C','T',T_w_evap,'P',P_atm,SecFld)*1E-3; %kJ/kg.K
    K_w_evap=refpropm('L','T',T_w_evap,'P',P_atm,SecFld)*1E-3; %kW/m.K
    Pr_w_evap=Cp_w_evap*Meu_w_evap/K_w_evap;
    % B-2 Calculations
    A_w_evap=pi*D_evap_i^2/4;
    N_t_evap_pass=N_t_evap/N_pass_evap;
    VeL_t_evap=m_dot_w_evap/(Roh_w_evap*A_w_evap*N_t_evap_pass);
    Re_w_evap=Roh_w_evap*VeL_t_evap*D_evap_i/Meu_w_evap;
    
    if(Re_w_evap<3000)
        f_evap_i=64/Re_w_evap;
        Nus_w_evap=4.36;
    elseif(Re_w_evap>=3000)
        f_evap_i=(-1.8*log10((6.9/Re_w_evap)+(epslon_evap/(3.7*D_evap_i))^1.11))^-2;
        %f_evap_i=((0.790*log(Re_w_evap))-1.64)^-2
        Num2=(Pr_w_evap*f_evap_i/8)*(Re_w_evap-1000);
        Den2=1+(12.7*((f_evap_i/8)^0.5)*((Pr_w_evap^(2/3))-1));
        Nus_w_evap=Num2/Den2;
    end
    htc_w_evap=Nus_w_evap*K_w_evap/D_evap_i;
    R_w_evap=1/(htc_w_evap*pi*D_evap_i*L_t_evap*N_t_evap);
    %==========================================================================
    % C- Tube wall conductance resistance
    K_t_evap=310*1E-3; %kW/m.K
    D_evap_1=A_evap_seg/(pi*Pt_evap);
    R_t_evap=(log(D_evap_1/D_evap_i))/(2*pi*K_t_evap*L_t_evap*N_t_evap);
    %==========================================================================
    % D- Overall thermal conductance
    UA_evap=1/(R_w_evap+R_t_evap+R_ref_evap);
    Tweo=T_evap+(T_w_evap_in-T_evap)*exp(-UA_evap/(m_dot_w_evap*Cp_w_evap));
    Tol_evap_w_out=abs(Tweo-T_w_evap_out);
    if(Tol_evap_w_out<=1E-5)
        break
    else
        T_w_evap_out=Tweo;
    end
end
T_w_evap_out=Tweo;
end
