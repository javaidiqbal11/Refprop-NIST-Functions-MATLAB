function [T_w_cond_out]=CondTwo(T_cond,T_w_cond_in,m_dot_w_cond,Dtime)

Ref='water'; 
SecFld='water';
% Input geometrical characteristics
D_cond_o=0.75*0.0254;
t_cond_f=0.8E-3;
D_cond_i=D_cond_o-(2*t_cond_f);
N_t_cond=457;
L_t_cond=2654E-3;
N_pass_cond=3;
N_row_cond=11;
epslon_cond=0.0015E-3;
% Geometrical characteristics calculations
A_cond=pi*D_cond_o*L_t_cond*N_t_cond;
%==========================================================================
ItrC=0;
T_w_cond_out=T_w_cond_in;
for n = 1:10
    ItrC = ItrC+1;
    % A- refrigerant side convective resistance
    % A-1 Thermo-physical properties
    T_ref_cond_sat=T_cond;
    T_w_cond=(T_w_cond_in+T_w_cond_out)/2;
    T_ref_cond_surf=T_w_cond;
    T_ref_cond_v=(T_ref_cond_surf+T_ref_cond_sat)/2;
    Roh_ref_cond_L=refpropm('D','T',T_ref_cond_sat,'Q',0,Ref); %kg/m3
    Roh_ref_cond_v=refpropm('D','T',T_ref_cond_v,'Q',1,Ref); %kg/m3
    K_ref_cond_L=refpropm('L','T',T_ref_cond_sat,'Q',0,Ref)*1E-3; %kW/m K
    h_ref_cond_L=refpropm('H','T',T_ref_cond_sat,'Q',0,Ref)*1E-3; %kJ/kg
    h_ref_cond_v=refpropm('H','T',T_ref_cond_v,'Q',1,Ref)*1E-3; %kJ/kg
    h_ref_cond_fg=abs(h_ref_cond_v-h_ref_cond_L);
    Cp_ref_cond_v=refpropm('C','T',T_ref_cond_v,'Q',1,Ref)*1E-3; %kJ/kg K
    Meu_ref_cond_L=refpropm('V','T',T_ref_cond_sat,'Q',0,Ref); %Pa s
    % A-2 Calculations
    gr=9.81; C1_ref_cond=0.729;
    DeltaT_surf_sat=abs(T_ref_cond_sat-T_ref_cond_surf);
    h_ref_cond_fg_bar=h_ref_cond_fg+0.68*Cp_ref_cond_v*abs(T_ref_cond_surfT_ref_cond_sat);
    Num3=gr*Roh_ref_cond_L*(Roh_ref_cond_LRoh_ref_cond_v)*K_ref_cond_L^3*h_ref_cond_fg_bar;
    Den3=N_row_cond*Meu_ref_cond_L*DeltaT_surf_sat*D_cond_o;
    htc_ref_cond=C1_ref_cond*((Num3/Den3)^(1/4));
    R_ref_cond=1/(htc_ref_cond*A_cond);
    %--------------------------------------------------------------------------% A-3 Mass balance
    m_dot_ref_cond=htc_ref_cond*pi*D_cond_o*DeltaT_surf_sat/h_ref_cond_fg_bar;
    m_dot_ref_cond_tot=m_dot_ref_cond*L_t_cond*N_t_cond;
    m_ref_cond_tot_out=m_dot_ref_cond_tot*Dtime;
    %==========================================================================
    % B- Water side convective resistance
    % B-1 Thermo-physical properties
    P_atm=refpropm('P','T',373.15,'Q',0,'water');
    Roh_w_cond=refpropm('D','T',T_w_cond,'P',P_atm,SecFld);  %Kg/m3
    Meu_w_cond=refpropm('V','T',T_w_cond,'P',P_atm,SecFld);  %Pa s
    Cp_w_cond=refpropm('C','T',T_w_cond,'P',P_atm,SecFld)*1E-3; %kJ/kg.K
    K_w_cond=refpropm('L','T',T_w_cond,'P',P_atm,SecFld)*1E-3;  %kW/m.K
    Pr_w_cond=Cp_w_cond*Meu_w_cond/K_w_cond;
    % B-2 Calculations
    A_w_cond=pi*D_cond_i^2/4;
    N_t_cond_pass=N_t_cond/N_pass_cond;
    VeL_t_cond=m_dot_w_cond/(Roh_w_cond*A_w_cond*N_t_cond_pass);
    Re_w_cond=Roh_w_cond*VeL_t_cond*D_cond_i/Meu_w_cond;
    if(Re_w_cond<3000)
        f_cond_i = 64/Re_w_cond
        Nus_w_cond=4.36
    elseif(Re_w_cond>=3000)
        f_cond_i=(-1.8*log10((6.9/Re_w_cond)+(epslon_cond/(3.7*D_cond_i))^1.11))^-2;
        Num4=(Pr_w_cond*f_cond_i/8)*(Re_w_cond-1000);
        Den4=1+(12.7*((f_cond_i/8)^0.5)*((Pr_w_cond^(2/3))-1));
        Nus_w_cond=Num4/Den4;
    end
    htc_w_cond=Nus_w_cond*K_w_cond/D_cond_i;
    R_w_cond=1/(htc_w_cond*pi*D_cond_i*L_t_cond*N_t_cond);
    %==========================================================================
    % C- Tube wall conductance resistance
    K_t_cond=310*1E-3; %kW/m.K
    R_t_cond=(log(D_cond_o/D_cond_i))/(2*pi*K_t_cond*L_t_cond*N_t_cond);
    %==========================================================================
    % D- Overall thermal conductance
    UA_cond=1/(R_w_cond+R_t_cond+R_ref_cond);
    Twco=T_cond-((T_cond-T_w_cond_in)*exp(-UA_cond/(m_dot_w_cond*Cp_w_cond)));
    Tol_cond_w_out=abs(Twco-T_w_cond_out);
    if(Tol_cond_w_out<=1E-5)
        break
    else
        T_w_cond_out=Twco;
    end
end