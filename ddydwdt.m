function dy=ddydwdt(~,y,M_ref_L_cond,m_dot_HW,m_dot_CW,m_dot_CHW,THW_in,TCW_in,TCHW_in,FLAG1,FLAG3,Dtime,FLAG4,FLAG5,fin_pitch_bed_mm,Metal,Pge)

%==========================================================================
% This sub-program for solving energy balance equations for different heat
% exchangers (adsorbent bed, evaporator and condenser)
%Initial vector
dy=zeros(7,1);
T_bed2=y(2); 
T_cond=y(3); 
T_bed=y(5); 
T_evap=y(6); 
M_ref_L_evap=y(7);
%==========================================================================
FLAG2=1; Cv=1.9; fin_pitch_bed=fin_pitch_bed_mm/1000; 
[MC_bed, MC_SBed, MC_Ad, M_bed_ads]=MCs(fin_pitch_bed, Metal, Pge);
MCp_bed_met=MC_Ad+MC_bed+(616.816*3.83-616.816);
MCp_bed_ads=MC_SBed; MCp_bed_w=407.8*4.18;
MCp_cond_met=191.3277*1.25; MCp_evap_met=185.3718;
D_so=2.54E-4; Rp=0.15E-3; Ea=4.2E4; R=8.3145; H_ads=2.51E3;
%==========================================================================
%Desorber temperature constant
%==========================================================================
if((FLAG1==1)&&(FLAG3==1))
    T_w_bed_in2=THW_in; m_dot_water_bed2=m_dot_HW;
elseif((FLAG1==1)&&(FLAG3==0))
    T_w_bed_in2=THW_in; m_dot_water_bed2=m_dot_HW;
elseif((FLAG1==0)&&(FLAG3==1))
    T_w_bed_in2=TCW_in; m_dot_water_bed2=m_dot_CW;
end
%--------------------------------------------------------------------------
%Bed heating water stream
[T_w_bed_out2, UA_bed2]= BedTwo(T_bed2,T_w_bed_in2,m_dot_water_bed2,fin_pitch_bed,Metal,Pge)
%--------------------------------------------------------------------------
P_sat_ref_cond=refpropm('P','T',T_cond,'Q',0,'water');  %P_sat_ref
P_sat_bed2=refpropm('P','T',T_bed2,'Q',0,'water');  %P_sat_ads
%--------------------------------------------------------------------------
if((FLAG1==1)&(FLAG3==0))
    w_star_bed2=Uptake_sat(T_bed,T_bed2);
    W_bed_const1=15*D_so/Rp^2;
elseif((FLAG1==0)&(FLAG3==1))
    w_star_bed2=Uptake_sat(T_bed2,T_bed2);
    W_bed_const1=0;
else
    w_star_bed2=Uptake_sat(T_cond,T_bed2);
    W_bed_const1=15*D_so/Rp^2;
end
%--------------------------------------------------------------------------
Cp_w_T_bed2=refpropm('C','T',T_bed2,'Q',1,'water')*1E-3;  %Cp_w_T_bed
%--------------------------------------------------------------------------
T_bed_const1=FLAG1*M_bed_ads*H_ads*(15*D_so/Rp^2);
T_bed_const2=(1-FLAG4)*m_dot_water_bed2*4.18*(T_w_bed_in2-T_w_bed_out2);
T_bed_const3=MCp_bed_ads+MCp_bed_met+(FLAG4*MCp_bed_w);
T_bed_const4=M_bed_ads*Cp_w_T_bed2;
%==========================================================================
%Adsorber temperature constant
%==========================================================================
if((FLAG1==1)&&(FLAG3==1))
    T_w_bed_in=TCW_in; m_dot_water_bed=m_dot_CW;
elseif((FLAG1==1)&&(FLAG3==0))
    T_w_bed_in=TCW_in; m_dot_water_bed=m_dot_CW;
elseif((FLAG1==0)&&(FLAG3==1))
    T_w_bed_in=T_w_bed_out2; m_dot_water_bed=m_dot_CW;
end
%--------------------------------------------------------------------------
%Bed cooling water stream
[T_w_bed_out,UA_bed]=BedTwo(T_bed,T_w_bed_in,m_dot_water_bed,fin_pitch_bed,Metal,Pge)
%--------------------------------------------------------------------------
P_sat_ref_evap=refpropm('P','T',T_evap,'Q',0,'water');
P_sat_bed=refpropm('P','T',T_bed,'Q',0,'water');
%--------------------------------------------------------------------------
h_ref_g_T_hex=refpropm('H','T',T_evap,'Q',1,'water')*1E-3;
h_ref_P_hex_T_bed=refpropm('H','T',T_bed,'P',P_sat_bed,'water')*1E-3;
Cp_w_T_bed=refpropm('C','T',T_bed,'Q',1,'water')*1E-3;
%--------------------------------------------------------------------------
if((FLAG1==1)&&(FLAG3==0))
    w_star_bed=Uptake_sat(T_evap,T_bed);
    W_bed_const2=0;
elseif((FLAG1==0)&&(FLAG3==1))
    w_star_bed=Uptake_sat(T_bed,T_bed);
    W_bed_const2=0;
else
    w_star_bed=Uptake_sat(T_evap,T_bed);
    W_bed_const2=15*D_so/Rp^2;
end
%--------------------------------------------------------------------------
T_bed_const5=FLAG1*M_bed_ads*H_ads*(15*D_so/Rp^2);
T_bed_const6=FLAG1*M_bed_ads*((FLAG3*(h_ref_g_T_hex-h_ref_P_hex_T_bed))...
+((1-FLAG3)*Cv*(T_bed2-T_bed)))*(15*D_so/Rp^2);
T_bed_const7=(1-FLAG4)*m_dot_water_bed*4.18*(T_w_bed_in-T_w_bed_out);
T_bed_const8=MCp_bed_ads+MCp_bed_met+(FLAG4*MCp_bed_w);
T_bed_const9=M_bed_ads*Cp_w_T_bed;
%==========================================================================
%Condenser constants
%==========================================================================
if((FLAG1==1)&&(FLAG3==1))
    T_w_cond_in=T_w_bed_out; m_dot_water_cond=m_dot_CW;
elseif((FLAG1==1)&&(FLAG3==0))
    T_w_cond_in=TCW_in; m_dot_water_cond=m_dot_CW;
elseif((FLAG1==0)&&(FLAG3==1))
    T_w_cond_in=T_w_bed_out; m_dot_water_cond=m_dot_CW;
end
%--------------------------------------------------------------------------
%Condenser water outlet temperature
[T_w_cond_out]=CondTwo(T_cond,T_w_cond_in,m_dot_water_cond,Dtime);
%--------------------------------------------------------------------------
h_ref_cond_out=refpropm('H','T',T_cond,'Q',0,'water')*1E-3;
h_ref_bed2_out=refpropm('H','T',T_bed2,'P',P_sat_bed2,'water')*1E-3;
Cp_ref_L_cond=refpropm('C','T',T_cond,'Q',0,'water')*1E-3;
h_ref_cond_in=refpropm('H','T',T_cond,'Q',1,'water')*1E-3;
C_ref_bed2_out=refpropm('C','T',T_bed2,'Q',1,'water')*1E-3;
%--------------------------------------------------------------------------
T_cond_const1=-FLAG5*M_bed_ads*(h_ref_cond_inh_ref_cond_out)*(15*D_so/Rp^2);
T_cond_const2=m_dot_water_cond*4.18*(T_w_cond_in-T_w_cond_out);
T_cond_const3=(M_ref_L_cond*Cp_ref_L_cond)+MCp_cond_met;
T_cond_const4=-FLAG5*M_bed_ads*C_ref_bed2_out*(T_bed2-T_cond)*(15*D_so/Rp^2);
%==========================================================================
%Evaporator constants
%==========================================================================
T_w_evap_in=TCHW_in;
m_dot_water_evap=m_dot_CHW;
%--------------------------------------------------------------------------
%Evaporator cooling water stream
[T_w_evap_out]=EvapTwo(T_evap,T_w_evap_in,m_dot_water_evap);
%--------------------------------------------------------------------------
h_ref_evap_in=refpropm('H','T',T_evap,'Q',0,'water')*1E-3;
h_ref_evap_out=refpropm('H','T',T_evap,'Q',1,'water')*1E-3;
Cp_ref_L_evap=refpropm('C','T',T_evap,'Q',0,'water')*1E-3;
%--------------------------------------------------------------------------
T_evap_const1=FLAG5*M_bed_ads*(h_ref_evap_inh_ref_evap_out)*(15*D_so/Rp^2);
T_evap_const2=m_dot_water_evap*4.18*(T_w_evap_in-T_w_evap_out);
T_evap_const3=(M_ref_L_evap*Cp_ref_L_evap)+MCp_evap_met;
%==========================================================================
dy(1)=W_bed_const1*(w_star_bed2-y(1))*exp(-Ea/(R*y(2)));
dy(2)=(T_bed_const1*exp(-Ea/(R*y(2)))*(w_star_bed2-y(1))+T_bed_const2)/(T_bed_const3+(y(1)*T_bed_const4));
dy(3)=((T_cond_const1+T_cond_const4)*exp(-Ea/(R*y(2)))*(w_star_bed2-y(1))+T_cond_const2)/T_cond_const3;
%==========================================================================
dy(4)=(FLAG3*W_bed_const2*exp(-Ea/(R*y(5)))*(w_star_bed-y(4)))-((1-FLAG3)*W_bed_const1*(w_star_bed2-y(1))*exp(-Ea/(R*y(2))));
dy(5)=((T_bed_const5+T_bed_const6)*exp(-Ea/(R*y(5)))*(w_star_bedy(4))+T_bed_const7)/(T_bed_const8+(T_bed_const9*y(4)));
dy(6)=(T_evap_const1*exp(-Ea/(R*y(5)))*(w_star_bedy(4))+T_evap_const2+(0.04*(h_ref_evap_in-h_ref_evap_out)))/T_evap_const3;
%==========================================================================
dy(7)=-M_bed_ads*FLAG5*((W_bed_const1*(w_star_bed2-y(1))*exp(-Ea/(R*y(2))))+(W_bed_const2*exp(-Ea/(R*y(5)))*(w_star_bed-y(4))));
end

