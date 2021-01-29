%This is the home screen of evaluating the performance of the adsorption
%chiller ADCM1-180
%==========================================================================
clear all; 
close all; 
clc;
% Input data
fin_pitch_bed_mm=1.1;  % Minimum fin spacing
Metal='Al'; 
Pge=15;  % Type of metal additives and the percentage
for i=1:19
    fin_pitch_bed_mm=fin_pitch_bed_mm+0.1;
    fin_pitch_bed=fin_pitch_bed_mm/1000;
    N_cyc=8; 
    Dtime=5;
    t_normal=430; t_regen=30; t_htrec=20;
    t_cycle=(t_normal+t_regen+t_htrec);
    [MC_bed, MC_SBed, MC_Ad, M_bed_ads]=MCs(fin_pitch_bed, Metal, Pge);
    M_ref_tot=185; M_ref_L_cond_MAX=50; 
    THW_in=(88.6+273); TCW_in=(29.5+273); TCHW_in=(11.1+273);
    m_dot_HW=18.3; m_dot_CW=66.6; m_dot_CHW=19.7;
    %==========================================================================
    % Initial conditions
    Time_i=0;
    t=Time_i;
    T_bed_i=27+273; 
    T_bed2_i=27+273; 
    T_evap_i=27+273; 
    T_cond_i=27+273; 
    W_bed2_i=0.05;
    W_bed_i=0.05; 
    M_ref_L_cond=50*0.9;
    M_ref_L_evap_i=M_ref_tot;
    %==========================================================================
    step=Dtime; S=0;
    %--------------------------------------------------------------------------
    for M=1:N_cyc 
        Q_in=0;
        Q_evap=0;
        t_round=0;
        %--------------------------------------------------------------------------
        for L=1:step:t_cycle
            S=S+1;
            t=t+step;
            Time(S)=t;
            %--------------------------------------------------------------------------
            t_round=t_round+step;
            if(t_round<=t_normal)
                FLAG1=1; FLAG3=1; FLAG4=0; FLAG5=1;
            elseif((t_round>t_normal)&&(t_round<=(t_normal+t_regen)))
                FLAG1=1; FLAG3=0; FLAG4=1; FLAG5=0;
            elseif((t_round>(t_normal+t_regen))&&(t_round<=(t_normal+t_regen+t_htrec)))
                FLAG1=0; FLAG3=1; FLAG4=0; FLAG5=0;
            end
    
            %==========================================================================
            Timerange2=[Time_i t];
            Initialbed2=[W_bed2_i T_bed2_i T_cond_i W_bed_i T_bed_i T_evap_i 
            M_ref_L_evap_i];
            %==========================================================================
            option2=odeset('RelTol',1E-4,'AbsTol',1E-4); 

            Y=ode45(@ddydwdt,Timerange2,Initialbed2,option2,M_ref_L_cond,m_dot_HW,m_dot,CW,m_dot_CHW,THW_in,TCW_in,TCHW_in,FLAG1,FLAG3,Dtime,FLAG4,FLAG5,fin_pitch, bed_mm,Metal,Pge);
            %-----------------Hot Bed Parameters---------------------------------------
            Y2_t=deval(Y,t);
            W_bed2(S)=Y2_t(1);
            T_bed2(S)=Y2_t(2);
            T_cond(S)=Y2_t(3);
            %-----------------Cold Bed Parameters--------------------------------------
            W_bed(S)=Y2_t(4);
            T_bed(S)=Y2_t(5);
            T_evap(S)=Y2_t(6);
            %-----------------Evaporator Parameters------------------------------------
            M_ref_L_evap(S)=Y2_t(7);
            %==========================================================================
            TB=T_bed(S); TB2=T_bed2(S); TE=T_evap(S); TC=T_cond(S);
            %--------------------------------------------------------------------------
            if((FLAG1==1)&&(FLAG3==1))
                T_w_bed_in2=THW_in; m_dot_water_bed2=m_dot_HW;
                T_w_bed_out2=BedTwo(TB2,T_w_bed_in2,m_dot_water_bed2,fin_pitch_bed,Metal,Pge);
                T_w_bed_in=TCW_in; m_dot_water_bed=m_dot_CW;
                T_w_bed_out=BedTwo(TB,T_w_bed_in,m_dot_water_bed,fin_pitch_bed,Metal,Pge);
            elseif((FLAG1==0)&&(FLAG3==1))
                T_w_bed_in2=TCW_in; m_dot_water_bed2=m_dot_CW;
                T_w_bed_out2=BedTwo(TB2,T_w_bed_in2,m_dot_water_bed2,fin_pitch_bed,Metal,Pge);
                T_w_bed_in=T_w_bed_out2; m_dot_water_bed=m_dot_CW;
                T_w_bed_out=BedTwo(TB,T_w_bed_in,m_dot_water_bed,fin_pitch_bed,Metal,Pge);
            end
            %-----------------Cooling Water Stream Temp--------------------------------
            if((FLAG1==1)&&(FLAG3==1))
                T_w_cond_in=T_w_bed_out; m_dot_water_cond=m_dot_CW;
            elseif((FLAG1==1)&&(FLAG3==0))
                T_w_cond_in=TCW_in; m_dot_water_cond=m_dot_CW;
            elseif((FLAG1==0)&&(FLAG3==1))
                T_w_cond_in=T_w_bed_out; m_dot_water_cond=m_dot_CW;
            end
            TCW_inn(S)=TCW_in; 
            TCW_out(S)=CondTwo(TC,T_w_cond_in,m_dot_water_cond,Dtime);
            %-----------------Hot Water Steam Temp-------------------------------------
            if((FLAG1==1)&&(FLAG3==1))
                THW_out(S)=T_w_bed_out2; THW_inn(S)=THW_in; 
            elseif((FLAG1==1)&&(FLAG3==0))
                THW_out(S)=THW_in; THW_inn(S)=THW_in;
            elseif((FLAG1==0)&&(FLAG3==1))
                THW_out(S)=THW_in; THW_inn(S)=THW_in;
            end
            %-----------------Chilled Water Stream Temp--------------------------------
            T_w_evap_in=TCHW_in; m_dot_water_evap=m_dot_CHW;
            TCHW_inn(S)=TCHW_in; TCHW_out(S)=EvapTwo(TE,T_w_evap_in,m_dot_water_evap); 
            Chilled_out=TCHW_out(S)-273;
            %==========================================================================
            % Adsorption / Desorption Switching 

            if ((M==1)||(M==3)||(M==5)||(M==7)||(M==9)||(M==11)||(M==13)||(M==15)||(M==17)||(M==19))
                T_des(S)=T_bed2(S); W_des(S)=W_bed2(S); T_ads(S)=T_bed(S); 
                W_ads(S)=W_bed(S);
            elseif((M==2)||(M==4)||(M==6)||(M==8)||(M==10)||(M==12)||(M==14)||(M==16)||(M==18)||(M==20))
                T_des(S)=T_bed(S); W_des(S)=W_bed(S); T_ads(S)=T_bed2(S); 
                W_ads(S)=W_bed2(S);
            end
            %==========================================================================
            m_water_ads(S)=W_bed(S)*M_bed_ads;
            m_water_des(S)=W_bed2(S)*M_bed_ads;
            m_water_tot(S)=m_water_ads(S)+m_water_des(S);
            %==========================================================================
            % Instant Performance Indicators
            dq_in=(m_dot_HW*4.18*Dtime/t_cycle)*(THW_inn(S)-THW_out(S)); 
            Q_in=Q_in+dq_in;
            dq_evap=(m_dot_water_evap*4.18*Dtime/t_cycle)*(TCHW_inn(S)-TCHW_out(S)); 
            Q_evap=Q_evap+dq_evap;
            %------------------------------Data Plot-----------------------------------
            subplot(3,2,1); PLT1=plot(Time,THW_inn,'k--',Time,TCW_inn,'k--',Time,TCHW_inn,...
            'k--',Time,T_des,'b-',Time,T_ads,'c-',Time,T_evap,'m-',Time,T_cond,'r-');
            set(PLT1,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]','fontsize',...
            12); ylabel('Temperature [C]','fontsize',12); title('Heat Exchangers Temperature Profile','fontsize',14);
            subplot(3,2,2); PLT2=plot(Time,THW_inn,'r-',Time,THW_out,'m-',Time,TCW_inn,...
            'k-',Time,TCW_out,'g-',Time,TCHW_inn,'b-',Time,TCHW_out,'c-'); 
            set(PLT2,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]','fontsize',...
            12); ylabel('Temperatures [C]','fontsize',12); title('Heat Exchangers Outlet Temperature','fontsize',14);
            subplot(3,2,3); PLT3=plot(Time,M_ref_L_evap,'k-'); 
            set(PLT3,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]','fontsize',...
            12); ylabel('Evaporator Refrigerant [kg]','fontsize',12); 
            title('Evaporator Refrigerant mass','fontsize',14);
            subplot(3,2,4); PLT3=plot(Time,m_water_tot,'k-'); 
            set(PLT3,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]', 'fontsize',...
            12); ylabel('Ads & Des Refrigerant [kg]','fontsize',12); title('Ads & Des Refrigerant mass','fontsize',14);
            subplot(3,2,5); PLT3=plot(Time,m_water_ads,'k-'); 
            set(PLT3,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]', 'fontsize',...
            12); ylabel('Adsorber Refrigerant [kg]','fontsize',12); title('Adsorber Refrigerant mass','fontsize',14);
            subplot(3,2,6); PLT3=plot(Time,m_water_des,'k-');  

            set(PLT3,'linewidth',2); set(gca,'fontsize',10); xlabel('Time [Sec]','fontsize',...
            12); ylabel('Desorber Refrigerant [kg]','fontsize',12); title('Desorber Refrigerant mass','fontsize',14);
            drawnow;
            %==========================================================================
            Time_i=t;
            W_bed2_i=Y2_t(1);
            T_bed2_i=Y2_t(2);
            T_cond_i=Y2_t(3);
            %--------------------------------------------------------------------------
            W_bed_i=Y2_t(4);
            T_bed_i=Y2_t(5);
            T_evap_i=Y2_t(6);
            %--------------------------------------------------------------------------
            M_ref_L_evap_i=Y2_t(7);
        end
        % Cyclic Performance Indicators
        Q_Cooling(M)=Q_evap;
        Q_Heating(M)=Q_in;
        COP(M)=Q_evap/Q_in;
        SCP(M)=Q_evap/(M_bed_ads);
        %--------------------------------------------------------------------------
        W_bed2_i=Y2_t(4);
        T_bed2_i=Y2_t(5);
        %--------------------------------------------------------------------------
        W_bed_i=Y2_t(1);
        T_bed_i=Y2_t(2);
    end 
    %--------------------------------------------------------------------------
    % Complete Run Performance Indicators
    [MC_bedi(i),MC_SBedi(i),MC_Adi(i),M_bed_adsi(i)]= MCs(fin_pitch_bed,Metal,Pge);
    MCSM(i)=MC_SBedi(i)/(MC_bedi(i)+MC_Adi(i));
    P(i)=fin_pitch_bed_mm;
    Q_Coolingi(i)=Q_Cooling(N_cyc);
    Q_Heatingi(i)=Q_Heating(N_cyc);
    COPi(i)=COP(N_cyc);
    SCPi(i)=SCP(N_cyc);
end

