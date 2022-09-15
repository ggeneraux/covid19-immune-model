function rates =  calculate_rates(s,data_dictionary)
  	

  	% unpack plasma and lung volume
  	vol_plasma = data_dictionary.vol_plasma;
  	vol_alv_ml = data_dictionary.vol_alv_ml;

	p = data_dictionary.parameters;

    IC.AT1 = data_dictionary.initial_condition(2);
    IC.AT2 = data_dictionary.initial_condition(3);

	kdc = 1;

	% unpack parameters
	A_V = p.A_V;
	k_int = p.k_int;
	km_int_IFNb = p.km_int_IFNb;
	k_V_innate = p.k_V_innate;
	b_V = p.b_V;
	f_int = p.f_int;
	mu_AT2 = p.mu_AT2;
	k_ROS_AT2 = p.k_ROS_AT2;
	km_ROS_AT2 = p.km_ROS_AT2;
	k_AT1_AT2 = p.k_AT1_AT2;
	km_AT1_AT2 = p.km_AT1_AT2;
	k_IFNb_kill = p.k_IFNb_kill;
	k_kill = p.k_kill;
	km_kill = p.km_kill;
	km_ROS_AT2 = p.km_ROS_AT2;
	b_AT2 = p.b_AT2;
	b_I = p.b_I;
	mu_AT1 = p.mu_AT1;
	k_ROS_AT1 = p.k_ROS_AT1;
	km_ROS_AT1 = p.km_ROS_AT1;
	b_AT1 = p.b_AT1;
	b_dAT1 = p.b_dAT1;
	k_damage_TNFa = p.k_damage_TNFa;
	km_damage_TNFa = p.km_damage_TNFa;
	k_damage_IL6 = p.k_damage_IL6;
	km_damage_IL6 = p.km_damage_IL6;
	k_damage_IL1b = p.k_damage_IL1b;
	km_damage_IL1b = p.km_damage_IL1b;
	k_damage_IFNg = p.k_damage_IFNg;
	km_damage_IFNg = p.km_damage_IFNg;
	k_damage_cyt = p.k_damage_cyt;
	a_DC = p.a_DC;
	kbasal_DC = p.kbasal_DC;
	b_DC = kdc*p.b_DC;
	k_DC_TNFa = p.k_DC_TNFa;
	km_DC_TNFa = p.km_DC_TNFa;
	k_DC_IFNg = p.k_DC_IFNg;
	km_DC_IFNg = p.km_DC_IFNg;
	k_DC_IL6 = p.k_DC_IL6;
	km_DC_IL6 = p.km_DC_IL6;
	km_DC_IL10 = p.km_DC_IL10;
	a_M1 = p.a_M1;
	kbasal_M1 = p.kbasal_M1;
	k_v = p.k_v;
	k_I = p.k_I;
	km_v = p.km_v;
	km_I = p.km_I;
	k_M1_IL6 = p.k_M1_IL6;
	km_M1_IL6 = p.km_M1_IL6;
	b_M1 = kdc*p.b_M1;
	k_M1_TNFa = p.k_M1_TNFa;
	km_M1_TNFa = p.km_M1_TNFa;
	k_M1_GMCSF = p.k_M1_GMCSF;
	km_M1_GMCSF = p.km_M1_GMCSF;
	k_M1_IFNg = p.k_M1_IFNg;
	km_M1_IFNg = p.km_M1_IFNg;
	km_M1_IL10 = p.km_M1_IL10;
	a_N = p.a_N;
	k_N_IFNg = p.k_N_IFNg;
	km_N_IFNg = p.km_N_IFNg;
	k_N_TNFa = p.k_N_TNFa;
	km_N_TNFa = p.km_N_TNFa;
	k_N_GMCSF = p.k_N_GMCSF;
	km_N_GMCSF = p.km_N_GMCSF;
	k_N_IL17c = p.k_N_IL17c;
	km_N_IL17c = p.km_N_IL17c;
	b_N = kdc*p.b_N;
	a_Th1 = p.a_Th1;
	b_Th1 = kdc*p.b_Th1;
	k_Th1_IL2 = p.k_Th1_IL2;
	K_Th1_IL12 = p.K_Th1_IL12;
	k_Th1_IL12IL2 = p.k_Th1_IL12IL2;
	K_Th1_IL12IL2 = p.K_Th1_IL12IL2;
	K_Th1_IL10 = p.K_Th1_IL10;
	K_Th1_TGFb = p.K_Th1_TGFb;
	k_Th1_IFNg = p.k_Th1_IFNg;
	K_IFNg_Th1 = p.K_IFNg_Th1;
	K_Th1_IL6 = p.K_Th1_IL6;
	k_Th1_Th17 = p.k_Th1_Th17;
	K_Th1_Th17 = p.K_Th1_Th17;
	k_Th1_Treg = p.k_Th1_Treg;
	K_Th1_Treg = p.K_Th1_Treg;
	a_Th17 = p.a_Th17;
	b_Th17 = kdc*p.b_Th17;
	k_Th17_TGFb = p.k_Th17_TGFb;
	K_Th17_TGFb = p.K_Th17_TGFb;
	K_Th17_IL2 = p.K_Th17_IL2;
	K_Th17_IFNg = p.K_Th17_IFNg;
	K_Th17_IL10 = p.K_Th17_IL10;
	k_Th17_IL6 = p.k_Th17_IL6;
	km_Th17_IL6 = p.km_Th17_IL6;
	k_Th17_IL1b = p.k_Th17_IL1b;
	km_Th17_IL1b = p.km_Th17_IL1b;
	a_CTL = p.a_CTL;
	b_CTL = kdc*p.b_CTL;
	k_CTL_IL2 = p.k_CTL_IL2;
	K_CTL_IL12 = p.K_CTL_IL12;
	k_CTL_IL12IL2 = p.k_CTL_IL12IL2;
	K_CTL_IL12IL2 = p.K_CTL_IL12IL2;
	K_CTL_IL10 = p.K_CTL_IL10;
	K_CTL_TGFb = p.K_CTL_TGFb;
	K_CTL_IL6	 = p.K_CTL_IL6	;
	k_CTL_IFNg = p.k_CTL_IFNg;
	K_CTL_IFNg = p.K_CTL_IFNg;
	kmax_MHC1 = p.kmax_MHC1;
	km_MHC1_IFNb = p.km_MHC1_IFNb;
	a_Treg = p.a_Treg;
	b_Treg = kdc*p.b_Treg;
	k_Treg_IL2 = p.k_Treg_IL2;
	K_Treg_IL2 = p.K_Treg_IL2;
	K_Treg_IL17 = p.K_Treg_IL17;
	K_Treg_IL6 = p.K_Treg_IL6;
	k_Treg_TGFb = p.k_Treg_TGFb;
	K_Treg_TGFb = p.K_Treg_TGFb;
	kbasal_SPD = p.kbasal_SPD;
	a_SPD = p.a_SPD;
	b_SPD = p.b_SPD;
	kbasal_FER = p.kbasal_FER;
	a_FER = p.a_FER;
	b_FER = p.b_FER;
	a_tnf = p.a_tnf;
	a_tnf_at1 = p.a_tnf_at1;
	a_tnf_i = p.a_tnf_i;
	a_tnf_at2 = p.a_tnf_at2;
	a_tnf_m1 = p.a_tnf_m1;
	a_tnf_th1 = p.a_tnf_th1;
	a_tnf_th17 = p.a_tnf_th17;
	b_tnf = p.b_tnf;
	a_il6 = p.a_il6;
	b_il6 = p.b_il6;
	a_il6_at1 = p.a_il6_at1;
	a_il6_i = p.a_il6_i;
	a_il6_at2 = p.a_il6_at2;
	a_il6_m1 = p.a_il6_m1;
	a_il6_th17 = p.a_il6_th17;
	a_il6_neu = p.a_il6_neu;
	a_ifng = p.a_ifng;
	b_ifng = p.b_ifng;
	a_ifng_dc = p.a_ifng_dc;
	a_ifng_th1 = p.a_ifng_th1;
	a_ifng_ctl = p.a_ifng_ctl;
	a_ifnb = p.a_ifnb;
	b_ifnb = p.b_ifnb;
	a_ifnb_at1 = p.a_ifnb_at1;
	a_ifnb_i = p.a_ifnb_i;
	a_ifnb_d = p.a_ifnb_d;
	a_ifnb_dc = p.a_ifnb_dc;
	a_il2_dc = p.a_il2_dc;
	a_il2_th1 = p.a_il2_th1;
	b_il2 = p.b_il2;
	a_il2 = p.a_il2;
	a_il12_m1 = p.a_il12_m1;
	a_il12_dc = p.a_il12_dc;
	b_il12 = p.b_il12;
	a_il12 = p.a_il12;
	a_il17_th17 = p.a_il17_th17;
	a_il17_ctl = p.a_il17_ctl;
	b_il17 = p.b_il17;
	a_il17 = p.a_il17;
	a_il10_treg = p.a_il10_treg;
	b_il10 = p.b_il10;
	a_il10 = p.a_il10;
	a_tgfb_th17 = p.a_tgfb_th17;
	a_tgfb_treg = p.a_tgfb_treg;
	b_tgfb = p.b_tgfb;
	a_tgfb = p.a_tgfb;
	a_gmcsf_m1 = p.a_gmcsf_m1;
	a_gmcsf_th1 = p.a_gmcsf_th1;
	a_gmcsf_th17 = p.a_gmcsf_th17;
	b_gmcsf = p.b_gmcsf;
	a_gmcsf = p.a_gmcsf;
	a_il1b = p.a_il1b;
	b_il1b = p.b_il1b;
	a_il1b_at1 = p.a_il1b_at1;
	a_il1b_i = p.a_il1b_i;
	a_il1b_at2 = p.a_il1b_at2;
	a_il1b_m1 = p.a_il1b_m1;
	a_il1b_dc = p.a_il1b_dc;
	kbasal_ROS = p.kbasal_ROS;
	b_ROS = p.b_ROS;
	ktr_TNFa = p.ktr_TNFa;
	ktr_IL6 = p.ktr_IL6;
	ktr_IL1b = p.ktr_IL1b;
	ktr_IFNb = p.ktr_IFNb;
	ktr_IFNg = p.ktr_IFNg;
	ktr_IL2 = p.ktr_IL2;
	ktr_IL12 = p.ktr_IL12;
	ktr_IL17 = p.ktr_IL17;
	ktr_IL10 = p.ktr_IL10;
	ktr_TGFb = p.ktr_TGFb;
	ktr_GMCSF = p.ktr_GMCSF;
	ktr_SPD = p.ktr_SPD;
	ktr_FER = p.ktr_FER;
	k_CTL_I_SPD = p.k_CTL_I_SPD;
	k_CTL_I_Fer = p.k_CTL_I_Fer;
% disp(t)

	%%%% LUNG TO CENTRAL TRANSPORT RATES
	tr_TNFa = ktr_TNFa*s.TNFa;
	tr_IL6 = ktr_IL6*s.IL6;
	tr_IL1b = ktr_IL1b*s.IL1b;
	tr_IFNb = ktr_IFNb*s.IFNb;
	tr_IFNg = ktr_IFNg*s.IFNg;
	tr_IL2 = ktr_IL2*s.IL2;
	tr_IL12 = ktr_IL12*s.IL12;
	tr_IL17 = ktr_IL17*s.IL17;
	tr_IL10 = ktr_IL10*s.IL10;
	tr_TGFb = ktr_TGFb*s.TGFb;
	tr_GMCSF = ktr_GMCSF*s.GMCSF;
	tr_SPD = ktr_SPD*s.SPD;
	tr_FER = ktr_FER*s.FER;


	tr_pDC = p.ktr_pDC*(s.pDC/vol_alv_ml -s.pDC_c);
	tr_M1 =  p.ktr_M1*(s.M1/vol_alv_ml -s.M1_c);
	tr_N = p.ktr_N*(s.N/vol_alv_ml -s.N_c);
	tr_Th1 = p.ktr_Th1*(s.Th1/vol_alv_ml -s.Th1_c);
	tr_Th17 = p.ktr_Th17*(s.Th17/vol_alv_ml -s.Th17_c);
	tr_CTL = p.ktr_CTL*(s.CTL/vol_alv_ml -s.CTL_c);
	tr_Treg = p.ktr_Treg*(s.Treg/vol_alv_ml -s.Treg_c);

    %% viral dynamics

	%%  Virus (# viral mRNA/mL)
	prod_virus_shedding = A_V * s.I;
	virus_endocytosis = k_int * s.AT2 * s.V * (1 - s.IFNb / (km_int_IFNb + s.IFNb));
	innate_clearance = k_V_innate * s.V * (s.M1 + s.pDC);
	deg_virus = b_V * s.V;

	% not in use for now
	ab_clearance = 0*p.k_Ab_V*s.Ab*s.V;
	


	%% healthy alveolar type 2 cells (AT2) (# cells)
	ICAT1 = IC.AT1;%
	ICAT2 = IC.AT2;%
    damage_cyt_AT = 1*k_damage_cyt*(k_damage_TNFa*(s.TNFa/(km_damage_TNFa+s.TNFa)) +  k_damage_IL6*(s.IL6/(km_damage_IL6+s.IL6)) ...
        + k_damage_IL1b*(s.IL1b/(km_damage_IL1b + s.IL1b)) + k_damage_IFNg*(s.IFNg/(km_damage_IFNg+s.IFNg)));
	damage_ROS_AT2 = k_ROS_AT2*s.N/(km_ROS_AT2+s.N);	
	growth_AT2 = mu_AT2*(1 + p.k_mu_AT2*(max((ICAT1+ICAT2)-(s.AT2+s.AT1),0)/(km_AT1_AT2*(ICAT1+ICAT2) + max((ICAT1+ICAT2)-(s.AT2+s.AT1),0))))*s.AT2; % induced differentiation and growth induced by AT1 decrease from baseline
	deg_AT2 = b_AT2*s.AT2;
	diff_AT2 = k_AT1_AT2 * s.AT2 * (1 + p.k_diff_AT1*max(ICAT1-s.AT1,0)/(km_AT1_AT2*ICAT1 + max(ICAT1-s.AT1,0)));

	%% Infected AT2 (# cells)
	kill_CTL_I = k_kill*(1+k_IFNb_kill*s.IFNb/(km_kill+s.IFNb))*s.I*s.CTL;
	deg_I = b_I*s.I;


	%% health alveolar type 1 cells (AT1) (# cells)
	growth_AT1 = mu_AT1*(ICAT1-s.AT1); % no growth, no division, source is from AT2 
	damage_ROS_AT1 = damage_ROS_AT2*s.AT1;
	deg_AT1 = b_AT1*s.AT1;


	%% damaged AT1 cells (# cells)
	deg_dAT1 = b_dAT1*s.dAT1;

	%% damaged AT2 cells (# cells)
	deg_dAT2 = b_dAT1*s.dAT2;


	% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
	% % % % % Immune cells. % % % % % 

	%% lung pulmonary dendritic cells (# cells)
    act_pDC_TNFa = k_DC_TNFa*(s.TNFa/(km_DC_TNFa+s.TNFa));
	act_pDC_IFNg = k_DC_IFNg*(s.IFNg/(km_DC_IFNg+s.IFNg));
	act_pDC_IL6 = k_DC_IL6*(s.IL6/(km_DC_IL6+s.IL6));
    act_pDC_GMCSF = k_M1_GMCSF*(s.GMCSF/(km_M1_GMCSF+s.GMCSF));
	inh_pDC_IL10 = (km_DC_IL10/(km_DC_IL10+s.IL10));


	%% lung M1 macrophages (# cells)
	act_M1_TNFa = k_M1_TNFa*(s.TNFa/(km_M1_TNFa+s.TNFa));
	act_M1_GMCSF = k_M1_GMCSF*(s.GMCSF/(km_M1_GMCSF+s.GMCSF));
	act_M1_IFNg = k_M1_IFNg*(s.IFNg/(km_M1_IFNg+s.IFNg));
	inh_M1_IL10 = (km_M1_IL10/(km_M1_IL10+s.IL10));




	%% lung Neutrophils (# cells)
    act_N_IFNg = k_N_IFNg*s.IFNg/(s.IFNg+km_N_IFNg); 
    act_N_TNFa= k_N_TNFa*s.TNFa/(s.TNFa+km_N_TNFa); 
    act_N_GMCSF = k_N_GMCSF*s.GMCSF/(s.GMCSF+km_N_GMCSF);
    rec_N_IL17c = k_N_IL17c*s.IL17_c/(s.IL17_c + km_N_IL17c);


	%% lung Th1 Cells (# cells)
	act_Th1_IL12 = k_Th1_IL2*(s.IL12/(K_Th1_IL12+s.IL12))*(1+k_Th1_IL12IL2*(s.IL2/(K_Th1_IL12IL2+s.IL2)));
	act_Th1_IFNg = k_Th1_IFNg*(s.IFNg/(K_IFNg_Th1 + s.IFNg)) * (K_Th1_IL6/(K_Th1_IL6+s.IL6));
	inh_Th1_IL10_TGFb = (K_Th1_IL10/(K_Th1_IL10+s.IL10)) * (K_Th1_TGFb/(K_Th1_TGFb+s.TGFb));
    diff_Th1_Th17 = k_Th1_Th17*s.Th17*(s.IL12/(K_Th1_Th17+s.IL12))*(K_Th1_TGFb/(K_Th1_TGFb+s.TGFb));
	diff_Th1_Treg = k_Th1_Treg*s.Treg*(s.IL12/(K_Th1_Treg+s.IL12));


	%% lung TH17 cells (# cells)
	inh_Th17 = (K_Th17_IL2/(K_Th17_IL2 + s.IL2)) * (K_Th17_IFNg/(K_Th17_IFNg + s.IFNg)) * (K_Th17_IL10/(K_Th17_IL10 + s.IL10));
    act_TH17_TGFb = k_Th17_TGFb*(s.TGFb/(K_Th17_TGFb + s.TGFb)) * inh_Th17;
	act_Th17_IL6 = k_Th17_IL6*(s.IL6/(km_Th17_IL6 + s.IL6)) * inh_Th17;
	act_Th17_IL1b = k_Th17_IL1b*(s.IL1b/(km_Th17_IL1b + s.IL1b)) * inh_Th17;

	%% lung CTL (# cells)
	act_CTL_IFNb = kmax_MHC1 * (s.IFNb / (km_MHC1_IFNb + s.IFNb));
	act_CTL_IL12 = k_CTL_IL2 * (s.IL12 / (K_CTL_IL12 + s.IL12)) * (1 + k_CTL_IL12IL2 * (s.IL12 / (K_CTL_IL12IL2 + s.IL12))) * (K_CTL_IL10 / (K_CTL_IL10 + s.IL10)) * (K_CTL_TGFb / (K_CTL_TGFb+s.TGFb));
	act_CTL_IFNg = k_CTL_IFNg * (s.IFNg / (K_CTL_IFNg + s.IFNg)) * (K_CTL_IL10 / (K_CTL_IL10 + s.IL10)) * (K_CTL_TGFb / (K_CTL_TGFb + s.TGFb)) * (K_CTL_IL6 / (K_CTL_IL6 + s.IL6));


	%% lung Treg (# cells)
	act_Treg_IL2 = k_Treg_IL2 * (s.IL2 / (K_Treg_IL2 + s.IL2)) * (K_Treg_IL17 / (K_Treg_IL17 + s.IL17)) * (K_Treg_IL6 / (K_Treg_IL6 + s.IL6));
	act_Treg_TGFb = k_Treg_TGFb * (s.TGFb / (K_Treg_TGFb + s.TGFb)) * (K_Treg_IL17 / (K_Treg_IL17 + s.IL17)) * (K_Treg_IL6 / (K_Treg_IL6 + s.IL6));


	%% lung Ferritin (pmol)

	% Liver CRP Production (pmol)
    Liver_CRP = p.k_livercrp * p.Liver * ((s.IL6_c * p.Liver));    
    prod_CRP_liver = p.kCRPSecretion * Liver_CRP;
    tr_CRP = p.kCRP_BloodtoLiver * s.Blood_CRP - p.kCRP_LivertoBlood * Liver_CRP / p.Liver;
    prod_CRP_blood = p.kbasal_CRP;


	%% Cytokine Dynamics

	%% TNF-A (pmol)
	prod_tnf_dat1 = a_tnf_at1 * s.dAT1;
	prod_tnf_i = a_tnf_i * s.I;
	prod_tnf_dat2 = a_tnf_at2 * s.dAT2;
	prod_tnf_m1 = a_tnf_m1 * s.M1;
	prod_tnf_th1 = a_tnf_th1 * s.Th1;
	prod_tnf_th17 = a_tnf_th17 * s.Th17;
	deg_tnf = b_tnf * s.TNFa;
   
	%% IL-6 (pmol)
	prod_il6_dat1 = a_il6_at1 * s.dAT1;
	prod_il6_i =a_il6_i * s.I;
	prod_il6_dat2 = a_il6_at2 * s.dAT2;
	prod_il6_m1 = a_il6_m1 * s.M1;
	prod_il6_th17 = a_il6_th17 * s.Th17;
	prod_il6_neu =a_il6_neu * s.N;
	deg_il6 = b_il6 * s.IL6;

	%% IL1-B (pmol)
	prod_il1b_dat1 = a_il1b_at2 * s.dAT1;
	prod_il1b_i = a_il1b_i * s.I;
	prod_il1b_dat2 = a_il1b_at2 * s.dAT2;
	prod_il1b_m1 = a_il1b_m1 * s.M1;
    a_il1b_dc = a_il1b_m1;
	prod_il1b_dc = a_il1b_dc * s.pDC;
	deg_il1b = b_il1b * s.IL1b;

	%% IFN-G (pmol)
	prod_ifng_dc = a_ifng_dc * s.pDC;
	prod_ifng_th1 = a_ifng_th1 * s.Th1;
	prod_ifng_ctl = a_ifng_ctl * s.CTL;
	del_ifng = b_ifng * s.IFNg;


	%% IFN-B (pmol)
	prod_ifnb_i = a_ifnb_i * s.I;
	prod_ifnb_dc = a_ifnb_dc * s.pDC;
	del_ifnb = b_ifnb * s.IFNb;


	%% IL-2 (pmol)
	prod_il2_dc = a_il2_dc * s.pDC;
	prod_il2_th1 = a_il2_th1 * s.Th1;
	deg_il2 = b_il2 * s.IL2;


	%% IL-12 (pmol)
	prod_il12_m1 = a_il12_m1 * s.M1;
	prod_il12_dc = a_il12_dc * s.pDC;
	deg_il12 = b_il12 * s.IL12;


	%% IL-17 (pmol)
	prod_il17_th17 = a_il17_th17 * s.Th17;
	prod_il17_ctl = a_il17_ctl * s.CTL;
	deg_il17 = b_il17 * s.IL17;


	%% IL-10 (pmol)
	prod_il10_treg = a_il10_treg * s.Treg;
	deg_il10 = b_il10 * s.IL10;


	%% TGF-B (pmol)
	prod_tgfb_th17 = a_tgfb_th17 * s.Th17;
	prod_tgfb_treg = a_tgfb_treg * s.Treg;
	deg_tgfb = b_tgfb * s.TGFb;


	%% GM-CSF (pmol)
	prod_gmcsf_m1 = a_gmcsf_m1 * s.M1;
	prod_gmcsf_th1 = a_gmcsf_th1 * s.Th1;
	prod_gmcsf_th17 = a_gmcsf_th17 * s.Th17;
	deg_gmcsf = b_gmcsf * s.GMCSF;


	%% Assign Output

    rates = transpose([
    tr_TNFa  
	tr_IL6
	tr_IL1b 
	tr_IFNb 
	tr_IFNg 
	tr_IL2
	tr_IL12 
	tr_IL17 
	tr_IL10 
	tr_TGFb 
	tr_GMCSF
	tr_SPD
	tr_FER
	tr_pDC
	tr_M1 
	tr_N
	tr_Th1
	tr_Th17
	tr_CTL
	tr_Treg
	prod_virus_shedding
	virus_endocytosis
	innate_clearance
	deg_virus
	ab_clearance
    damage_cyt_AT
	damage_ROS_AT2
	growth_AT2
	deg_AT2
	diff_AT2
	kill_CTL_I
	deg_I
	growth_AT1
	damage_ROS_AT1
	deg_AT1
	deg_dAT1
	deg_dAT2
    act_pDC_TNFa
	act_pDC_IFNg
	act_pDC_IL6
    act_pDC_GMCSF
	inh_pDC_IL10
	act_M1_TNFa
	act_M1_GMCSF
	act_M1_IFNg 
	inh_M1_IL10 
    act_N_IFNg
    act_N_TNFa
    act_N_GMCSF
    rec_N_IL17c
	act_Th1_IL12
	act_Th1_IFNg
	inh_Th1_IL10_TGFb
    diff_Th1_Th17
	diff_Th1_Treg
	inh_Th17
    act_TH17_TGFb
	act_Th17_IL6
	act_Th17_IL1b 
	act_CTL_IFNb
	act_CTL_IL12
	act_CTL_IFNg
	act_Treg_IL2
	act_Treg_TGFb 
    Liver_CRP
    prod_CRP_liver
    tr_CRP
    prod_CRP_blood
	prod_tnf_dat1
	prod_tnf_i
	prod_tnf_dat2
	prod_tnf_m1
	prod_tnf_th1
	prod_tnf_th17
	deg_tnf
	prod_il6_dat1
	prod_il6_i
	prod_il6_dat2
	prod_il6_m1
	prod_il6_th17
	prod_il6_neu
	deg_il6
	prod_il1b_dat1
	prod_il1b_i
	prod_il1b_dat2
	prod_il1b_m1
    a_il1b_dc
	prod_il1b_dc
	deg_il1b
	prod_ifng_dc
	prod_ifng_th1
	prod_ifng_ctl
	del_ifng
	prod_ifnb_i
	prod_ifnb_dc
	del_ifnb 
	prod_il2_dc
	prod_il2_th1
	deg_il2
	prod_il12_m1
	prod_il12_dc
	deg_il12 
	prod_il17_th17
	prod_il17_ctl
	deg_il17 
	prod_il10_treg
	deg_il10 
	prod_tgfb_th17
	prod_tgfb_treg
	deg_tgfb 
	prod_gmcsf_m1
	prod_gmcsf_th1
	prod_gmcsf_th17
	deg_gmcsf]);
end
