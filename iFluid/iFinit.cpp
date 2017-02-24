#include <iFluid.h>
#include <ifluid_basic.h>

static double kh_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static double zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static void initRayleiTaylorIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initKHIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initRSRVIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initRSSYIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initTCIntfc(Front*,LEVEL_FUNC_PACK*,char*);
static void initCirclePlaneIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
static void initChannelFlow(Front*,LEVEL_FUNC_PACK*,char*);

extern void setInitialIntfc(
	Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	iFparams->m_comp1 = LIQUID_COMP1;
	iFparams->m_comp2 = LIQUID_COMP2;
        level_func_pack->is_RS_RV = NO;
        level_func_pack->is_RS_SY = NO;
        switch (prob_type)
        {
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            initCirclePlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case FLUID_SOLID_CIRCLE:
            iFparams->m_comp1 = SOLID_COMP;
            initCirclePlaneIntfc(front,level_func_pack,inname,prob_type);
            break;
        case TWO_FLUID_RT:
            initRayleiTaylorIntfc(front,level_func_pack,inname);
            break;
        case TWO_FLUID_KH:
            initKHIntfc(front,level_func_pack,inname);
            break;
        case TWO_FLUID_RS_RV:
            level_func_pack->is_RS_RV = YES;
            initRSRVIntfc(front,level_func_pack,inname);
            break;
        case TWO_FLUID_RS_SY:
            level_func_pack->is_RS_SY = YES;
            initRSSYIntfc(front,level_func_pack,inname);
            // TODO && FIXME: computeSourceTerm
            iFparams->ifluid_type = TWO_FLUID_RS_SY;
            break;
        case CHANNEL_FLOW:
            initChannelFlow(front,level_func_pack,inname);
            break;
        case TWO_FLUID_TC:
            initTCIntfc(front,level_func_pack,inname);
            break;
	default:
	    (void) printf("In setInitialIntfc unknown type: %d\n",prob_type);
        }
}       /* end setInitialIntfc */

static void initTCIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	static FOURIER_POLY level_func_params;
	static double L[MAXD],U[MAXD];
	FILE *infile = fopen(inname,"r");
	int i,j,dim,num_modes;
	char mesg[100],s[100];

	dim = level_func_params.dim = front->rect_grid->dim;
	level_func_params.L = L;
	level_func_params.U = U;
	for (i = 0; i < dim; ++i)
	{
	    level_func_params.L[i] = front->rect_grid->L[i];
	    level_func_params.U[i] = front->rect_grid->U[i];
	}

	level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params.z0);
	(void) printf("%f\n",level_func_params.z0);
        CursorAfterString(infile,"Enter the minimum mode numbers:");
        fscanf(infile,"%d",&level_func_params.min_n);
        (void) printf("%d\n",level_func_params.min_n);
        CursorAfterString(infile,"Enter the maximum mode numbers:");
        fscanf(infile,"%d",&level_func_params.max_n);
        (void) printf("%d\n",level_func_params.max_n);
        CursorAfterString(infile,"Enter the amplitude standard deviation:");
        fscanf(infile,"%lf",&level_func_params.A_sd);
        (void) printf("%f\n",level_func_params.A_sd);
        CursorAfterString(infile,"Enter the average phase:");
        fscanf(infile,"%lf",&level_func_params.av_phase);
        (void) printf("%f\n",level_func_params.av_phase);
        CursorAfterString(infile,"Enter the phase standard deviation:");
        fscanf(infile,"%lf",&level_func_params.P_sd);
        (void) printf("%f\n",level_func_params.P_sd);
	CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
	(void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
        sprintf(mesg,"Type y to use couette initial velocity:");
        iFparams->use_couette_init_vel = NO;
        if (CursorAfterStringOpt(infile,mesg))
        {
            fscanf(infile,"%s",s);
            if (s[0] == 'y' || s[0] == 'Y')
                iFparams->use_couette_init_vel = YES;
        }
        CursorAfterString(infile,"Enter constant initial velocoty:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->init_vel[i]);
            (void) printf("%f ",iFparams->init_vel[i]);
        }
        CursorAfterString(infile,"Enter inner velocoty:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->bvel[0][i]);
            (void) printf("%f ",iFparams->bvel[0][i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter outer velocoty:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->bvel[1][i]);
            (void) printf("%f ",iFparams->bvel[1][i]);
        }
        (void) printf("\n");
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	    (void) printf("%f ",iFparams->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);

	level_func_pack->func_params = (POINTER)&level_func_params;
        /* After this assignment, using level_func_pack->func &
         * level_func_pack->func_params to make level surface */

	if (front->coordinate == 'R' || front->coordinate == 'r')
	    level_func_pack->func = level_wave_func;
	else if(front->coordinate == 'C' || front->coordinate == 'c')
	    level_func_pack->func = level_wave_func_cylindrical;
	fclose(infile);
}	/* end initTCIntfc */


static void initRSSYIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        static FOURIER_POLY level_func_params;
        static double L[MAXD],U[MAXD];
        FILE *infile = fopen(inname,"r");
        int i,j,dim,num_modes;
        char mesg[100],s[Gets_BUF_SIZE],sbdry[200];
        double surfacetension, wv_len,grav,rhop1,rhop2;//for the calculation of wavelength
        int perturbed;// switch for perturbation theory, 1 is enabled; 0 is disabled. There is no default value. Need to fill in one.
        //variable naming rule and reference here:
        //surfacetension is surface tension
        //wv_len is most unstable wavelength from Smeeton Youngs' paper
        //grav is experimental gravity, thanks to very specific experiment, it's a scalar(only exists on Z direction)
        //rhop1, rhop2 are densities of 1 and 2

        if (debugging("trace"))
            (void) printf("Enter initRSSYIntfc()!\n");

        dim = level_func_params.dim = front->rect_grid->dim;
        level_func_params.L = L;
        level_func_params.U = U;
        for (i = 0; i < dim; ++i)
        {
            level_func_params.L[i] = front->rect_grid->L[i];
            level_func_params.U[i] = front->rect_grid->U[i];
        }

        level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
        // TODO && FIXME: REFLECTION BOUNDARY CONDITION
        level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
        CursorAfterString(infile,"Enter mean position of fluid interface:");
        fscanf(infile,"%lf",&level_func_params.z0);
        (void) printf("%12.8g\n",level_func_params.z0);
        CursorAfterString(infile,"Enter Perturbation Switch(1 or 0):");
        fscanf(infile,"%d",&perturbed);
        if (perturbed == 1)
            level_func_params.perturbed = YES;
        else
            level_func_params.perturbed = NO;
        (void) printf("%d\n",level_func_params.perturbed);
        CursorAfterString(infile,"Enter the amplitude standard deviation:");
        fscanf(infile,"%lf",&level_func_params.A_sd);
        (void) printf("%12.8g\n",level_func_params.A_sd);
        CursorAfterString(infile,"Enter the average phase:");
        fscanf(infile,"%lf",&level_func_params.av_phase);
        (void) printf("%12.8g\n",level_func_params.av_phase);
        CursorAfterString(infile,"Enter the phase standard deviation:");
        fscanf(infile,"%lf",&level_func_params.P_sd);
        (void) printf("%12.8g\n",level_func_params.P_sd);
        CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
        (void) printf("%12.8g %12.8g\n",iFparams->rho1,iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%12.8g %12.8g\n",iFparams->rho2,iFparams->mu2);
        // for vd
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->c1,&iFparams->Dcoef1);
        (void) printf("%12.8g %12.8g\n",iFparams->c1,iFparams->Dcoef1);
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->c2,&iFparams->Dcoef2);
        (void) printf("%12.8g %12.8g\n",iFparams->c2,iFparams->Dcoef2);

        CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->gravity[i]);
            (void) printf("%12.8g ",iFparams->gravity[i]);
        }
        rhop1 = iFparams->rho1;
        rhop2 = iFparams->rho2;
        grav = iFparams->gravity[2];
        //printf("Density 1 is %f; Density 2 is %f; Gravity is %f in function %s\n",rhop1, rhop2, grav,__func__);
        (void) printf("\n");
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%12.8g\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%12.8g\n",iFparams->smoothing_radius);
        surfacetension = iFparams->surf_tension;

        // Add Contact Angle for Smeeton Youngs Experiment 105
        // level_func_params->contact_angle == 0 => NO Mensicus
        CursorAfterString(infile,"Enter the degree of Contact Angle:");
        fscanf(infile,"%lf",&level_func_params.contact_angle);
        (void) printf("%12.8g\n", level_func_params.contact_angle);

        CursorAfterString(infile,"Enter the starting point of Meniscus:");
        fscanf(infile,"%lf",&level_func_params.Meniscus);
        (void) printf("%12.8g\n",level_func_params.Meniscus);
        //TODO && FIXME: Need to check contact angle and meniscus : if contact angle is 0, then no meniscus.
        if (level_func_params.contact_angle == 90.0)
            level_func_params.Meniscus = 0.0;
        front->contactangle = level_func_params.contact_angle;
        printf("in function %s Contact Angle is %24.24f\n", __func__, front->contactangle);


        // TODO && FIXME: This was missing.
        if (CursorAfterStringOpt(infile,"Enter volume fraction threshold:"))
        {
            fscanf(infile,"%lf",&iFparams->vol_frac_threshold);
            (void) printf("%lf\n",iFparams->vol_frac_threshold);
        }
        // TODO && FIXME: This was missing.
        for (i = 0; i < dim-1; ++i)
        {
            sprintf(mesg,"Enter the boundary type of perturbation in direction %d:",i);
            CursorAfterString(infile,mesg);
            fscanf(infile,"%s",sbdry);
            (void) printf("%s\n",sbdry);
            switch (sbdry[0])
            {
            case 'P':
            case 'p':
                level_func_params.pert_bdry_type[i] = PERIODIC;
                break;
            case 'S':
            case 's':
                level_func_params.pert_bdry_type[i] = SYMMETRIC;
                break;
            case 'U':
            case 'u':
            default:
                level_func_params.pert_bdry_type[i] = UNMODIFIED;
                break;
            }
        }
        // Introduce critical wavelength because of the existence of surface tension.
        if (surfacetension > 0.0)
        {
             wv_len = 2*PI*sqrt(3*surfacetension/(grav*fabs(rhop1-rhop2)));
            level_func_params.wv_len = wv_len;
            //printf("most unstable wavelength is %f in func %s\n", wv_len,__func__);
        }
        else
        {
            printf("ERROR! NO surface tension was considered! More Implementation/Consideration Needed.\n");
             clean_up(ERROR);
        }

        level_func_pack->func_params = (POINTER)&level_func_params;
        /* After this assignment, using level_func_pack->func &
         * level_func_pack->func_params to make level surface */

        level_func_pack->func = level_wave_func_Meniscus;
        fclose(infile);
}       /* end initRSSYIntfc */


static void initRayleiTaylorIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	static FOURIER_POLY level_func_params;
	static double L[MAXD],U[MAXD];
	FILE *infile = fopen(inname,"r");
	int i,j,dim,num_modes;
	char mesg[100];

	dim = level_func_params.dim = front->rect_grid->dim;
	level_func_params.L = L;
	level_func_params.U = U;
	for (i = 0; i < dim; ++i)
	{
	    level_func_params.L[i] = front->rect_grid->L[i];
	    level_func_params.U[i] = front->rect_grid->U[i];
	}

	level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params.z0);
	(void) printf("%f\n",level_func_params.z0);

	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	(void) printf("%d\n",num_modes);
	level_func_params.num_modes = num_modes;
	FT_MatrixMemoryAlloc((POINTER*)&level_func_params.nu,num_modes,
			dim-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.phase,num_modes,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.A,num_modes,sizeof(double));
	for (i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (j = 0; j < dim; ++j)
	    {
	    	fscanf(infile,"%lf",&level_func_params.nu[i][j]);
		(void) printf("%f ",level_func_params.nu[i][j]);
	    }
	    (void) printf("\n");
	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.A[i]);
	    (void) printf("%f\n",level_func_params.A[i]);
	    sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.phase[i]);
	    (void) printf("%f\n",level_func_params.phase[i]);
	}
	CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
	(void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

        //for vd
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->c1,&iFparams->Dcoef1);
        (void) printf("%f %f\n",iFparams->c1,iFparams->Dcoef1);
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->c2,&iFparams->Dcoef2);
        (void) printf("%f %f\n",iFparams->c2,iFparams->Dcoef2);

	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	    (void) printf("%f ",iFparams->gravity[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);

        iFparams->width_idl = 0;
        if (CursorAfterStringOpt(infile,"Enter width of initial diffusion layer:"))
        {
            fscanf(infile,"%lf",&iFparams->width_idl);
            (void) printf("%lf\n",iFparams->width_idl);
        }

        iFparams->vol_frac_threshold = 0.05;
        if (CursorAfterStringOpt(infile,"Enter volume fraction threshold:"))
        {
            fscanf(infile,"%lf",&iFparams->vol_frac_threshold);
            (void) printf("%lf\n",iFparams->vol_frac_threshold);
        }

	level_func_pack->func_params = (POINTER)&level_func_params;
	level_func_pack->func = level_wave_func;
	fclose(infile);
}	/* end initRayleiTaylor */

static void initKHIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	static FOURIER_POLY level_func_params;
	static double L[MAXD],U[MAXD];
	FILE *infile = fopen(inname,"r");
	int i,j,dim,num_modes;
	char mesg[100];

	dim = level_func_params.dim = front->rect_grid->dim;
	level_func_params.L = L;
	level_func_params.U = U;
	for (i = 0; i < dim; ++i)
	{
	    level_func_params.L[i] = front->rect_grid->L[i];
	    level_func_params.U[i] = front->rect_grid->U[i];
	}

	level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	CursorAfterString(infile,"Enter mean position of fluid interface:");
	fscanf(infile,"%lf",&level_func_params.z0);
	CursorAfterString(infile,"Enter number of sine modes:");
	fscanf(infile,"%d",&num_modes);
	level_func_params.num_modes = num_modes;
	FT_MatrixMemoryAlloc((POINTER*)&level_func_params.nu,num_modes,
			dim-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.phase,num_modes,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&level_func_params.A,num_modes,sizeof(double));
	for (i = 0; i < num_modes; ++i)
	{
	    sprintf(mesg,"Enter frequency of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    for (j = 0; j < dim; ++j)
	    	fscanf(infile,"%lf",&level_func_params.nu[i][j]);
	    sprintf(mesg,"Enter amplitude of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.A[i]);
	    sprintf(mesg,"Enter phase of mode %d:",i+1);
	    CursorAfterString(infile,mesg);
	    fscanf(infile,"%lf",&level_func_params.phase[i]);
	}
	CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        CursorAfterString(infile,"Enter bottom horizontal velocity:");
	for (j = 0; j < dim-1; ++j)
            fscanf(infile,"%lf",&iFparams->U1[j]);
        CursorAfterString(infile,"Enter top horizontal velocity:");
	for (j = 0; j < dim-1; ++j)
            fscanf(infile,"%lf",&iFparams->U2[j]);

	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);

	level_func_pack->func_params = (POINTER)&level_func_params;
	level_func_pack->func = level_wave_func;
	fclose(infile);
}	/* end initKH */


static void initRSRVIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        static FOURIER_POLY level_func_params;
        static double L[MAXD],U[MAXD];
        FILE *infile = fopen(inname,"r");
        int i,j,dim,num_modes;
        char mesg[200],s[Gets_BUF_SIZE],sbdry[200];

        if (debugging("trace"))
            (void) printf("Enter initRSRVIntfc()!\n");

        dim = level_func_params.dim = front->rect_grid->dim;
        level_func_params.L = L;
        level_func_params.U = U;
        for (i = 0; i < dim; ++i)
        {
            level_func_params.L[i] = front->rect_grid->L[i];
            level_func_params.U[i] = front->rect_grid->U[i];
        }

        level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
        level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
        CursorAfterString(infile,"Enter mean position of fluid interface:");
        fscanf(infile,"%lf",&level_func_params.z0);
        (void) printf("%12.8g\n",level_func_params.z0);
        CursorAfterString(infile,"Enter the minimum mode numbers:");
        fscanf(infile,"%d",&level_func_params.min_n);
        (void) printf("%d\n",level_func_params.min_n);        CursorAfterString(infile,"Enter the maximum mode numbers:");        fscanf(infile,"%d",&level_func_params.max_n);
        (void) printf("%d\n",level_func_params.max_n);        CursorAfterString(infile,"Enter the amplitude standard deviation:");
        fscanf(infile,"%lf",&level_func_params.A_sd);        (void) printf("%12.8g\n",level_func_params.A_sd);        CursorAfterString(infile,"Enter the average phase:");        fscanf(infile,"%lf",&level_func_params.av_phase);
        (void) printf("%12.8g\n",level_func_params.av_phase);        CursorAfterString(infile,"Enter the phase standard deviation:");
        fscanf(infile,"%lf",&level_func_params.P_sd);        (void) printf("%12.8g\n",level_func_params.P_sd);

        for (i = 0; i < dim-1; ++i)
        {
            sprintf(mesg,"Enter the boundary type of perturbation in direction %d:",i);
            CursorAfterString(infile,mesg);
            fscanf(infile,"%s",sbdry);
            (void) printf("%s\n",sbdry);
            switch (sbdry[0])
            {
            case 'P':
            case 'p':
                level_func_params.pert_bdry_type[i] = PERIODIC;
                break;
            case 'S':
            case 's':
                level_func_params.pert_bdry_type[i] = SYMMETRIC;
                break;
            case 'U':
            case 'u':
            default:
                level_func_params.pert_bdry_type[i] = UNMODIFIED;
                break;
            }
        }

        CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
        (void) printf("%12.8g %12.8g\n",iFparams->rho1,iFparams->mu1);
        CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%12.8g %12.8g\n",iFparams->rho2,iFparams->mu2);

        //for vd
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 1:");
        fscanf(infile,"%lf %lf",&iFparams->c1,&iFparams->Dcoef1);
        (void) printf("%12.8g %12.8g\n",iFparams->c1,iFparams->Dcoef1);
        CursorAfterString(infile,"Enter concentration and binary diffusivity of fluid 2:");
        fscanf(infile,"%lf %lf",&iFparams->c2,&iFparams->Dcoef2);
        (void) printf("%12.8g %12.8g\n",iFparams->c2,iFparams->Dcoef2);

        CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&iFparams->gravity[i]);
            (void) printf("%12.8g ",iFparams->gravity[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%12.8g\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%12.8g\n",iFparams->smoothing_radius);

        iFparams->width_idl = 0;
        if (CursorAfterStringOpt(infile,"Enter width of initial diffusion layer:"))
        {
            fscanf(infile,"%lf",&iFparams->width_idl);
            (void) printf("%lf\n",iFparams->width_idl);
        }

        iFparams->vol_frac_threshold = 0.05;
        if (CursorAfterStringOpt(infile,"Enter volume fraction threshold:"))
        {
            fscanf(infile,"%lf",&iFparams->vol_frac_threshold);
            (void) printf("%lf\n",iFparams->vol_frac_threshold);
        }

        level_func_pack->func_params = (POINTER)&level_func_params;
        //level_func_pack->func = random_pert_vd_func;
        level_func_pack->func = level_wave_func_cylindrical;
        fclose(infile);
} /* end initRSRVIntfc */


static void initCirclePlaneIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	IF_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	static CIRCLE_PARAMS *circle_params;
	int i,dim;

	iFparams = (IF_PARAMS*)front->extra1;
	FT_ScalarMemoryAlloc((POINTER*)&circle_params,sizeof(CIRCLE_PARAMS));
        circle_params->dim = dim = front->rect_grid->dim;
        circle_params->add_plan_surf = NO;
        CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&circle_params->cen[i]);
            (void) printf("%f ",circle_params->cen[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter radius of the circle:");
        fscanf(infile,"%lf",&circle_params->R);
        (void) printf("%f\n",circle_params->R);

        if (prob_type == BUBBLE_SURFACE)
        {
            CursorAfterString(infile,"Enter height of the surface:");
            fscanf(infile,"%lf",&circle_params->H);
            (void) printf("%f\n",circle_params->H);
            circle_params->add_plan_surf = YES;
        }
	level_func_pack->func_params = (POINTER)circle_params;

	switch (prob_type)
	{
	case TWO_FLUID_BUBBLE:
        case BUBBLE_SURFACE:
            level_func_pack->neg_component = LIQUID_COMP1;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	    CursorAfterString(infile,"Enter density and viscosity of fluid 1:");
            fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
            (void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
            CursorAfterString(infile,"Enter density and viscosity of fluid 2:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
        case FLUID_SOLID_CIRCLE:
	    iFparams->m_comp1 = SOLID_COMP;
            level_func_pack->neg_component = SOLID_COMP;
            level_func_pack->pos_component = LIQUID_COMP2;
            level_func_pack->func = level_circle_func;
            level_func_pack->wave_type = NEUMANN_BOUNDARY;
            CursorAfterString(infile,
			"Enter density and viscosity of the fluid:");
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
            break;
	}
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf\n",&iFparams->gravity[i]);
	CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
        (void) printf("%f\n",iFparams->surf_tension);
	CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
        (void) printf("%f\n",iFparams->smoothing_radius);
	if (prob_type == TWO_FLUID_BUBBLE &&
	    iFparams->surf_tension != 0.0)
	{
	    double mag_g = mag_vector(iFparams->gravity,dim);
	    double rho_b = iFparams->rho1;
	    double rho_o = iFparams->rho2;
	    double mu_b = iFparams->mu1;
	    double mu_o = iFparams->mu2;
	    double sigma = iFparams->surf_tension;
	    double de = 2.0*circle_params->R;
	    double V_morton = mag_g*sqr(mu_o)*sqr(mu_o)/rho_o/sqr(sigma)/sigma;
	    double V_eotvos = rho_o*mag_g*sqr(de)/sigma;
	    printf("Bubble simulation\n");
	    printf("Morton number   = %g\n",V_morton);
	    printf("Eotvos number   = %g\n",V_eotvos);
	    printf("Density ratio   = %g\n",rho_o/rho_b);
	    printf("Viscosity ratio = %g\n",mu_o/mu_b);
	}

	fclose(infile);
}	/* end initCirclePlaneIntfc */

extern void init_fluid_state_func(
	Incompress_Solver_Smooth_Basis *l_cartesian,
	IF_PROB_TYPE prob_type)
{
	switch (prob_type)
	{
	case TWO_FLUID_KH:
	    l_cartesian->getInitialState = kh_state;
	    break;
	default:
	    l_cartesian->getInitialState = zero_state;
	}
}	/* end init_fluid_state_func */

static double kh_state(
	COMPONENT comp,
	double *coords,
	L_STATE& state,
	int dim,
	IF_PARAMS *iFparams)
{
	int i;
	for (i = 0; i < dim-1; ++i)
	{
	    state.m_U[i] = (comp == LIQUID_COMP1) ? iFparams->U1[i] :
				iFparams->U2[i];
	}
	state.m_U[dim-1] = 0.0;
}	/* end kh_state */

static double zero_state(
	COMPONENT comp,
	double *coords,
	L_STATE& state,
	int dim,
	IF_PARAMS *iFparams)
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = 0.0;
}	/* end zero_state */


extern void read_iF_prob_type(
	char *inname,
	IF_PROB_TYPE *prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	*prob_type = ERROR_TYPE;
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	*prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
            {
                if (string[11] == 'T' || string[11] == 't')
	    	    *prob_type = TWO_FLUID_RT;
                else if ((string[11] == 'S' || string[11] == 's') && (string[13] == 'R' || string[13] == 'r'))
                    *prob_type = TWO_FLUID_RS_RV;
                else if ((string[11] == 'S' || string[11] == 's') && (string[13] == 'S' || string[13] == 's'))
                    *prob_type = TWO_FLUID_RS_SY;
            }
	    else if (string[10] == 'K' || string[10] == 'k')
	    	*prob_type = TWO_FLUID_KH;
            else if (string[10] == 'T' || string[10] == 't')
                *prob_type = TWO_FLUID_TC;
	}
	else if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[6] == 'S' || string[6] == 's')
                *prob_type = FLUID_SOLID_CIRCLE;
            else if (string[6] == 'R' || string[6] == 'r')
                *prob_type = FLUID_RIGID_BODY;
	}
	else if (string[0] == 'B' || string[0] == 'b')
	{
	    *prob_type = BUBBLE_SURFACE;
	}
        else if (string[0] == 'R' || string[0] == 'r')
	{
            if (string[6] == 'O' || string[6] == 'o')
                *prob_type = ROTOR_ONE_FLUID;
            else if (string[6] == 'T' || string[6] == 't')
                *prob_type = ROTOR_TWO_FLUID;
	}
	else if (string[0] == 'C' || string[0] == 'c')
	{
            *prob_type = CHANNEL_FLOW;
	}

	assert(*prob_type != ERROR_TYPE);
	fclose(infile);
}	/* end read_iFparams */

static void initChannelFlow(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	int i,dim;

	iFparams = (IF_PARAMS*)front->extra1;

        level_func_pack->neg_component = LIQUID_COMP1;
        level_func_pack->pos_component = LIQUID_COMP2;
        level_func_pack->func = NULL;
        CursorAfterString(infile,"Enter density and viscosity of the fluid:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
        (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf",&iFparams->gravity[i]);

	fclose(infile);
}	/* end initChannelFlow */

