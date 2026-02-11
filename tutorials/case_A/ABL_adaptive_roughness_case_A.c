#include "udf.h"
#include "mem.h"
#include "models.h"

// Definition of minimum and maximum functions
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define Max(a, b) ((a) > (b) ? (a) : (b))

// Universal constants
#define Ec 9.793
real kappa = KAPPA;
real Cmu;

/*******************************************************************************/
/* INPUT DATA                                                                  */
/*******************************************************************************/

// Parameters for the power-law velocity distribution
// delta (m) - Height of the Atmospheric Boundary Layer (ABL) for vertical blending
#define delta 0.48191953E+00

// alpha_b (-) - Blending function coefficient
#define alpha_b 0.75905552E+01

// z_alpha (-) - Reference height for power-law velocity profile
#define z_alpha 0.14932874E-04

// u_alpha (ms^(-1)) - Wind speed at z_alpha
#define u_alpha 0.44724820E+00

// alpha (-) - Power law exponent
#define alpha 0.27073645E+00

/*******************************************************************************/

// Parameters of the logarithmic-polynomial TKE distribution, see Equation (4)
// Constants A, B, C and D, specific to kappa = 0.4187, as they are presented in Table (3)
#define A_k  0.2339825
#define B_k -0.1110081
#define C_k  0.0116020
#define D_k  0.49123987

// z_0k (m) - Roughness height for inflow TKE distribution
#define z0k 0.50000000E-02

/*******************************************************************************/

// Parameters for the roughness calculations
// u_star_ABL (ms^(-1)) - Friction velocity for the near-wall logarithmic velocity distibution
// u_star_ABL = u_ref*kappa/log((z_ref + z0)/z0), where kappa = 0.4187, see Equation (2)
#define u_star_ABL 0.38389056E+00

// Set flow direction vector (x, y, z)
real flow_dir[3] = {1.0, 0.0, 0.0};

// Number of aerodynamically rough surfaces
int n_rough_walls = 1;

// Thread ID-s of aerodynamically rough surfaces as a 1D array
int rough_wall_ID[1] = {10};

// z_0 (m) - Aerodynamic roughness height values as a 1D array
real z0_array[1] = {0.0010604866};

// Ks+ (-) - Selected dimensionless sand-grain roughness height values as a 1D array
real ks_plus_array[1] = {1000.0};

// Cs_min (-) - Lower limit for the Cs roughness constant to avoid singularities
real Cs_min = 1e-3;

/*******************************************************************************/
/* FUNCTION PROTOTYPES 														   */
/*******************************************************************************/

// Blending function and its derivative
real f_blend(real z);
real df_blend_dz(real z);

// zmod and its derivative
real zmod(real z);
real dzmod_dz(real z);

// Velocity power law profile and its partial and total derivative
real vel_in(real z);
real partial_du_dz(real z);
real du_dz(real z);

// Turbulent kinetic energy profile
real tke_in(real z);

// Turbulent dissipation rate profile
real epsilon_in(real z);

// Specific dissipation rate profile
real omega_in(real z);

/*******************************************************************************/

DEFINE_EXECUTE_ON_LOADING(set_C_mu, libname)
{
	Cmu = GVAR_TURB(coeff, ke_Cmu);
}

/*******************************************************************************/

DEFINE_EXECUTE_ON_LOADING(name_udmi, libname)
{
	Set_User_Memory_Name(0,"Ks");
	Set_User_Memory_Name(1,"zp");
	Set_User_Memory_Name(2,"zp_plus");
	Set_User_Memory_Name(3,"u_star");
	Set_User_Memory_Name(4,"u_star_ABL_calc");
	Set_User_Memory_Name(4,"Cs");
}

/*******************************************************************************/
/* PATCH IN INLET PROFILES (K-EPSILON MODEL)                                   */
/*******************************************************************************/

DEFINE_ON_DEMAND(Patch_ABL_k_e)
{
  real x[ND_ND];
  real z, zm, umag;

  cell_t c;
  Thread *t;
  Domain *domain;

  Message0("\nINTERPOLATING PROFILES INTO THE CELL CENTERS ...\n");
  domain = Get_Domain(1);
  // Loop through the threads of the domain
  thread_loop_c(t,domain)
  {
	// Loop through the cells of the thread
    begin_c_loop(c,t)
    {
      // Calculate vertical coordinate
      C_CENTROID(x,c,t);
      z = x[ND_ND-1];
	  // Apply blending function to vertical coordinate
      zm = zmod(z);
      umag = vel_in(zm);
	  // Assign values u, v, w, k and epsilon to the cell from inflow distribution
      C_U(c,t) = flow_dir[0]*umag;
      C_V(c,t) = flow_dir[1]*umag;
      C_W(c,t) = flow_dir[2]*umag;
      C_K(c,t) = tke_in(zm);
      C_D(c,t) = epsilon_in(zm);
    }
    end_c_loop(c,t)
  }
  Message0("DONE\n");
}

/*******************************************************************************/
/* PATCH IN INLET PROFILES (K-OMEGA SST MODEL)                                 */
/*******************************************************************************/

DEFINE_ON_DEMAND(Patch_ABL_k_w)
{
  real x[ND_ND];
  real z, zm, umag;

  cell_t c;
  Thread *t;
  Domain *domain;

  Message0("\nINTERPOLATING PROFILES INTO THE CELL CENTERS ...\n");
  domain = Get_Domain(1);
  // Loop through the threads of the domain
  thread_loop_c(t,domain)
  {
	// Loop through the cells of the thread
    begin_c_loop(c,t)
    {
      // Calculate vertical coordinate
      C_CENTROID(x,c,t);
      z = x[ND_ND-1];
	  // Apply blending function to vertical coordinate
      zm = zmod(z);
      umag = vel_in(zm);
	  // Assign values u, v, w, k and omega to the cell from inflow distribution
      C_U(c,t) = flow_dir[0]*umag;
      C_V(c,t) = flow_dir[1]*umag;
      C_W(c,t) = flow_dir[2]*umag;
      C_K(c,t) = tke_in(zm);
      C_O(c,t) = omega_in(zm);
    }
    end_c_loop(c,t)
  }
  Message0("DONE\n");
}

/*******************************************************************************/
/* DEFINE INLET STREAMWISE VELOCITY PROFILE                                    */
/*******************************************************************************/

DEFINE_PROFILE(inlet_U,t,i)
{
  real x[ND_ND];
  real z, zm, umag;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  {
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
    umag = vel_in(zm);
	// Assign velocity value to the cell face
    F_PROFILE(f,t,i) = flow_dir[0]*umag;
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* DEFINE INLET LATERAL VELOCITY PROFILE                                       */
/*******************************************************************************/

DEFINE_PROFILE(inlet_V,t,i)
{
  real x[ND_ND];
  real z, zm, umag;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  {
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
    umag = vel_in(zm);
	// Assign velocity value to the cell face
    F_PROFILE(f,t,i) = flow_dir[1]*umag;
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* DEFINE INLET VERTICAL VELOCITY PROFILE                                      */
/*******************************************************************************/

DEFINE_PROFILE(inlet_W,t,i)
{
  real x[ND_ND];
  real z, zm, umag;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  {
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
    umag = vel_in(zm);
	// Assign velocity value to the cell face
    F_PROFILE(f,t,i) = flow_dir[2]*umag;
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* DEFINE INLET TKE PROFILE                                                    */
/*******************************************************************************/

DEFINE_PROFILE(inlet_k,t,i)
{
  real x[ND_ND];
  real z, zm;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  {
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
	// Assign TKE value to the cell face
    F_PROFILE(f,t,i) = tke_in(zm);
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* DEFINE INLET EPSILON PROFILE                                                */
/*******************************************************************************/

DEFINE_PROFILE(inlet_epsilon,t,i)
{
  real x[ND_ND];
  real z, zm;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  { 
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
	// Assign epsilon value to the cell face
    F_PROFILE(f,t,i) = epsilon_in(zm);
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* DEFINE INLET OMEGA PROFILE                                                  */
/*******************************************************************************/

DEFINE_PROFILE(inlet_omega,t,i)
{
  real x[ND_ND];
  real z, zm;
  face_t f;

  // Loop through the faces of the thread
  begin_f_loop(f,t)
  { 
    F_CENTROID(x,f,t);
    z = x[ND_ND-1];
    zm = zmod(z);
	// Assign omega value to the cell face
    F_PROFILE(f,t,i) = omega_in(zm);
  }
  end_f_loop(f,t)
}

/*******************************************************************************/
/* CAlCULATE Ks FROM THE FIXED Ks+                                             */
/*******************************************************************************/

DEFINE_PROFILE(wall_ks,t,i)
{
  real z0;
  real ks_plus;
  face_t f;
  cell_t c0;
  Thread *t0;
  real u_star;
  
  int zone_ID; 
  int index;

  zone_ID = THREAD_ID(t);

  // Loop through list of rough wall ID-s
  for(index = 0; index < n_rough_walls; index++)
  {
	if (zone_ID == rough_wall_ID[index])
	{
		// Assign z0 and Ks+ to the thread based on its ID
		z0 = z0_array[index];
		ks_plus = ks_plus_array[index];
	}
  }	
  
  // Loop through the faces of the thread
  begin_f_loop(f,t)
  { 
	c0 = F_C0(f,t);
	t0 = THREAD_T0(t);
	// Calculate u*, see Equation (10)
	u_star = pow(Cmu,0.25)*pow(C_K(c0,t0),0.5);
	// Assign calculated Ks value to the face
    F_PROFILE(f,t,i) = ks_plus*C_MU_L(c0,t0)/u_star;
	C_UDMI(c0,t0,0) = ks_plus*C_MU_L(c0,t0)/u_star;
  }
  end_f_loop(f,t)

}

/*******************************************************************************/
/* CALCULATE Cs BASED ON Ks+ ACCORDING TO METHOD 1                             */
/*******************************************************************************/

DEFINE_PROFILE(wall_cs_M1,t,i)
{
	face_t f;
	cell_t c0;
	Thread *t0;
	real zp;
	real zp_plus;
	real z0;
	real ks_plus;
	real Cs;
	real u_p[ND_ND]; 
	real u_p_mag;
	real u_star;
	real u_star_ABL_calc;
	real tau_w_per_rho;
	real phi_wall;
	real theta_wall;
	
	real xf[ND_ND];
	real dx[ND_ND]; 
	real dx_mag;
	real xc[ND_ND];
	
	int zone_ID;
	int index;
	
	zone_ID = THREAD_ID(t);
	
	// Loop through list of rough wall ID-s
	for(index = 0; index < n_rough_walls; index++)
	{
		if (zone_ID == rough_wall_ID[index])
		{
			// Assign z0 and Ks+ to the thread based on its ID
			z0 = z0_array[index];
			ks_plus = ks_plus_array[index];
		}
	}	

	// Loop through the faces of the thread
	begin_f_loop(f,t)
	{
		c0 = F_C0(f,t);
		t0 = THREAD_T0(t);
		
		F_CENTROID(xf, f, t);
		C_CENTROID(xc, c0, t0);
				
		dx[0] = xc[0] - xf[0];
		dx[1] = xc[1] - xf[1];
		dx[2] = xc[2] - xf[2];
		dx_mag = NV_MAG(dx);
		zp = dx_mag;
		// Calculate u* - see Equation (10)
		u_star = pow(Cmu,0.25)*pow(C_K(c0,t0),0.5);
		// Calculate zp+ - see Section 2.4
		zp_plus = zp*u_star/C_MU_L(c0,t0);
		
		u_p[0] = C_U(c0,t0);
		u_p[1] = C_V(c0,t0);
		u_p[2] = C_W(c0,t0);
		u_p_mag = NV_MAG(u_p);
		
		u_star_ABL_calc = u_p_mag*kappa/log((zp+z0)/z0);
		// Approximate tau_w/rho using Method 1 - see Equation (18)
		tau_w_per_rho = C_MU_L(c0,t0)*u_p_mag/zp+pow(Cmu,0.5)*C_K(c0,t0);
		// Calculate exponent phi using Method 1 - see Equation (19)
		phi_wall = u_star*u_star_ABL/tau_w_per_rho;
		// Calculate theta - see Equation (25)
		theta_wall = 1./Ec*pow((zp+z0)/z0,phi_wall);
		
		// Calculate Cs - see Equation (24)
		if (ks_plus > zp_plus)
		{
			Cs = Max((zp_plus+ks_plus/2-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		else
		{
			Cs = Max((zp_plus-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		
		// Store diagnostic variables in UDMI
		C_UDMI(c0,t0,1) = phi_wall;
		C_UDMI(c0,t0,2) = zp_plus;
		C_UDMI(c0,t0,3) = u_star;
		C_UDMI(c0,t0,4) = u_star_ABL_calc;
		C_UDMI(c0,t0,5) = Cs;
		
		// Return Cs as a profile
		F_PROFILE(f,t,i) = Cs;
	}
	end_f_loop(f,t)

}

/*******************************************************************************/
/* CALCULATE Cs BASED ON Ks+ ACCORDING TO METHOD 2                             */
/*******************************************************************************/

DEFINE_PROFILE(wall_cs_M2,t,i)
{
	face_t f;
	cell_t c0;
	Thread *t0;
	real zp;
	real zp_plus;
	real z0;
	real ks_plus;
	real Cs;
	real u_p[ND_ND]; 
	real u_p_mag;
	real u_star;
	real u_star_ABL_calc;
	real tau_w_per_rho;
	real phi_wall;
	real theta_wall;
	
	real xf[ND_ND];
	real dx[ND_ND]; 
	real dx_mag;
	real xc[ND_ND];
	
	int zone_ID;
	int index;
	
	zone_ID = THREAD_ID(t);
	
	// Loop through list of rough wall ID-s
	for(index = 0; index < n_rough_walls; index++)
	{
		if (zone_ID == rough_wall_ID[index])
		{
			// Assign z0 and Ks+ to the thread based on its ID
			z0 = z0_array[index];
			ks_plus = ks_plus_array[index];
		}
	}	

	// Loop through the faces of the thread
	begin_f_loop(f,t)
	{
		c0 = F_C0(f,t);
		t0 = THREAD_T0(t);
		
		F_CENTROID(xf, f, t);
		C_CENTROID(xc, c0, t0);
				
		dx[0] = xc[0] - xf[0];
		dx[1] = xc[1] - xf[1];
		dx[2] = xc[2] - xf[2];
		dx_mag = NV_MAG(dx);
		zp = dx_mag;
		// Calculate u* - see Equation (10)
		u_star = pow(Cmu,0.25)*pow(C_K(c0,t0),0.5);
		// Calculate zp+ - see Section 2.4
		zp_plus = zp*u_star/C_MU_L(c0,t0);
		
		u_p[0] = C_U(c0,t0);
		u_p[1] = C_V(c0,t0);
		u_p[2] = C_W(c0,t0);
		u_p_mag = NV_MAG(u_p);
		
		u_star_ABL_calc = u_p_mag*kappa/log((zp+z0)/z0);
		// Approximate tau_w/rho using Method 2 - see Equation (20)
		tau_w_per_rho = (C_MU_L(c0,t0)+C_MU_T(c0,t0))*u_p_mag/zp;
		// Calculate exponent phi using Method 2 - see Equation (21)
		phi_wall = u_star*u_star_ABL/tau_w_per_rho;
		// Calculate theta - see Equation (25)
		theta_wall = 1./Ec*pow((zp+z0)/z0,phi_wall);
		
		// Calculate Cs - see Equation (24)
		if (ks_plus > zp_plus)
		{
			Cs = Max((zp_plus+ks_plus/2-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		else
		{
			Cs = Max((zp_plus-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		
		// Store diagnostic variables in UDMI
		C_UDMI(c0,t0,1) = phi_wall;
		C_UDMI(c0,t0,2) = zp_plus;
		C_UDMI(c0,t0,3) = u_star;
		C_UDMI(c0,t0,4) = u_star_ABL_calc;
		C_UDMI(c0,t0,5) = Cs;
		
		// Return Cs as a profile
		F_PROFILE(f,t,i) = Cs;
	}
	end_f_loop(f,t)

}

/*******************************************************************************/
/* CALCULATE Cs BASED ON Ks+ ACCORDING TO METHOD 3                             */
/*******************************************************************************/

DEFINE_PROFILE(wall_cs_M3,t,i)
{
	face_t f;
	cell_t c0;
	Thread *t0;
	real zp;
	real zp_plus;
	real z0;
	real ks_plus;
	real Cs;
	real u_p[ND_ND]; 
	real u_p_mag;
	real u_star;
	real u_star_ABL_calc;
	real phi_wall;
	real theta_wall;
	
	real xf[ND_ND];
	real dx[ND_ND]; 
	real dx_mag;
	real xc[ND_ND];
	
	int zone_ID;
	int index;
	
	zone_ID = THREAD_ID(t);
	
	// Loop through list of rough wall ID-s
	for(index = 0; index < n_rough_walls; index++)
	{
		if (zone_ID == rough_wall_ID[index])
		{
			// Assign z0 and Ks+ to the thread based on its ID
			z0 = z0_array[index];
			ks_plus = ks_plus_array[index];
		}
	}	

	// Loop through the faces of the thread
	begin_f_loop(f,t)
	{
		c0 = F_C0(f,t);
		t0 = THREAD_T0(t);
		
		F_CENTROID(xf, f, t);
		C_CENTROID(xc, c0, t0);
				
		dx[0] = xc[0] - xf[0];
		dx[1] = xc[1] - xf[1];
		dx[2] = xc[2] - xf[2];
		dx_mag = NV_MAG(dx);
		zp = dx_mag;
		// Calculate u* - see Equation (10)
		u_star = pow(Cmu,0.25)*pow(C_K(c0,t0),0.5);
		// Calculate zp+ - see Section 2.4
		zp_plus = zp*u_star/C_MU_L(c0,t0);
		
		u_p[0] = C_U(c0,t0);
		u_p[1] = C_V(c0,t0);
		u_p[2] = C_W(c0,t0);
		u_p_mag = NV_MAG(u_p);
		
		u_star_ABL_calc = u_p_mag*kappa/log((zp+z0)/z0);
		// Calculate exponent phi using Method 2 - see Equation (21)
		phi_wall = u_star_ABL/u_star;
		// Calculate theta - see Equation (25)
		theta_wall = 1./Ec*pow((zp+z0)/z0,phi_wall);
		
		// Calculate Cs - see Equation (24)
		if (ks_plus > zp_plus)
		{
			Cs = max((zp_plus+ks_plus/2-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		else
		{
			Cs = max((zp_plus-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		
		// Store diagnostic variables in UDMI
		C_UDMI(c0,t0,1) = phi_wall;
		C_UDMI(c0,t0,2) = zp_plus;
		C_UDMI(c0,t0,3) = u_star;
		C_UDMI(c0,t0,4) = u_star_ABL_calc;
		C_UDMI(c0,t0,5) = Cs;
		
		// Return Cs as a profile
		F_PROFILE(f,t,i) = Cs;
	}
	end_f_loop(f,t)

}

/*******************************************************************************/
/* CALCULATE Cs BASED ON Ks+ ACCORDING TO METHOD 4                             */
/*******************************************************************************/

DEFINE_PROFILE(wall_cs_M4,t,i)
{
	face_t f;
	cell_t c0;
	Thread *t0;
	real zp;
	real zp_plus;
	real z0;
	real ks_plus;
	real Cs;
	real u_p[ND_ND]; 
	real u_p_mag;
	real u_star;
	real u_star_ABL_calc;
	real phi_wall;
	real theta_wall;
	
	real xf[ND_ND];
	real dx[ND_ND]; 
	real dx_mag;
	real xc[ND_ND];
	
	int zone_ID;
	int index;
	
	zone_ID = THREAD_ID(t);
	
	// Loop through list of rough wall ID-s
	for(index = 0; index < n_rough_walls; index++)
	{
		if (zone_ID == rough_wall_ID[index])
		{
			// Assign z0 and Ks+ to the thread based on its ID
			z0 = z0_array[index];
			ks_plus = ks_plus_array[index];
		}
	}	

	// Loop through the faces of the thread
	begin_f_loop(f,t)
	{
		c0 = F_C0(f,t);
		t0 = THREAD_T0(t);
		
		F_CENTROID(xf, f, t);
		C_CENTROID(xc, c0, t0);
				
		dx[0] = xc[0] - xf[0];
		dx[1] = xc[1] - xf[1];
		dx[2] = xc[2] - xf[2];
		dx_mag = NV_MAG(dx);
		zp = dx_mag;
		// Calculate u* - see Equation (10)
		u_star = pow(Cmu,0.25)*pow(C_K(c0,t0),0.5);
		// Calculate zp+ - see Section 2.4
		zp_plus = zp*u_star/C_MU_L(c0,t0);
		
		u_p[0] = C_U(c0,t0);
		u_p[1] = C_V(c0,t0);
		u_p[2] = C_W(c0,t0);
		u_p_mag = NV_MAG(u_p);
		
		u_star_ABL_calc = u_p_mag*kappa/log((zp+z0)/z0);
		// Calculate exponent phi using Method 4 - see Equation (23)
		phi_wall = 1;
		// Calculate theta - see Equation (25)
		theta_wall = 1./Ec*pow((zp+z0)/z0,phi_wall);
		
		// Calculate Cs - see Equation (24)
		if (ks_plus > zp_plus)
		{
			Cs = Max((zp_plus+ks_plus/2-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		else
		{
			Cs = Max((zp_plus-theta_wall)/(theta_wall*ks_plus),Cs_min);
		}
		
		// Store diagnostic variables in UDMI
		C_UDMI(c0,t0,1) = phi_wall;
		C_UDMI(c0,t0,2) = zp_plus;
		C_UDMI(c0,t0,3) = u_star;
		C_UDMI(c0,t0,4) = u_star_ABL_calc;
		C_UDMI(c0,t0,5) = Cs;
		
		// Return Cs as a profile
		F_PROFILE(f,t,i) = Cs;
	}
	end_f_loop(f,t)

}


/*******************************************************************************/
/* BLENDING FUNCTION AND ITS DERIVATIVES - SEE EQUATION (8)                    */
/*******************************************************************************/

real f_blend(real z)
{
	return 0.5*(1.0+sin(M_PI/2.0*max(-1.0,min(alpha_b*(z-delta),1.0))));
}
/*******************************************************************************/

real df_blend_dz(real z)
{
	if (z > (delta-1/alpha_b) && z < (delta+1/alpha_b))
		{
			return M_PI/4.0*alpha_b*cos(M_PI/2.0*alpha_b*(z-delta));
		}
	else
		{
			return 0.0;
		}
}

/*******************************************************************************/
/* MODIFIED VERTICAL COORDINATE AND ITS DERIVATIVE - SEE EQUATION (7)          */
/*******************************************************************************/

real zmod(real z)
{
	return z+(delta-z)*f_blend(z);
}

/*******************************************************************************/

real dzmoddz(real z)
{
	return 1-f_blend(z)+(delta-z)*df_blend_dz(z);
}

/*******************************************************************************/
/* VELOCITY POWER LAW PROFILE AND ITS DERIVATIVE - SEE EQUATION (3)            */
/*******************************************************************************/

real vel_in(real z)
{
    return u_alpha*pow(((z+z_alpha)/z_alpha),alpha)-u_alpha;
}

/*******************************************************************************/

real pdudz(real z)
{
	return u_alpha*alpha/z_alpha*pow(((z+z_alpha)/z_alpha),alpha-1);
}

/*******************************************************************************/

real dudz(real z)
{
	real zm = zmod(z);
	return pdudz(zm)*dzmoddz(z);
}

/*******************************************************************************/
/* lOGARITHMIC-POLYNOMIAL TKE INFLOW PROFILE - SEE EQUATION (4)                */
/*******************************************************************************/

real tke_in(real z)
{
    real zf = (z + z0k)/z0k;
    real logzf = log(zf);
    return (A_k*pow(logzf,2) + B_k*pow(logzf,3) + C_k*pow(logzf,4) + 1)*D_k;
}

/*******************************************************************************/
/* EPSILON INFLOW PROFILE                                                      */
/*******************************************************************************/

real epsilon_in(real z)
{
    real du_dzc = alpha*u_alpha*pow(((z+z_alpha)/z_alpha),alpha)/(z+z_alpha);
    real k = tke_in(z);
    
    return sqrt(Cmu)*du_dzc*k;
}

/*******************************************************************************/
/* OMEGA INFLOW PROFILE - SEE EQUATION (6)                                     */
/*******************************************************************************/

real omega_in(real z)
{
    return pdudz(z)/sqrt(Cmu);
}