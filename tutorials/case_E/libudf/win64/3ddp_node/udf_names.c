/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_EXECUTE_ON_LOADING(set_C_mu, libname);
extern DEFINE_EXECUTE_ON_LOADING(name_udmi, libname);
extern DEFINE_ON_DEMAND(Patch_ABL_k_e);
extern DEFINE_ON_DEMAND(Patch_ABL_k_w);
extern DEFINE_PROFILE(inlet_U,t,i);
extern DEFINE_PROFILE(inlet_V,t,i);
extern DEFINE_PROFILE(inlet_W,t,i);
extern DEFINE_PROFILE(inlet_k,t,i);
extern DEFINE_PROFILE(inlet_epsilon,t,i);
extern DEFINE_PROFILE(inlet_omega,t,i);
extern DEFINE_PROFILE(wall_ks,t,i);
extern DEFINE_PROFILE(wall_cs_M1,t,i);
extern DEFINE_PROFILE(wall_cs_M2,t,i);
extern DEFINE_PROFILE(wall_cs_M3,t,i);
extern DEFINE_PROFILE(wall_cs_M4,t,i);
__declspec(dllexport) UDF_Data udf_data[] = {
{"set_C_mu", (void(*)())set_C_mu, UDF_TYPE_EXECUTE_ON_LOADING},
{"name_udmi", (void(*)())name_udmi, UDF_TYPE_EXECUTE_ON_LOADING},
{"Patch_ABL_k_e", (void(*)())Patch_ABL_k_e, UDF_TYPE_ON_DEMAND},
{"Patch_ABL_k_w", (void(*)())Patch_ABL_k_w, UDF_TYPE_ON_DEMAND},
{"inlet_U", (void(*)())inlet_U, UDF_TYPE_PROFILE},
{"inlet_V", (void(*)())inlet_V, UDF_TYPE_PROFILE},
{"inlet_W", (void(*)())inlet_W, UDF_TYPE_PROFILE},
{"inlet_k", (void(*)())inlet_k, UDF_TYPE_PROFILE},
{"inlet_epsilon", (void(*)())inlet_epsilon, UDF_TYPE_PROFILE},
{"inlet_omega", (void(*)())inlet_omega, UDF_TYPE_PROFILE},
{"wall_ks", (void(*)())wall_ks, UDF_TYPE_PROFILE},
{"wall_cs_M1", (void(*)())wall_cs_M1, UDF_TYPE_PROFILE},
{"wall_cs_M2", (void(*)())wall_cs_M2, UDF_TYPE_PROFILE},
{"wall_cs_M3", (void(*)())wall_cs_M3, UDF_TYPE_PROFILE},
{"wall_cs_M4", (void(*)())wall_cs_M4, UDF_TYPE_PROFILE},
};
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
