#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>
#include <ProbParm.H>

using namespace amrex;
namespace tr_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
        Real specDiff(int i, int j, int k,
                      int spec_id,
                Array4<Real> const& phi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                const Real time, ProbParm const& prob_parm)
        {
            Real solidfrac=phi(i,j,k,SOLID_ID);
            Real voidfrac=1.0-solidfrac;
            
            //https://en.wikipedia.org/wiki/Thermal_diffusivity
            Real dcoeff_solid=0.082e-4;  //m2/s for wood
                                         
            //https://www.engineeringtoolbox.com/air-thermal-diffusivity-d_2011.html
            Real dcoeff_gas=22.39e-6; //m2/s for air

            Real dcoeff_inv=solidfrac/dcoeff_solid+voidfrac/dcoeff_gas;
            return(1.0/dcoeff_inv);
        }
}
#endif
