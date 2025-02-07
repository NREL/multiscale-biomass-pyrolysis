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
            Real dcoeff_inv=phi(i,j,k,MAT1_ID)/prob_parm.dc_mat1;
            dcoeff_inv += phi(i,j,k,MAT2_ID)/prob_parm.dc_mat2;
            dcoeff_inv += phi(i,j,k,MAT3_ID)/prob_parm.dc_mat3;
            
            return(1.0/dcoeff_inv);
        }
}
#endif
