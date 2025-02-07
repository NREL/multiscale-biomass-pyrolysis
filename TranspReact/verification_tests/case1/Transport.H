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
            Real dcoeff=1.0;
            return(dcoeff);
        }
}
#endif
