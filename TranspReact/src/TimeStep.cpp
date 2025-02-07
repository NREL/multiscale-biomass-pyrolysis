#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <TranspReact.H>

// a wrapper for EstTimeStep
void TranspReact::ComputeDt(amrex::Real cur_time)
{
    dt[0]=fixed_dt;
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        dt[lev] = dt[lev - 1];
    }
}
