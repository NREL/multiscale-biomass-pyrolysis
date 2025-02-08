#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real left_temp = 500.0;
    amrex::Real right_temp = 300.0;
    int axial_dir=1;
    amrex::Real *xctdata = nullptr;
    amrex::Real xct_threshold=3300.0;
    int nx_xct=0;
    int ny_xct=0;
    int nz_xct=0;
    amrex::Real lo_xct_x=0.0;
    amrex::Real lo_xct_y=0.0;
    amrex::Real lo_xct_z=0.0;
    amrex::Real dx_xct_x=3e-6;
    amrex::Real dx_xct_y=3e-6;
    amrex::Real dx_xct_z=3e-6;
};
#endif
