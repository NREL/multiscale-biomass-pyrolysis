#ifndef _BOUNDARIES_H_
#define _BOUNDARIES_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Species.H>
#include <VarDefines.H>

using namespace amrex;
namespace tr_boundaries
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void bc_ib(IntVect face,
               int dir, int sgn,
               int solved_comp,
               Array4<Real> const& phi,
               Array4<Real> const& acoeff,
               Array4<Real> const& rhs,
               GpuArray<int, AMREX_SPACEDIM> domlo,
               GpuArray<int, AMREX_SPACEDIM> domhi,
               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
               GpuArray<Real, AMREX_SPACEDIM> dx,
               const Real time,
               ProbParm const& prob_parm,int conjsolve)
    {
        //sgn is outward normal with the usual conjsolve=0
        //solve
        IntVect lcell=face;
        IntVect rcell=face;
        lcell[dir]-=1;
        amrex::Real posvec[2];
        posvec[0]=prob_lo[0]+(lcell[0]+0.5)*dx[0];
        posvec[1]=prob_lo[1]+(lcell[1]+0.5)*dx[1];
        posvec[dir]+=0.5*dx[dir];

        //these are implicit_solve_scalar
        //mlabec constants used
        amrex::Real bscalar=1.0;
        amrex::Real betasign=1.0;

        //interior cell
        IntVect intcell;
        IntVect outcell;

        if(conjsolve==0)
        {
            intcell=(sgn==1)?lcell:rcell;
            outcell=(sgn==1)?rcell:lcell;
        }
        else
        {
            intcell=(sgn==1)?rcell:lcell;
            outcell=(sgn==1)?lcell:rcell;
        }

        IntVect intcell1=intcell;
        int intcell1_adjust=(conjsolve==0)?-sgn:sgn;
        intcell1[dir]+=intcell1_adjust;
        
        IntVect outcell1=outcell;
        int outcell1_adjust=(conjsolve==0)?sgn:-sgn;
        outcell1[dir]+=outcell1_adjust;
        IntVect outcell2=outcell1;
        outcell2 += outcell1_adjust;

        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=(conjsolve==0)?sgn:-sgn;
        amrex::Real dx2=dx[dir]*dx[dir];

        //
        amrex::Real Temp_D=2.0*phi(outcell1,TEMP_ID)-phi(outcell2,TEMP_ID);
        amrex::Real dTdn=(phi(outcell2,TEMP_ID)-phi(outcell1,TEMP_ID))/dx[dir]*(-sgn);
        amrex::Real dcoeff_inv=phi(outcell1,MAT1_ID)/prob_parm.dc_mat1;
        dcoeff_inv += phi(outcell1,MAT2_ID)/prob_parm.dc_mat2;
        dcoeff_inv += phi(outcell1,MAT3_ID)/prob_parm.dc_mat3;
        amrex::Real dcoeff=1.0/dcoeff_inv;
        amrex::Real outflux=betasign*dcoeff*dTdn;

        /*if(conjsolve==0)
        {
            amrex::Print()<<"Temp_D:"<<Temp_D<<"\t"<<conjsolve<<"\n";
            amrex::Print()<<"intcell,outcell,face:"<<intcell<<"\t"<<outcell<<"\t"<<face<<"\t"
            <<phi(intcell,CMASK_ID)<<"\t"<<phi(outcell,CMASK_ID)<<"\n";
        }
        else
        {
            amrex::Print()<<"outheatflux:"<<outheatflux<<"\t"<<conjsolve<<"\n";
            amrex::Print()<<"intcell,outcell,face:"<<intcell<<"\t"<<outcell<<"\t"<<face<<"\t"
            <<phi(intcell,CMASK_ID)<<"\t"<<phi(outcell,CMASK_ID)<<"\n";
        }*/

        if(solved_comp==TEMP_ID)
        {
            if(conjsolve==0)
            {
                amrex::Real dc_inv=phi(intcell,MAT1_ID)/prob_parm.dc_mat1;
                dc_inv += phi(intcell,MAT2_ID)/prob_parm.dc_mat2;
                dc_inv += phi(intcell,MAT3_ID)/prob_parm.dc_mat3;
                amrex::Real dc=1.0/dc_inv;
                acoeff(intcell) += bscalar*betasign*dc/(0.5*dx2);
                rhs(intcell) += bscalar*betasign*dc*Temp_D/(0.5*dx2);
            }
            else
            {
                rhs(intcell)   += bscalar*outflux/dx[dir];
            }
        }
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn, 
                    int spec_id, Array4<Real> const &phi, 
                    Array4<Real> const& bc_arr,
                    Array4<Real> const& robin_a,
                    Array4<Real> const& robin_b,
                    Array4<Real> const& robin_f,
                    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                    GpuArray<Real, AMREX_SPACEDIM> dx,
                    const Real time,
                    ProbParm const& prob_parm)

    {
        //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        //default
        robin_a(ghost_cell) = 0.0;
        robin_b(ghost_cell) = 1.0;
        robin_f(ghost_cell) = 0.0;
        bc_arr(ghost_cell) = 0.0;

        if(sgn == -1) 
        { 
            // lo sides
            if(spec_id==TEMP_ID)
            {
                robin_a(ghost_cell) = 1.0;
                robin_b(ghost_cell) = 0.0;
                robin_f(ghost_cell) = 1.0;
                bc_arr(ghost_cell) = 1.0;
            }
        }
        else
        {
            if(spec_id==TEMP_ID)
            {
                robin_a(ghost_cell) = 1.0;
                robin_b(ghost_cell) = 0.0;
                robin_f(ghost_cell) = 0.0;
                bc_arr(ghost_cell) = 0.0;
            }
        } 

    }
}
#endif
