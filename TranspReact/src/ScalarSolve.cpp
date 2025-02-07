#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <TranspReact.H>
#include <Species.H>
#include <Transport.H>
#include <UserBCs.H>
#include <Chemistry.H>
#include <AMReX_MLABecLaplacian.H>

void TranspReact::update_rxnsrc_at_all_levels(Vector<MultiFab>& Sborder,
                                         Vector<MultiFab>& rxn_src, 
                                         amrex::Real cur_time)
{
    amrex::Real time = cur_time;
    ProbParm const* localprobparm = d_prob_parm;

    // Zero out reactive source MFs
    for(int lev=0; lev <= finest_level; lev++)
    {
        rxn_src[lev].setVal(0.0);
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();

        for (MFIter mfi(rxn_src[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            Array4<Real> rxn_arr = rxn_src[lev].array(mfi);

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                tr_reactions::production_rate(i, j, k, sborder_arr, rxn_arr,
                 prob_lo, prob_hi, dx, time, *localprobparm);

            });
        }
    }
}

void TranspReact::implicit_solve_scalar(Real current_time, Real dt, int spec_id, 
                                   Vector<MultiFab>& Sborder, 
                                   Vector<MultiFab>& Sborder_old,
                                  Vector<MultiFab>& rxn_src) 
{
    BL_PROFILE("TranspReact::implicit_solve_species(" + std::to_string( spec_id ) + ")");


    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int linsolve_verbose = 1;
    int captured_spec_id=spec_id;
    int steady_solve=steadyspec[spec_id];

    //==================================================
    // amrex solves
    // read small a as alpha, b as beta

    //(A a - B del.(b del)) phi = f
    //
    // A and B are scalar constants
    // a and b are scalar fields
    // f is rhs
    // in this case: A=0,a=0,B=1,b=conductivity
    // note also the negative sign
    //====================================================
    ProbParm const* localprobparm = d_prob_parm;

    const Real tol_rel = linsolve_reltol;
    const Real tol_abs = linsolve_abstol;

    // set A and B, A=1/dt, B=1
    Real ascalar = 1.0;
    Real bscalar = 1.0;

#ifdef AMREX_USE_HYPRE
    if(use_hypre)
    {
        amrex::Print()<<"using hypre\n";
    }
#endif

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_lo 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_hi 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (all_bcs_lo[spec_id][idim] == PERBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (all_bcs_lo[spec_id][idim] == DIRCBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Dirichlet;
        }
        if (all_bcs_lo[spec_id][idim] == HNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }
        if (all_bcs_lo[spec_id][idim] == IHNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::inhomogNeumann;
        }
        if (all_bcs_lo[spec_id][idim] == ROBINBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
        if (all_bcs_lo[spec_id][idim] == AXISBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }

        //higher side bcs
        if (all_bcs_hi[spec_id][idim] == PERBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Periodic;
        }
        if (all_bcs_hi[spec_id][idim] == DIRCBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Dirichlet;
        }
        if (all_bcs_hi[spec_id][idim] == HNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
        if (all_bcs_hi[spec_id][idim] == IHNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::inhomogNeumann;
        }
        if (all_bcs_hi[spec_id][idim] == ROBINBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
        if (all_bcs_hi[spec_id][idim] == AXISBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
    }

    Vector<MultiFab> specdata(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<MultiFab> bcoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        specdata[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        bcoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    linsolve_ptr.reset(new MLABecLaplacian(Geom(0,finest_level), 
                                           boxArray(0,finest_level), 
                                           DistributionMap(0,finest_level), info));
    MLMG mlmg(*linsolve_ptr);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(linsolve_verbose);

#ifdef AMREX_USE_HYPRE
        if (use_hypre)
        {
            mlmg.setHypreOptionsNamespace("tr.hypre");
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        }
#endif
    linsolve_ptr->setDomainBC(bc_linsolve_lo, bc_linsolve_hi);
    linsolve_ptr->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        // Copy args (FabArray<FAB>& dst, FabArray<FAB> const& src, 
        // int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
        specdata[ilev].setVal(0.0);
        amrex::Copy(specdata[ilev], Sborder_old[ilev], captured_spec_id, 
                    0, 1, num_grow);

        acoeff[ilev].setVal(0.0);
        if(!steady_solve)
        {
           acoeff[ilev].setVal(1.0/dt);
        }

        bcoeff[ilev].setVal(1.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        rhs[ilev].setVal(0.0);
        amrex::MultiFab::Saxpy(rhs[ilev], 1.0, rxn_src[ilev], spec_id, 0, 1, 0);

        if(!steady_solve)
        {
            amrex::MultiFab::Saxpy(rhs[ilev], 1.0/dt, specdata[ilev], 0, 0, 1, 0);
        }

        amrex::Copy(specdata[ilev], Sborder[ilev], captured_spec_id, 
                0, 1, num_grow);

        solution[ilev].setVal(0.0);
        amrex::MultiFab::Copy(solution[ilev], specdata[ilev], 0, 0, 1, 0);

        // fill cell centered diffusion coefficients and rhs
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Real time = current_time; // for GPU capture

            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Array4<Real> bcoeff_arr = bcoeff[ilev].array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    bcoeff_arr(i,j,k)=tr_transport::specDiff(i,j,k,captured_spec_id,sb_arr,
                            dx,prob_lo,prob_hi,time,*localprobparm); 

                    });
        }



        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoeff[ilev].boxArray(), 
                    IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, bcoeff[ilev].DistributionMap(), 1, 0);
        }
        // true argument for harmonic averaging
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoeff), 
                bcoeff[ilev], geom[ilev], true);


        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> bc_arr = specdata[ilev].array(mfi);
            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Real time = current_time; // for GPU capture

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (!geom[ilev].isPeriodic(idim))
                {
                    //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                    //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                    //so the ghost cell index at left side is i-1 while it is i on the right
                    if (bx.smallEnd(idim) == domain.smallEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                tr_boundaries::species_bc(i, j, k, idim, -1, 
                                        captured_spec_id, sb_arr, bc_arr, robin_a_arr,
                                        robin_b_arr, robin_f_arr, 
                                        prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                    if (bx.bigEnd(idim) == domain.bigEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                tr_boundaries::species_bc(i, j, k, idim, +1, 
                                                             captured_spec_id, sb_arr, bc_arr, robin_a_arr, 
                                                             robin_b_arr, robin_f_arr,
                                                             prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                }
            }
        }

        linsolve_ptr->setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_ptr->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells
        if(mixedbc)
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]), &(robin_a[ilev]), 
                                     &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]));
        }
    }

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    //bound species density
    if(bound_specden)
    { 
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            amrex::Real minconc=min_species_conc; 
            for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> soln_arr = solution[ilev].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                      if(soln_arr(i,j,k) < minconc)
                      {
                        soln_arr(i,j,k)=minconc;
                      } 
                });
            }
        }
    }

    // copy solution back to phi_new
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, spec_id, 1, 0);
    }
    
    Print()<<"Solved species:"<<allvarnames[spec_id]<<"\n";

    //clean-up
    specdata.clear();
    acoeff.clear();
    bcoeff.clear();
    solution.clear();
    rhs.clear();

    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}

void TranspReact::transform_variables(Vector<MultiFab>& Sborder,amrex::Real cur_time)
{
    amrex::Real time = cur_time;
    ProbParm const* localprobparm = d_prob_parm;
    
    for(int lev=0;lev<=finest_level;lev++)
    {
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();

        for (MFIter mfi(phi_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            Array4<Real> phi_arr = phi_new[lev].array(mfi);

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                tr_reactions::transform(i, j, k, sborder_arr, phi_arr,
                 prob_lo, prob_hi, dx, time, *localprobparm);

            });
        }
    }
}


