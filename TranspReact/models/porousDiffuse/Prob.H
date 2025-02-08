#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

    AMREX_INLINE
void amrex_probinit (ProbParm & h_prob_parm, ProbParm & d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    pp.query("left_temp", h_prob_parm.left_temp);
    pp.query("right_temp", h_prob_parm.right_temp );
    pp.query("axial_dir",h_prob_parm.axial_dir);
    
    std::string xctfile="pine_box.txt";
    pp.query("xctfile",xctfile);

    std::ifstream infile(xctfile.c_str());

    infile>>h_prob_parm.nx_xct>>h_prob_parm.ny_xct>>h_prob_parm.nz_xct;
    int ndata=h_prob_parm.nx_xct*h_prob_parm.ny_xct*h_prob_parm.nz_xct;

    amrex::Vector<amrex::Real> xctdatavec;
    xctdatavec.resize(ndata);

    for(int i=0;i<ndata;i++)
    {
        infile>>xctdatavec[i];
    }
    infile.close();
    
    h_prob_parm.xctdata  = (amrex::Real*) amrex::The_Arena()->alloc(ndata*sizeof(amrex::Real));
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                xctdatavec.begin(),
                                xctdatavec.end(),
                       h_prob_parm.xctdata);

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif

}

AMREX_GPU_DEVICE
    AMREX_INLINE
void initdomaindata(Box const& bx, Array4<Real> const& phi, 
        GeometryData const& geomdata,
        ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

    int i_xct,j_xct,k_xct;
    int axdir=prob_parm->axial_dir;
    Real dist_axial=prob_hi[axdir]-prob_lo[axdir];


#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real z = prob_lo[2] + (0.5+k) * dx[2];
            Real y = prob_lo[1] + (0.5+j) * dx[1];
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) 
                {
                    IntVect iv(i,j,k);
                    Real x = prob_lo[0] + (0.5+i) * dx[0];
                    RealVect centroid_coord(x,y,z);
                    
                    Real vfrac=0.0;
                    for(int kk=0;kk<2;kk++)
                    {
                        for(int jj=0;jj<2;jj++)
                        {
                            for(int ii=0;ii<2;ii++)
                            {
                                Real xx[3]={prob_lo[0]+(i+ii)*dx[0],
                                            prob_lo[1]+(j+jj)*dx[1],           
                                            prob_lo[2]+(k+kk)*dx[2]};

		                i_xct = floor((xx[0]-prob_parm->lo_xct_x)/prob_parm->dx_xct_x);
		                j_xct = floor((xx[1]-prob_parm->lo_xct_y)/prob_parm->dx_xct_y);
		                k_xct = floor((xx[2]-prob_parm->lo_xct_z)/prob_parm->dx_xct_z);

		                 if( (i_xct >= 0) && (i_xct < prob_parm->nx_xct)
			         && (j_xct >= 0)  && (j_xct < prob_parm->ny_xct)
			         && (k_xct >= 0)  && (k_xct < prob_parm->nz_xct)  )
		                 {
			            int index=k_xct*(prob_parm->nx_xct*prob_parm->ny_xct)
			            + j_xct*prob_parm->nx_xct + i_xct;
			
                                    if(prob_parm->xctdata[index] < prob_parm->xct_threshold)
			            {
		                        vfrac+=1.0;	                
			            }
		                }
                             }
                         }
                     }

                    vfrac*=0.125;
                    phi(i,j,k,SOLID_ID) = 1.0-vfrac;


                    phi(i,j,k,TEMP_ID)  = prob_parm->left_temp-(prob_parm->left_temp-prob_parm->right_temp)/dist_axial*centroid_coord[axdir];
                }
        }
    }
}
#endif
