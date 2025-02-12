#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#define S1_ID 0
#define PHI_ID 1
#define EFX_ID 2
#define EFY_ID 3
#define EFZ_ID 4
#define NUM_SPECIES 5

namespace tr_species
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
}

#endif
