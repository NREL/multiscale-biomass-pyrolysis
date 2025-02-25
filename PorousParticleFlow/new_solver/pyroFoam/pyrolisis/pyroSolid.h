#pragma once

#include "fvCFD.H"
#include <utility>

namespace Foam
{

namespace pyrolisis
{

class irreversibleArrheniusReaction;

class pyroSolid
{

private:

    const fvMesh& m_mesh;
    IOdictionary m_dict;
    volScalarField m_porosity;
    label m_nSubTimeSteps;
    scalar m_poreSize;
    List<word> m_speciesName;
    PtrList<volScalarField> m_species;
    PtrList<volScalarField> m_reactionRates;
    scalarField m_rho;
    scalarField m_cp;
    scalarField m_molWeight;
    scalarField m_kappa;
    List<std::pair<bool,label>> m_is_gas; 
    PtrList<irreversibleArrheniusReaction>  m_reactions;


    label getSpecieId(const word& name);

public:

    pyroSolid(const fvMesh& mesh);
    ~pyroSolid();

    void evolve();

    const volScalarField& porosity();

    bool isGasTransferSpecie(const word& name);

    const volScalarField& RR(const word& name);

    const scalar& poreSize();
};

}
}
