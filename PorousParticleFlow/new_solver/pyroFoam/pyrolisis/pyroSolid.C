#include "pyroSolid.h"
#include "irreversibleArrheniusReaction.h"

using namespace Foam;
using namespace pyrolisis;

pyroSolid::pyroSolid(const fvMesh& mesh)
:
    m_mesh(mesh),
    m_dict (
        IOobject
        (
            "pyrolisisDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    m_porosity
    (
        IOobject
        (
            "porosity",
            m_mesh.time().timeName(),
            m_mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        m_mesh
    ),
    m_nSubTimeSteps(m_dict.lookupOrDefault<label>("nSubTimeSteps",1)),
    m_poreSize(readScalar(m_dict.lookup("poreSize"))),
    m_speciesName(m_dict.lookup("species")),
    m_species(m_speciesName.size()),
    m_rho(m_speciesName.size()),
    m_cp(m_speciesName.size()),
    m_molWeight(m_speciesName.size()),
    m_kappa(m_speciesName.size()),
    m_is_gas(m_speciesName.size())

{
    Info << "Reading pyrolisis model" << endl;
    Info << "Found species: " << m_speciesName << endl;

    forAll(m_speciesName, specieI)
    {

        dictionary& speciesDict(m_dict.subDict("speciesCoeffs").subDict(m_speciesName[specieI]));
        m_species.set (
            specieI,
            new volScalarField
            (
                IOobject
                (
                    m_speciesName[specieI] + ".solid",
                    m_mesh.time().timeName(),
                    m_mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                m_mesh,
                dimensionedScalar("Y",dimless,0.)
            )
        );

        m_rho[specieI] = scalar(readScalar(speciesDict.lookup("rho")));
        m_cp[specieI] = scalar(readScalar(speciesDict.lookup("cp")));
        m_kappa[specieI] = scalar(readScalar(speciesDict.lookup("kappa")));
        m_molWeight[specieI] = scalar(readScalar(speciesDict.lookup("molWeight")));
        m_molWeight[specieI] *= 1e-3;
        m_is_gas[specieI].first = speciesDict.lookupOrDefault<bool>("gas",false);

        if (m_is_gas[specieI].first)
        {
            m_is_gas[specieI].second = m_reactionRates.size();

            m_reactionRates.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        m_speciesName[specieI]  + ".RR",
                        m_mesh.time().timeName(),
                        m_mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    m_mesh,
                    dimensionedScalar("Rdot",dimDensity/dimTime,0.)
                )            
            );
        }
    }


    List<word> reactionList(m_dict.lookup("reactions"));

    Info << "Found reactions: " << reactionList << endl;
    
    m_reactions.resize(reactionList.size());

    forAll(reactionList, reactI)
    {
        if( !m_dict.subDict("reactionCoeffs").found(reactionList[reactI]) )
        {
            FatalErrorInFunction << "Cannot find reaction " << reactionList[reactI] << "\n" << abort(FatalError);            
        }

        m_reactions.set (
            reactI,
            new irreversibleArrheniusReaction
            (
                reactionList[reactI],
                m_dict.subDict("reactionCoeffs").subDict(reactionList[reactI]),
                m_speciesName
            )
        );        
    }
}

pyroSolid::~pyroSolid()
{}

void pyroSolid::evolve()
{
    Info << "Updating solid composition" << endl;

    const volScalarField& T = m_mesh.lookupObject<volScalarField>("T");

    scalarField Y(m_species.size(),0.);

    /* Reset reaction rate */
    forAll(m_species, specieI)
    {
        if( m_is_gas[specieI].first )
        {
            m_reactionRates[m_is_gas[specieI].second] *= 0.;
        }
    }

    /* Go cell-by-cell */
    forAll(m_mesh.C(), cellI)
    {
        // Skip if outside the solid
        if ( m_porosity[cellI] > 0.999 )
        {
            continue;
        }

        scalar dt = m_mesh.time().deltaT().value();

        scalar sub_dt = dt / m_nSubTimeSteps;

        forAll(m_species, specieI)
        {
            //- Convert from [vol_specie/m3] to [mol_specie/m3] by multiplying by density and dividing
            //  by molar weight.
            //  Always use oldTime to make it work with pimple outer iterations.
            Y[specieI] =  m_species[specieI].oldTime()[cellI] * ( m_rho[specieI] / m_molWeight[specieI] );
        }

        /* Use Euler marching */
        for (label ts = 0; ts < m_nSubTimeSteps; ts++)
        {
            scalarField ndot(m_species.size(),0.);

            forAll(m_reactions, reactI)
            {
                ndot += m_reactions[reactI].computeMolarSources(T[cellI], Y);
            }

            Y += sub_dt * ndot;

            Y = max(Y,0.);
        }

        /* Update reaction rate and concentration*/
        forAll(m_species, specieI)
        {
            scalar Y_0 = m_species[specieI].oldTime()[cellI] * m_rho[specieI];
            m_species[specieI][cellI] = Y[specieI] * m_molWeight[specieI] / m_rho[specieI];
            
            if( m_is_gas[specieI].first )
            {
                m_reactionRates[m_is_gas[specieI].second][cellI] = ( (Y[specieI] *  m_molWeight[specieI] ) - Y_0 )  / dt;
                m_species[specieI][cellI] = 0.;             
            }
        }

        /* Update local porosity */
        scalar vol_species(0.);
        forAll(m_species, specieI)
        {
            vol_species += m_species[specieI][cellI];
        }
        
        m_porosity[cellI] = 1.0 - vol_species;
    }
}

const volScalarField& pyroSolid::porosity()
{
    return m_porosity;
}

label pyroSolid::getSpecieId(const word& name)
{
    label id = -1;
    forAll(m_speciesName, specieI)
    {
        if( m_speciesName[specieI] == name )
        {
            id = specieI;
        }
    }

    return id;
}


bool pyroSolid::isGasTransferSpecie(const word& name)
{
    label id = getSpecieId(name);

    if (id == -1)
    {
        FatalErrorInFunction << "Error: unknown gas specie " << name << "\n" << abort(FatalError);
    } 

    return m_is_gas[id].first;
}

const volScalarField& pyroSolid::RR(const word& name)
{
    label id = getSpecieId(name);

    if (id == -1)
    {
        FatalErrorInFunction << "Error: unknown gas specie " << name << "\n" << abort(FatalError);
    } 

    return m_reactionRates[m_is_gas[id].second];
}

const scalar& pyroSolid::poreSize()
{
    return m_poreSize;
}