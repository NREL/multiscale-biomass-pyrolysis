/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

References
    The overall approach is based on:
    Gidaspow, D. (1994). Multiphase flow and fluidisation: continuum and 
    kinetic theory descriptions. Academic press, New York.

    And on the method presented in:
    L.Lu, X. Gao, M. Shahnam, W. A. Rogers, "Coarse grained computational 
    fluid dynamic simulation of sands and biomass fluidization with a hybrid 
    drag", Particle Technology and Fluidization, vol 66 4 (2020)

    With extension to other drag laws in the low particle volume fraction
    regime.

Author
    Federico Municchi NLR (2026)


\*---------------------------------------------------------------------------*/

#ifndef HybridErgunBiomass_H
#define HybridErgunBiomass_H

#include "dispersedDragModel.H"
#include "Ergun.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                     Class HybridErgunBiomass Declaration
\*---------------------------------------------------------------------------*/

class HybridErgunBiomass
:
    public dispersedDragModel
{
    // Private Data

        //- Ergun drag model
        Ergun Ergun_;

        //- Model for biomass
        autoPtr<dispersedDragModel> Biomass_;

public:

    //- Runtime type information
    TypeName("HybridErgunBiomass");


    // Constructors

        //- Construct from a dictionary and an interface
        HybridErgunBiomass
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~HybridErgunBiomass();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
