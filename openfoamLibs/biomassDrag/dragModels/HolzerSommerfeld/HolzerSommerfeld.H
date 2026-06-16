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
    A. Holzer and M. Sommerfeld, "New simple correlation formula for the drag 
    coefficient of non-spherical particles", Powder Technology, 184 (2008) 361-365

Author
    Federico Municchi NLR (2026) 

\*---------------------------------------------------------------------------*/

#ifndef HolzerSommerfeld_H
#define HolzerSommerfeld_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                            Class HolzerSommerfeld Declaration
\*---------------------------------------------------------------------------*/

class HolzerSommerfeld
:
    public dispersedDragModel
{

private:

    //- Sphericity
    scalar sphericity_;

    // This is the ratio between the cross section of a spherical particle
    // and the actual cross section in the flow direction
    scalar projectedArea_;

public:

    //- Runtime type information
    TypeName("HolzerSommerfeld");


    // Constructors

        //- Construct from a dictionary and an interface
        HolzerSommerfeld
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~HolzerSommerfeld();


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
