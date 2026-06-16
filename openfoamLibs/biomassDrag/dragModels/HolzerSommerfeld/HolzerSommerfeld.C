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

\*---------------------------------------------------------------------------*/

#include "HolzerSommerfeld.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(HolzerSommerfeld, 0);
    addToRunTimeSelectionTable(dragModel, HolzerSommerfeld, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::HolzerSommerfeld::HolzerSommerfeld
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    sphericity_(readScalar(dict.lookup("sphericity"))),
    projectedArea_(readScalar(dict.lookup("projectedArea")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::HolzerSommerfeld::~HolzerSommerfeld()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::HolzerSommerfeld::CdRe() const
{
    const volScalarField alpha2
    (
        max(interface_.continuous(), interface_.continuous().residualAlpha())
    );

    const volScalarField Res(alpha2*interface_.Re());

    const scalar A = 8.0/sqrt(projectedArea_) + 16.0/sqrt(sphericity_);

    const scalar B = 3.0/pow(sphericity_, 0.75);

    const scalar C =
        pow
        (
            0.421,
            0.4*pow(-log(sphericity_), 0.2)
        )/projectedArea_;

    const volScalarField CdRe
    (
          A + B*sqrt(Res) + C*Res
    );


    return CdRe*pow(alpha2, -2.65);

}


// ************************************************************************* //
