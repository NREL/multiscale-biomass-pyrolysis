volTensorField Kp
(
    IOobject
    (
        "Kp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("Kp",dimensionSet(1,-3,-1,0,0,0,0),tensor::zero)
);

IOdictionary Kdict
(
    IOobject
    (
        "permeabilityDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


word permeabilityModel(Kdict.lookup("Model"));

if(permeabilityModel == "constant")
{
    dimensionedTensor Kvalue("Kvalue",Kp.dimensions(),Kdict.lookup("value"));

    Kp = Kvalue;
}