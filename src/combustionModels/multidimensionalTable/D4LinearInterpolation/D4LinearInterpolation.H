/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

Class
    Foam::combustionModels::linearInter

Description
    Linear Table Interpolation

\*---------------------------------------------------------------------------*/

#ifndef D4LINEARINTERPOLATION_H_
#define D4LINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "multidimensionalTable.H"

namespace Foam
{

class fvMesh;

class D4LinearInterpolation
:
    public multidimensionalTable
{

private:

	//- Disallow default bitwise copy construct
    D4LinearInterpolation(const D4LinearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const D4LinearInterpolation&);

    //- Table Values
    List<List<List<scalarList> > > tableValues_ ;

public:

    static const List<List<List<scalarList> > > defaultList;

    //- Runtime type information
    TypeName("D4LinearInterpolation");

    // Constructors

        //- Construct from components
        D4LinearInterpolation(const fvMesh& mesh, const word& tableName);

        //- Destructor
        virtual ~D4LinearInterpolation();

    //- Member functions

        //- Interpolate tableValues
        virtual inline scalar interpolate(const List<int>& ub, const scalarList& pos) const;

        //- Return table Size
        virtual inline List<int> tableSize() const;
};

} // End Foam namespace

#endif /* D4LINEARINTERPOLATION_H_ */
