#pragma once
////////////////////////////////////////////////////////////////////////////////
/*
Copyright (c) 2013, Ugo Varetto
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL UGO VARETTO BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cassert>
#include <string>
#include <map>

namespace chr {

typedef real_t double;
typedef std::string String

struct Point3D {
    real_t x, y, z;
};

enum UnitOfLength {
    ATOMIC_UNIT,
    ANGSTROM
};

enum AtomicOrbitalType {
	s, sp, p, d, f
};

enum BasisFunctionType {
	GUSSIAN,
	SLATER
};

struct AtomicOrbital {
	AtomicOrbitalType type; //s, sp, p...
	BaisFunctionType basisType; // gaussian, slater
	String basisSetId; //e.g. 3-21G (6D, 7F) 
	std::vector< real_t > coeff;
	std::vector< real_t > exp;
};

struct AtomicOrbitalShell {
    std::vector< AtomicOrbital > orbitals;
    int atomid;
};

enum MolecularOrbitalType {
	MO_SIGMA, MO_PI, MO_UNKNOWN
};

enum SpinType {
    SPIN_ALPHA,
    SPIN_BETA,
    SPIN_ALPHA_BETA
};

struct MolecularOrbital {
    MolecularOrbitalType type; //MO_SIGMA, MO_PI, MO_UNKNOWN
    String name;    
    real_t energy;
    SpinType spin;
    real_t occupancy;
    std::vector< real_t > coeffs;
};

struct Bond {
    int atomId1;
    int atomId2;
    String type;
    String text;
    int multiplicity;
};

// single time/optimization step data
struct DataFrame {
	String dataFormat; //molden, gaussian09 ...
	String source; //file 
	///@todo add support of per-atom scalar values
	std::vector< int > atomElements;
	std::vector< String > atomNames;
	UnitOfLength unit; //AU or Angstrom
	std::vector< Point3D > atomCoord;
	std::vector< Bond > bonds;
	std::vector< AtomicOrbitalShell > atomicOrbitalShell;
	std::vector< MolecularOrbital > molecularOrbitals;
};

template < typename SequenceT >
struct DataSet {
    std::map< SequenceT, DataFrame > data; 
};


} //namespace
