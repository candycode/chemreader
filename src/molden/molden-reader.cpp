#include <parsley/types.h>
#include <parsley/InStream.h>
#include <parsley/parsers.h>
#include <parsley/parser_operators.h>
#include <parsley/ParserManager.h>
#include "../../include/typedefs.h"
#include "../../include/chemreader.h"

namespace chr {

using namespace parsley;

static SkipBlankParser __;
static BlankParser _;
static SkipToNextLineParser endl_;
static EofParser eof_;
static const bool SKIP_BLANKS = true;
static const bool DONT_SKIP_BLANKS = !SKIP_BLANKS;
static const bool SKIP_TO_NEXT_LINE = true;
static const bool DONT_SKIP_TO_NEXT_LINE = !SKIP_TO_NEXT_LINE;

enum { START, 
	       SKIP_LINE, 
	       ATOMS, 
	           ATOM_DATA, 
	       GTO,
               GTO_ATOM, 
                   GTO_ATOM_ORBITAL, 
                       GTO_ATOM_ORBITAL_COEFF_EXP,
           MOL_ORBITAL,
               MOL_ORBITAL_SYM,
               MOL_ORBITAL_ENE,
               MOL_ORBITAL_SPIN,
               MOL_ORBITAL_OCC,
               MOL_ORBITAL_COEFFS   
       INVALID_STATE };

///Shared state
struct State {
	std::vector< std::pair< StateID, StateID > >
		transitions;
    unsigned basisFunctionNum_;
    unsigned basisFunctionCounter_;
    unsigned line_;
    unsigned atoms_;
    unsigned shells_;
    State() : basisFunctionNum_( 0 ), 
              basisFunctionCounter_( 0 ), 
              line_( 0 ),
              atoms_( 0 ),
              shells_( 0 ) {}
} 

struct DefaultCBack : TransitionCBackDefault {
    State* state_;
    DefaultCBack( State* s ) : state_( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
       state_->transitions.push_back( std::make_pair( prev, cur ) );       
    }
    DefaultCBack* Clone() const { return new DefaultCBack( *this ); }
};
struct AtomDataCBack : DefaultCBack {
    AtomDataCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
    	state_->dataFrame.atomElements.push_back( 
    		AtomicNumFromElementName( v[ "element name" ] ) );
    	state_->dataFrame.atomCoord.push_back( 
    		Vec3D( v[ "x" ], v[ "y"], v[ "z" ] ) );
        DefaultCBack::Apply( v, prev, cur );
    }
    AtomDataCBack* Clone() const { return new AtomDataCBack( *this ); }
};
struct AtomGTOCBack : DefaultCBack {
    AtomGTOCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
    	if( state_->atomicOrbitalShells.empty() ) {
    		state_->atomicOrbitalShells.resize( state_.atomElements.size() ); 
    	}
    	state_->curAtomicOrbitalShell = v[ "num" ];
        DefaultCBack::Apply( v, prev, cur );
    }
    AtomGTOCBack* Clone() const { return new AtomGTOCBack( *this ); }
};
struct AddGTOAtomOrbitalCBack : DefaultCBack {
    AddGTOAtomOrbitalCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
    	state_->curAtomicOrbitalShell = v[ "num" ];
        DefaultCBack::Apply( v, prev, cur );
    }
    AddGTOAtomOrbitalCBack* Clone() const { 
    	return new AddGTOAtomOrbitalCBack( *this ); 
    }
};

DataFrame ReadMolden( const char* fname ) {

    std::ifstream ifs( fname );
    InStream is( ifs );


    ParserManager pM;
    pM.SetBeginEndStates( START, EOF_ );   
    state;
    
    //DEFINE STATE TRANSITIONS 
    pM.AddStates< INVALID_STATE >()
       ( START, 
             ATOMS )
       ( ATOMS,
             ATOM_DATA)
       ( ATOM_DATA,
       	     ATOM_DATA, GTO )
       ( GTO, 
             GTO_ATOM )
       ( GTO_ATOM,
       	     GTO_ATOM_ORBITAL )
       ( GTO_ATOM_ORBITAL,
       	     GTO_ATOM_ORBITAL_COEFF_EXP )
       ( GTO_ATOM_ORBITAL_COEFF_EXP,
       	     GTO_ATOM_ORBITAL_COEFF_EXP, 
       	     GTO_ATOM_ORBITAL
       	     MOL_ORBITAL )
       ( MOL_ORBITAL,
       	     MOL_ORBITAL_SYM )
       ( MOL_ORBITAL_SYM,
       	     MOL_ORBITAL_ENE )
       ( MOL_ORBITAL_ENE,
             MOL_ORBITAL_SPIN )
       ( MOL_ORBITAL_SPIN,
       	     MOL_ORBITAL_OCC )
       ( MOL_ORBITAL_OCC,
       	     MOL_ORBITAL_COEFFS )
       ( MOL_ORBITAL_COEFFS,
       	     MOL_ORBITAL_COEFFS )
       ( MOL_ORBITAL_COEFFS,
       	     MOL_ORBITAL_SYM, EOF, EMPTY_LINE )

    pM.SetAllTransitionsCBack( new DefaultCBack( &state ) );
    // assign cback to *ALL* transitions whose target is GTO_COEFF_EXP   
    pM.SetCBacks()
       ( ATOM_DATA, new AddAtomCBack( state ) )
       ( GTO_ATOM, new AddGTOAtomCBack( state ) )
       ( GTO_ATOM_ORBITAL, new AddGTOAtomOrbitalCBack( state ) )
       ( GTO_ATOM_ORBITAL_COEFF_EXP, new AddOrbCoeffExp( state ) )
       ( MOL_ORBITAL_SYM, new AddNewMolOrbCBack( state ) )
       ( MOL_ORBITAL_ENE, new AddNewMolOrbCBack( state ) )
       ( MOL_ORBITAL_SPIN, new AddNewMolOrbCBack( state ) )
       ( MOL_ORBITAL_OCC, new AddNewMolOrbCBack( state ) )
       ( MOL_ORBITAL_COEFFS, new AddNewMolOrbCoeff( state ) )
#endif       
    //SPECIFY PER-STATE PARSERS   
    TupleParser<> coord( FloatParser(), "coord", __, _, endl_ );
    pM.SetParsers()
        ( MOLDEN_FORMAT, ( C("[") > C("Molden Format") > C("]") ) )
        ( START, NotParser< AL >( C("[") > C("Molden Format") > C("]") ) )
        ( SKIP_LINE, endl_ )
        ( TITLE, ( C("[") > C("Title") > C("]") ) )
        ( TITLE_DATA, ( __ > NC("\n") ) )
        ( ATOMS, AL(DONT_SKIP_BLANKS) >= __ >= C("[") >= __ >= C("Atoms") 
            >= __ >= C("]") >= __ >= A("unit") >= endl_ );
    AL atomDataParser(false);
    atomDataParser >= __ >= A("name") >= _ >= U("#") >= _  >= U("element") 
                   >= _ >= F("x") >= _ >= F("y") >= _ >= F("z") >= endl_;
    pM.SetParsers()
        ( ATOM_DATA, atomDataParser )
        ( GTO, ( C("[") > C("GTO") > C("]") ) )
        ( GTO_ATOM_NUM, AL( DONT_SKIP_BLANKS ) >= __ 
                        >= U("gto atom #") >= endl_  )
        ( GTO_SHELL_INFO, ( AL( DONT_SKIP_BLANKS ) >= __ 
                           >= FA("shell") >= _ >= U("num basis" ) >= endl_ ) )
        ( GTO_COEFF_EXP, ( AL( DONT_SKIP_BLANKS ) >= __ >= F("exp") 
                           >= _ >= F("coeff") >= endl_ ) )
        ( EOF_, ( __ > eof_ ) )
        ( EOL_, endl_ );
    
    //APPLY TO FIRST INPUT
    pM.Apply( is, START );
   
}

} //namespace