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

enum { START, MOLDEN_FORMAT, TITLE, TITLE_DATA, ATOMS, ATOM_DATA, GTO,
       GTO_ATOM_NUM, GTO_SHELL_INFO, GTO_COEFF_EXP, EOL_, EOF_, SKIP_LINE, 
       INVALID_STATE };

///Shared state
struct State {
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
        std::cout << state_->line_ << "> [" << cur << "]: ";
        Values::const_iterator i = v.begin();
        for( ; i != v.end(); ++i ) {
            std::cout << '(' << i->first << " = " << i->second << ')' << ' ';
        }
        std::cout << '\n';
    }
    DefaultCBack* Clone() const { return new DefaultCBack( *this ); }
};
struct AtomDataCBack : DefaultCBack {
    AtomDataCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
        DefaultCBack::Apply( v, prev, cur );
        ++state_->atoms_;
    }
    AtomDataCBack* Clone() const { return new AtomDataCBack( *this ); }
};
struct AtomNumCBack : DefaultCBack {
    AtomNumCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
        DefaultCBack::Apply( v, prev, cur );
        ++state_->shells_;
    }
    bool Enabled(  const Values&, StateID prev, StateID next ) {
        if( state_->atoms_ == state_->shells_ ) return false;
        return true; 
    }
    AtomNumCBack* Clone() const { return new AtomNumCBack( *this ); }
};
struct CoeffExpCBack : DefaultCBack { 
    CoeffExpCBack( State* s ) : DefaultCBack( s ) {}
    void Apply( const Values& v, StateID prev, StateID cur ) {
        DefaultCBack::Apply( v, prev, cur );
        ++state_->basisFunctionCounter_;
    }
    bool Enabled( const Values& , StateID prev, StateID next ) const {
        if( prev == GTO_COEFF_EXP 
            && state_->basisFunctionCounter_ == state_->basisFunctionNum_ ) {
            return false;
        }
        return true; 
    }
    CoeffExpCBack* Clone() const { return new CoeffExpCBack( *this ); } 
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
           MOLDEN_FORMAT, SKIP_LINE )
       ( SKIP_LINE,
           MOLDEN_FORMAT, START )
       ( MOLDEN_FORMAT,
           TITLE )
       ( TITLE, ATOMS, 
           TITLE_DATA )
       ( TITLE_DATA, 
           ATOMS )
       ( ATOMS, 
           ATOM_DATA )
       ( ATOM_DATA, 
           ATOM_DATA, GTO)
       ( GTO, 
           GTO_ATOM_NUM )
       ( GTO_ATOM_NUM, 
           GTO_SHELL_INFO )
       ( GTO_SHELL_INFO, 
           GTO_COEFF_EXP )
       ( GTO_COEFF_EXP, /* from */
         /*to*/ GTO_SHELL_INFO, /*or*/ GTO_COEFF_EXP, 
         /*or*/ GTO_ATOM_NUM, /*or*/ EOF_, EOL_ );
    pM.SetAllTransitionsCBack( new DefaultCBack( &state ) );
    // assign cback to *ALL* transitions whose target is GTO_COEFF_EXP   
    pM.SetCBacks()
       ( GTO_COEFF_EXP, new CoeffExpCBack( &state ) );
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