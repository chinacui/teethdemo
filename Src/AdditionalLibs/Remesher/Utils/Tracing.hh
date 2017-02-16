/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                       *
 *   $Author$                                                      *
 *   $Date$                   *
 *                                                                           *
\*===========================================================================*/




//=============================================================================
//
//  CLASS Tracing
//
//=============================================================================


#ifndef ACG_TRACING_HH
#define ACG_TRACING_HH


//== INCLUDES =================================================================

#include "StopWatch.hh"
#include <iostream>
#include <string>
#include "../Config/ACGDefines.hh"


//== NAMESPACES ===============================================================

namespace ACG {


//== HELPER MACROS ============================================================


#define ACG_TRACE(_text) ACG::Tracing acg_tracer__(_text)
#define ACG_TRACE_CMD(_text, _cmd) { ACG_TRACE(_text); cmd; }

#define ACG_TIMED_TRACE(_text) ACG::TimedTracing acg_tracer__(_text)
#define ACG_TIMED_TRACE_CMD(_text, _cmd) { ACG_TIMED_TRACE(_text); _cmd; }

#define ACG_TRACE_PROGRESS acg_tracer__.progress()


//== CLASS DEFINITION =========================================================


/** \class Tracing Tracing.hh <ACG/Utils/Tracing.hh>

    Tracing outputs starting and finishing messages for some lengthy
    procedures. Its constructor will output a given string, the
    destructor will output "finished\\n". That means you have to
    start a new block and define a local trancing object in
    there. When entering the block the constructor is called, the
    start message is printed. When leaving the block the descructor
    displays the finishing message.  For convenience there are some
    macros defining the local object so you can write:

    \code
    { // start a new block
      ACG_TRACE("What's the meaning of like?");

      for (years=0; years <= ACG::NumLimitsT<int>::max(); ++years)
      {
        // ... lengthy computation, show some progress
	ACG_TRACE_PROGRESS;
      }

      std::cout << "42";

    } // leaving the block will display "finished\n"
    \endcode

    If you want to trace only one command you can also use:
    \code
    ACG_TRACE_CMD("What's the meaning of like?", calc_meaning_of_life());
    \endcode

    \see TimedTracing
**/
class ACGDLLEXPORT Tracing
{
public:

  /// Constructor
  Tracing(const std::string& _text,
	  std::ostream&      _os = std::cerr)
    : os_(_os)
  { os_ << _text << "  " << std::flush; }


  /// Destructor
  ~Tracing() { os_ << "finished.\n" << std::flush; }


  /** Show some progress (the next frame of a rotating bar is
      displayed for every call of this function) */
  void progress() {
    os_ << progress_[((++idx_)&=0x3)] << "\b" << std::flush;
  }


protected:

  std::ostream&  os_;
  static char    progress_[4];
  static unsigned char idx_;
};


//== CLASS DEFINITION =========================================================



/** \class TimedTracing Tracing.hh <ACG/Utils/Tracing.hh>

    This class is similar to ACG::Tracing. It just starts a timer
    at the start of the process to be traced and outputs the
    time the computation took when leaving the block. Again
    we provide some convenience macros.

    \code
    { // start a new block
      ACG_TIMES_TRACE("What's the meaning of like?");

      for (years=0; years <= ACG::NumLimitsT<int>::max(); ++years)
      {
        // ... lengthy computation, show some progress
	ACG_TRACE_PROGRESS;
      }

      std::cout << "42";

    } // leaving the block will display " <...> secs, finished\n"
    \endcode

    If you want to trace only one command you can also use:
    \code
    ACG_TIMED_TRACE_CMD("What's the meaning of like?", calc_meaning_of_life());
    \endcode

    \see TimedTracing
**/
class ACGDLLEXPORT TimedTracing : public Tracing
{
public:

  /// default constructor
  TimedTracing( const std::string& _text,
		std::ostream&      _os = std::cerr )
    : Tracing(_text, _os)
  { timer_.start(); }


  /// destructor
  ~TimedTracing()
  {
    os_ << timer_.stop()*0.001 << " secs, ";
  }

private:

  StopWatch timer_;
};



//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_TRACING_HH defined
//=============================================================================

