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
//  CLASS StopWatch
//
//=============================================================================


#ifndef ACG_STOPWATCH_HH
#define ACG_STOPWATCH_HH


//== INCLUDES =================================================================


#ifdef _WIN32

#include <windows.h>

#else // Linux

#include <sys/time.h>

#endif

#include "../Config/ACGDefines.hh"


//== NAMESPACES ===============================================================


namespace ACG {


//== CLASS DEFINITION =========================================================



/** \class StopWatch StopWatch.hh <ACG/Utils/StopWatch.hh>

    This class can be used for measuring time, providing milli-second
    precision by using the gettimeofday() funtion.  It is e.g. used in
    the class ACG::TimedTracing.
**/

class StopWatch
{
public:

  /// Constructor
  StopWatch() {
    #ifdef _WIN32
      QueryPerformanceFrequency(&freq_);
    #else
      starttime_.tv_sec  = 0;
      starttime_.tv_usec = 0;
      endtime_.tv_sec    = 0;
      endtime_.tv_usec   = 0;
    #endif
  }

  /// Destructor
  ~StopWatch() {}

  /// Start time measurement
  void start() {
    #ifdef _WIN32
      QueryPerformanceCounter(&starttime_);
    #else
      starttime_ = current_time();
    #endif
  }

  /// Restart, return time elapsed until now.
  double restart() {
    #ifdef _WIN32
      QueryPerformanceCounter(&endtime_);
    #else
      endtime_ = current_time();
    #endif

    double t = elapsed();
    start();
    return t;
  }

  /// Stop time measurement, return time.
  double stop() {
    #ifdef _WIN32
      QueryPerformanceCounter(&endtime_);
    #else
      endtime_ = current_time();
    #endif

    return elapsed();
  }

  /// Get the total time in milli-seconds (watch has to be stopped).
  double elapsed() const {
    #ifdef _WIN32
      return (double)(endtime_.QuadPart - starttime_.QuadPart)
	   / (double)freq_.QuadPart * 1000.0f;
    #else
      return ((endtime_.tv_sec  - starttime_.tv_sec )*1000.0 +
	      (endtime_.tv_usec - starttime_.tv_usec)*0.001);
    #endif
  }


private:

  #ifdef _WIN32
    LARGE_INTEGER starttime_, endtime_;
    LARGE_INTEGER freq_;
  #else // Linux
    timeval current_time() const {
      struct timeval tv;
      gettimeofday(&tv, 0);
      return tv;
    }

    timeval starttime_, endtime_;
  #endif

};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_STOPWATCH_HH defined
//=============================================================================

