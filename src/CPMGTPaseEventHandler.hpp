/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef CPMGTPASEEVENTHANDLER_HPP_
#define CPMGTPASEEVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"

/**
 * An event class that can be used to calculate the time taken to
 * execute various parts of a CPM+GTPase simulation.
 */
class CPMGTPaseEventHandler : public GenericEventHandler<11, CPMGTPaseEventHandler>
{
public:

    /** Character array holding cell_based event names. There are eleven cell_based events. */
    static const char* EventName[11];

    /** Definition of cell_based event types. */
    typedef enum
    {
        VOLUME=0,
        SURFACE,
        ADHESION,
        SETUP_MESH,
        UPDATE_SOLUTION,
        SOLVE_PDE,
        UPDATE_CELL_DATA,
        IO,
        RHO_CONTRACTION,
        BARBED_ENDS,
        EVERYTHING
    } CPMGTPaseEventType;
};

#endif /*CPMGTPASEEVENTHANDLER_HPP_*/
