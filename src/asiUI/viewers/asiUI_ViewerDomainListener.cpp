//-----------------------------------------------------------------------------
// Created on: 07 November 2016 (99 years of October Revolution)
//-----------------------------------------------------------------------------
// Copyright (c) 2017, Sergey Slyadnev
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of Sergey Slyadnev nor the
//      names of all contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//-----------------------------------------------------------------------------

// Own include
#include <asiUI_ViewerDomainListener.h>

//-----------------------------------------------------------------------------

//! Constructor accepting all necessary facilities.
//! \param wViewerPart    [in] part viewer.
//! \param wViewerDomain  [in] domain viewer.
//! \param wViewerSurface [in] host geometry viewer.
//! \param model          [in] Data Model instance.
asiUI_ViewerDomainListener::asiUI_ViewerDomainListener(asiUI_ViewerPart*              wViewerPart,
                                                       asiUI_ViewerDomain*            wViewerDomain,
                                                       asiUI_ViewerHost*              wViewerSurface,
                                                       const Handle(asiEngine_Model)& model)
: QObject          (),
  m_wViewerPart    (wViewerPart),
  m_wViewerDomain  (wViewerDomain),
  m_wViewerSurface (wViewerSurface),
  m_model          (model)
{}

//-----------------------------------------------------------------------------

//! Destructor.
asiUI_ViewerDomainListener::~asiUI_ViewerDomainListener()
{}

//-----------------------------------------------------------------------------

//! Connects this listener to the target widget.
void asiUI_ViewerDomainListener::Connect()
{
  // Nothing yet
}
