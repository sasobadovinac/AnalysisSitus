//-----------------------------------------------------------------------------
// Created on: 17 December 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <asiVisu_InteractorStylePick2d.h>

// VTK includes
#include <vtkCallbackCommand.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>

// Qt includes
#pragma warning(push, 0)
#include <QPoint>
#pragma warning(pop)

vtkStandardNewMacro(asiVisu_InteractorStylePick2d);

//! Default constructor.
asiVisu_InteractorStylePick2d::asiVisu_InteractorStylePick2d()
: vtkInteractorStyleImage()
{
  m_pPickInput = new asiVisu_PickInput();
}

//! Default destructor.
asiVisu_InteractorStylePick2d::~asiVisu_InteractorStylePick2d() 
{
  delete m_pPickInput;
}

//! Prohibited copy constructor.
asiVisu_InteractorStylePick2d::asiVisu_InteractorStylePick2d(const asiVisu_InteractorStylePick2d&)
: vtkInteractorStyleImage(),
  m_pPickInput(NULL)
{
}

//! Prohibited assignment operator.
void asiVisu_InteractorStylePick2d::operator=(const asiVisu_InteractorStylePick2d&)
{
}

//! Sets Renderer instance.
//! \param theRenderer [in] Renderer instance to set.
void asiVisu_InteractorStylePick2d::SetRenderer(const vtkSmartPointer<vtkRenderer>& theRenderer)
{
  m_renderer = theRenderer;
}

//! Returns Renderer instance.
//! \return Renderer instance.
vtkSmartPointer<vtkRenderer> asiVisu_InteractorStylePick2d::GetRenderer() const
{
  return m_renderer;
}

//! Callback for "Mouse Move" event.
void asiVisu_InteractorStylePick2d::OnMouseMove()
{
  // Invoke basic method
  vtkInteractorStyleImage::OnMouseMove();

  const int pos[2] = { this->Interactor->GetEventPosition()[0],
                       this->Interactor->GetEventPosition()[1] };

  QPoint aPickPoint(pos[0], pos[1]);
  m_pPickInput->Start      = aPickPoint;
  m_pPickInput->Finish     = aPickPoint;
  m_pPickInput->IsMultiple = false;

  // Invoke observers
  this->InvokeEvent(EVENT_DETECT_DEFAULT, m_pPickInput);
}

//! Callback for "Left Button Down" event.
void asiVisu_InteractorStylePick2d::OnLeftButtonDown()
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0],
                          this->Interactor->GetEventPosition()[1]);

  if ( this->CurrentRenderer == NULL )
    return;

  const int pos[2] = { this->Interactor->GetEventPosition()[0],
                       this->Interactor->GetEventPosition()[1] };
  //
  QPoint aPickPoint(pos[0], pos[1]);
  m_pPickInput->Start  = aPickPoint;
  m_pPickInput->Finish = aPickPoint;

  if ( this->Interactor->GetShiftKey() )
    m_pPickInput->IsMultiple = true;
  else
    m_pPickInput->IsMultiple = false;

  // Invoke observers
  this->InvokeEvent(EVENT_PICK_DEFAULT, m_pPickInput);
}

//! Callback for "Key Press" event.
void asiVisu_InteractorStylePick2d::OnKeyPress()
{
  std::string key = this->Interactor->GetKeySym();
  if ( key == "Delete" )
    this->InvokeEvent(EVENT_DELETE);
  if ( key == "j" )
    this->InvokeEvent(EVENT_JOIN);
}