//-----------------------------------------------------------------------------
// Created on: 25 February 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

#ifndef visu_topo_graph_item_h
#define visu_topo_graph_item_h

// Visualization includes
#include <visu_common.h>

// VTK includes
#include <vtkAbstractArray.h>
#include <vtkColor.h>
#include <vtkContext2D.h>
#include <vtkContextMouseEvent.h>
#include <vtkContextScene.h>
#include <vtkDataSetAttributes.h>
#include <vtkGraph.h>
#include <vtkGraphItem.h>
#include <vtkStdString.h>
#include <vtkTextProperty.h>

// Qt includes
#include <QObject>

#define ARRNAME_LABELS          "Labels"
#define ARRNAME_COLORS          "Colors"
#define ARRNAME_GROUP           "Group"
#define ARRNAME_GROUP_ORDINARY  "Ordinary"
#define ARRNAME_GROUP_ADJACENT  "Adjacent"
#define ARRNAME_GROUP_COMPOUND  "Compound"

//! Item of topology graph.
class visu_topo_graph_item : public QObject,
                             public vtkGraphItem
{
  Q_OBJECT

public:

  static visu_topo_graph_item* New();
  vtkTypeMacro(visu_topo_graph_item, vtkGraphItem);

  virtual ~visu_topo_graph_item();

signals:

  void vertexPicked(const vtkIdType);

protected:

  vtkIdType focusedVertex;

  //---------------------------------------------------------------------------
  visu_topo_graph_item()
  {
    focusedVertex = -1;
  }

  //---------------------------------------------------------------------------
  virtual vtkColor4ub VertexColor(vtkIdType vertex)
  {
    if ( vertex == focusedVertex )
      return vtkColor4ub(255, 0, 0, 255);

    vtkAbstractArray* domain = this->GetGraph()->GetVertexData()->GetAbstractArray(ARRNAME_GROUP);
    //
    if ( domain && domain->GetVariantValue(vertex).ToString() == ARRNAME_GROUP_ADJACENT )
      return vtkColor4ub(128, 175, 128, 255);
    if ( domain && domain->GetVariantValue(vertex).ToString() == ARRNAME_GROUP_COMPOUND )
      return vtkColor4ub(215, 75, 75, 255);

    return vtkColor4ub(175, 128, 128, 255);
  }

  //---------------------------------------------------------------------------
  virtual vtkStdString VertexTooltip(vtkIdType vertex)
  {
    vtkAbstractArray* labels = this->GetGraph()->GetVertexData()->GetAbstractArray(ARRNAME_LABELS);
    if ( labels )
      return labels->GetVariantValue(vertex).ToString();

    return "";
  }

  //---------------------------------------------------------------------------
  virtual vtkColor4ub EdgeColor(vtkIdType /*line*/, vtkIdType /*point*/)
  {
    return vtkColor4ub(128, 128, 128, 128);
  }

  //---------------------------------------------------------------------------
  virtual float EdgeWidth(vtkIdType /*line*/, vtkIdType /*point*/)
  {
    return 1.0f;
  }

  //---------------------------------------------------------------------------
  virtual void PaintBuffers(vtkContext2D* painter)
  {
    // Turn off the tooltip if the superclass turned it on
    this->PlaceTooltip(-1);

    this->Superclass::PaintBuffers(painter);

    if ( focusedVertex >= 0 )
    {
      painter->GetTextProp()->SetColor(1, 1, 1);
      painter->GetTextProp()->SetJustificationToCentered();
      painter->GetTextProp()->BoldOff();

      vtkVector2f pos = this->VertexPosition(focusedVertex);
      vtkStdString label = this->VertexTooltip(focusedVertex);

      painter->GetTextProp()->SetFontSize(12);
      painter->DrawString(pos.GetX(), pos.GetY(), label);
    }
  }

  //---------------------------------------------------------------------------
  virtual bool MouseButtonPressEvent(const vtkContextMouseEvent& evt)
  {
    this->Superclass::MouseButtonPressEvent(evt);

    focusedVertex = this->HitVertex( evt.GetPos() );

    this->GetGraph()->Modified();
    this->GetScene()->SetDirty(true);

    emit vertexPicked(focusedVertex);
    return true;
  }

};

#endif