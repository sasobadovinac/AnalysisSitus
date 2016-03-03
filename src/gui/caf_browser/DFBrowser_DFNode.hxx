//-----------------------------------------------------------------------------
// Created on: 18 February 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

#ifndef _DFBrowser_DFNode_HeaderFile
#define _DFBrowser_DFNode_HeaderFile

#include <Handle_DFBrowser_DFNode.hxx>
#include <MMgt_TShared.hxx>

#include <Standard_CString.hxx>
#include <TCollection_AsciiString.hxx>
#include <DFBrowser_Colors.hxx>
#include <DFBrowser_Picture.hxx>
#include <DFBrowser_NodeType.hxx>
#include <Handle_DFBrowser_DFTree.hxx>

class DFBrowser_DFNode : public MMgt_TShared
{
 public:

  Standard_EXPORT DFBrowser_DFNode();

  Standard_EXPORT virtual DFBrowser_NodeType GetType() const = 0;

  Standard_EXPORT virtual void Update() = 0;

  Standard_EXPORT void Next(const Handle(DFBrowser_DFNode)& theNext);

  inline const Handle_DFBrowser_DFNode & Next() const
  {
    return myNext;
  }

  Standard_EXPORT void Parent(const Handle(DFBrowser_DFNode)& theParent);

  inline const Handle_DFBrowser_DFNode & Parent() const
  {
    return myParent;
  }

  Standard_EXPORT virtual void AddSub(Handle(DFBrowser_DFNode)& theNode) = 0;

  Standard_EXPORT virtual Handle_DFBrowser_DFNode Sub() const = 0;

  Standard_EXPORT virtual const TCollection_AsciiString & Name() = 0;

  inline const Handle_DFBrowser_DFTree & Tree() const
  {
    return myTree;
  }

  inline void Tree(const Handle(DFBrowser_DFTree)& theTree)
  {
    myTree = theTree;
  }

  inline Standard_Boolean Opened() const
  {
    return myIsOpened;
  }

  Standard_EXPORT void Opened(const Standard_Boolean theOpened);

  inline Standard_Boolean CanOpen() const
  {
    return myCanOpen;
  }

  Standard_EXPORT void CanOpen(const Standard_Boolean theCanOpen);

  inline Standard_Boolean Selected() const
  {
    return myIsSelected;
  }

  inline void Selected(const Standard_Boolean theIsSelected)
  {
    myIsSelected = theIsSelected;
  }

  inline Standard_Boolean Changed() const
  {
    return myIsChanged;
  }

  inline void Changed(const Standard_Boolean theIsChanged)
  {
    myIsChanged = theIsChanged;
  }

  inline Standard_Boolean Visible() const
  {
    return myIsVisible;
  }

  Standard_EXPORT void Visible(const Standard_Boolean theIsVisible);

  inline DFBrowser_Colors Color() const
  {
    return myColor;
  }

  inline DFBrowser_Picture Pixmap() const
  {
    return myPixmap;
  }

  Standard_EXPORT virtual void Del() = 0;

  DEFINE_STANDARD_RTTI_INLINE(DFBrowser_DFNode, MMgt_TShared)

 protected:

  TCollection_AsciiString myName;
  DFBrowser_Colors myColor;
  DFBrowser_Picture myPixmap;

 private:

  Handle_DFBrowser_DFTree myTree;
  Handle_DFBrowser_DFNode myNext;
  Handle_DFBrowser_DFNode myParent;
  Standard_Boolean myIsOpened;
  Standard_Boolean myCanOpen;
  Standard_Boolean myIsSelected;
  Standard_Boolean myIsChanged;
  Standard_Boolean myIsVisible;
};

#endif