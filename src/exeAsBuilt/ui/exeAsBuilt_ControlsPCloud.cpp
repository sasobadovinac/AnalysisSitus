//-----------------------------------------------------------------------------
// Created on: 21 June 2016
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/
//-----------------------------------------------------------------------------

// Own include
#include <exeAsBuilt_ControlsPCloud.h>

// exeAsBuilt includes
#include <exeAsBuilt_CommonFacilities.h>
#include <exeAsBuilt_DialogCloudify.h>

// asiAlgo includes
#include <asiAlgo_Timer.h>

// asiUI includes
#include <asiUI_Common.h>

// Qt include
#include <QGroupBox>

//-----------------------------------------------------------------------------

#define BTN_MIN_WIDTH 120

//-----------------------------------------------------------------------------
// Construction & destruction
//-----------------------------------------------------------------------------

//! Constructor.
//! \param parent [in] parent widget.
exeAsBuilt_ControlsPCloud::exeAsBuilt_ControlsPCloud(QWidget* parent) : QWidget(parent)
{
  // Main layout
  m_pMainLayout = new QVBoxLayout();

  // Buttons
  m_widgets.pLoadPoints   = new QPushButton("Load points");
  m_widgets.pEstimNormals = new QPushButton("Estimate normals");
  //
  m_widgets.pLoadPoints   -> setMinimumWidth(BTN_MIN_WIDTH);
  m_widgets.pEstimNormals -> setMinimumWidth(BTN_MIN_WIDTH);

  // Group for data exchange
  QGroupBox*   pDEGroup = new QGroupBox("Data exchange");
  QVBoxLayout* pDELay   = new QVBoxLayout(pDEGroup);
  //
  pDELay->addWidget(m_widgets.pLoadPoints);

  // Group for analysis
  QGroupBox*   pAnalysisGroup = new QGroupBox("Analysis");
  QVBoxLayout* pAnalysisLay   = new QVBoxLayout(pAnalysisGroup);
  //
  pAnalysisLay->addWidget(m_widgets.pEstimNormals);

  // Set layout
  m_pMainLayout->addWidget(pDEGroup);
  m_pMainLayout->addWidget(pAnalysisGroup);
  //
  m_pMainLayout->setAlignment(Qt::AlignTop);
  //
  this->setLayout(m_pMainLayout);

  // Connect signals to slots
  connect( m_widgets.pLoadPoints,   SIGNAL( clicked() ), SLOT( onLoadPoints   () ) );
  connect( m_widgets.pEstimNormals, SIGNAL( clicked() ), SLOT( onEstimNormals () ) );
}

//! Destructor.
exeAsBuilt_ControlsPCloud::~exeAsBuilt_ControlsPCloud()
{
  delete m_pMainLayout;
  m_widgets.Release();
}

//-----------------------------------------------------------------------------

//! Reaction on user request to load points.
void exeAsBuilt_ControlsPCloud::onLoadPoints()
{
  TIMER_NEW
  TIMER_GO

  // Select file for point cloud data
  QString
    qfilename = asiUI_Common::selectFile(QStringList() << "*.pts" << "*.xyz",
                                        "Choose point cloud to load", "",
                                         asiUI_Common::OpenSaveAction_Open);
  //
  TCollection_AsciiString filename = QStr2AsciiStr(qfilename);

  // Load point cloud
  Handle(asiAlgo_PointCloud)
    cloud = new asiAlgo_PointCloud(exeAsBuilt_CommonFacilities::Instance()->Notifier,
                                   exeAsBuilt_CommonFacilities::Instance()->Plotter);
  //
  if ( !cloud->Load(filename) )
  {
    std::cout << "Error: cannot load point cloud" << std::endl;
    return;
  }
  std::cout << "Point cloud was loaded successfully" << std::endl;

  TIMER_FINISH
  TIMER_COUT_RESULT_MSG("Load point cloud")

  ActAPI_PlotterEntry IV(exeAsBuilt_CommonFacilities::Instance()->Plotter);
  IV.DRAW_POINTS( cloud->GetPoints(), Color_White, "Scan data" );
}

//-----------------------------------------------------------------------------

//! Reaction on estimation of normal field.
void exeAsBuilt_ControlsPCloud::onEstimNormals()
{
  // TODO: the following is just weird
  Handle(ActAPI_INode)
    N = exeAsBuilt_CommonFacilities::Instance()->Model->GetIVPointSetPartition()->GetNode(1);
  //
  Handle(asiData_IVPointSetNode) PN = Handle(asiData_IVPointSetNode)::DownCast(N);
  //
  if ( N.IsNull() || !PN->IsWellFormed() )
  {
    std::cout << "Error: point cloud is not ready" << std::endl;
    return;
  }

  // Calculate normals
  //asiAlgo_PointCloudNormals::Calculate( PN->GetPoints() );
  // TODO
}

//-----------------------------------------------------------------------------

//! Builds point cloud from the model.
void exeAsBuilt_ControlsPCloud::onCloudify()
{
  Handle(exeAsBuilt_CommonFacilities) cf = exeAsBuilt_CommonFacilities::Instance();
  Handle(asiData_PartNode)            part_n;
  TopoDS_Shape                        part;
  //
  if ( !asiUI_Common::PartShape(cf->Model, part_n, part) ) return;

  // Run dialog for cloudification
  exeAsBuilt_DialogCloudify*
    wCloudify = new exeAsBuilt_DialogCloudify(cf->Model, cf->Notifier, cf->Plotter, this);
  //
  wCloudify->show();
}