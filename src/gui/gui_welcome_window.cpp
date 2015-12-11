//-----------------------------------------------------------------------------
// Created on: 08 December 2015
// Created by: Sergey SLYADNEV
//-----------------------------------------------------------------------------
// Web: http://quaoar.su/blog/
//-----------------------------------------------------------------------------

// Own include
#include <gui_welcome_window.h>

// A-Situs (GUI) includes
#include <gui_main_window_asitus.h>
#include <gui_main_window_skinning.h>

// Qt includes
#pragma warning(push, 0)
#include <QDesktopWidget>
#include <QVBoxLayout>
#pragma warning(pop)

//-----------------------------------------------------------------------------
// Construction & destruction
//-----------------------------------------------------------------------------

//! Constructor.
gui_welcome_window::gui_welcome_window() : QMainWindow()
{
  this->createControls();

  this->setCentralWidget(m_widgets.pMainPane);
  this->setWindowTitle("Analysis Situs: choose composition");
}

//! Destructor.
gui_welcome_window::~gui_welcome_window()
{}

//-----------------------------------------------------------------------------

//! Creates all widgets.
void gui_welcome_window::createControls()
{
  // Buttons
  m_widgets.pAnalysis = new QPushButton("Analysis");
  m_widgets.pSkinning = new QPushButton("Skinning");
  m_widgets.pAnalysis->setMaximumWidth(100);
  m_widgets.pSkinning->setMaximumWidth(100);
  //
  // Layout for buttons
  QVBoxLayout* pMainLayout = new QVBoxLayout();
  pMainLayout->addWidget(m_widgets.pAnalysis);
  pMainLayout->addWidget(m_widgets.pSkinning);
  //
  pMainLayout->setAlignment(Qt::AlignTop | Qt::AlignHCenter);
  pMainLayout->setContentsMargins(10, 10, 10, 10);

  // Sizing
  const int width  = 100;
  const int height = 80;
  //
  m_widgets.pMainPane = new QWidget();
  m_widgets.pMainPane->setMinimumSize(width, height);
  m_widgets.pMainPane->setLayout(pMainLayout);

  // Connect signals to slots
  connect( m_widgets.pAnalysis, SIGNAL( clicked() ), SLOT( onAnalysis() ) );
  connect( m_widgets.pSkinning, SIGNAL( clicked() ), SLOT( onSkinning() ) );
}

//-----------------------------------------------------------------------------

//! Reaction to clicking "Analysis" button.
void gui_welcome_window::onAnalysis()
{
  this->hide();
  //
  gui_main_window_asitus* pMainWindow = new gui_main_window_asitus();
  pMainWindow->show();
}

//! Reaction to clicking "Skinning" button.
void gui_welcome_window::onSkinning()
{
  this->hide();
  //
  gui_main_window_skinning* pMainWindow = new gui_main_window_skinning();
  pMainWindow->show();
}