#include "mainwindow.h"
#include <QLayout>
#include <QGroupBox>
#include <QDockWidget>
#include <QSlider>
#include <QFormLayout>
#include <QPushButton>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    if (this->objectName().isEmpty())
        this->setObjectName("window");
    this->resize(929, 891);

    sliderMax = 100;
    fibulaOffsetMax = 100;

    // The qglviewer
    //QString skullFilename = "C:\\Users\\Medmax\\Documents\\Project\\Mand_B.off";
    StandardCamera *sc = new StandardCamera();
    skullViewer = new Viewer(this, sc, sliderMax);

    // The fibula viewer
    //QString fibulaFilename = "C:\\Users\\Medmax\\Documents\\Project\\Fibula_G.off";
    StandardCamera *scFibula = new StandardCamera();
    fibulaViewer = new ViewerFibula(this, scFibula, sliderMax, fibulaOffsetMax);

    // Main widget
    QWidget *mainWidget = new QWidget(this);
    //QWidget *fibulaWidget = new QWidget(this);

    // Horizontal layout
    QHBoxLayout *windowLayout = new QHBoxLayout();

    // Add the viewer to the layout
    windowLayout->addWidget(skullViewer);
    windowLayout->addWidget(fibulaViewer);

    // Add the layout to the main widget
    mainWidget->setLayout(windowLayout);

    QGroupBox * viewerGroupBox = new QGroupBox();

    QGridLayout * gridLayout = new QGridLayout(viewerGroupBox);
    gridLayout->setObjectName("gridLayout");

    gridLayout->addWidget(mainWidget, 0, 1, 1, 1);

    viewerGroupBox->setLayout(gridLayout);

    this->setCentralWidget(viewerGroupBox);

    initDisplayDockWidgets();
    initFileMenu();
    initToolBars();

    this->setWindowTitle("MedMax");
}

MainWindow::~MainWindow()
{

}

void MainWindow::initDisplayDockWidgets(){

    skullDockWidget = new QDockWidget("Plane controls");

    QHBoxLayout* layout = new QHBoxLayout();

    // The contents of the dockWidget
    QWidget *contentsMand = new QWidget();
    QFormLayout *contentLayoutMand = new QFormLayout();

    QWidget *contentsFibula = new QWidget();
    QFormLayout *contentLayoutFibula = new QFormLayout();

    // Add the sliders (skull)
    QSlider *leftPlaneSlider = new QSlider(Qt::Horizontal);
    leftPlaneSlider->setMaximum(sliderMax);
    contentLayoutMand->addRow("Left slider", leftPlaneSlider);

    QSlider *rightPlaneSlider = new QSlider(Qt::Horizontal);
    rightPlaneSlider->setMaximum(sliderMax);
    contentLayoutMand->addRow("Right slider", rightPlaneSlider);

    // PlaneRotation sliders (skull)
    QSlider *leftPlaneRotationSlider = new QSlider(Qt::Horizontal);
    leftPlaneRotationSlider->setMaximum(sliderMax);
    contentLayoutMand->addRow("Left rotation", leftPlaneRotationSlider);

    QSlider *rightPlaneRotationSlider = new QSlider(Qt::Horizontal);
    rightPlaneRotationSlider->setMaximum(sliderMax);
    contentLayoutMand->addRow("Right rotation", rightPlaneRotationSlider);

    // Connect the skull sliders
    connect(leftPlaneSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), skullViewer, &Viewer::moveLeftPlane);
    connect(rightPlaneSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), skullViewer, &Viewer::moveRightPlane);
    connect(leftPlaneRotationSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), skullViewer, &Viewer::rotateLeftPlane);
    connect(rightPlaneRotationSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), skullViewer, &Viewer::rotateRightPlane);
    connect(leftPlaneRotationSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), fibulaViewer, &ViewerFibula::rotateLeftPlane);
    connect(rightPlaneRotationSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), fibulaViewer, &ViewerFibula::rotateRightPlane);

    // Allow the viewers to change the values
    connect(skullViewer, &Viewer::setLRSliderValue, leftPlaneRotationSlider, &QSlider::setValue);
    connect(skullViewer, &Viewer::setRRSliderValue, rightPlaneRotationSlider, &QSlider::setValue);
    connect(leftPlaneSlider, &QSlider::sliderReleased, skullViewer, &Viewer::onLeftSliderReleased);
    connect(rightPlaneSlider, &QSlider::sliderReleased, skullViewer, &Viewer::onRightSliderReleased);
    connect(skullViewer, &Viewer::setLMSliderValue, leftPlaneSlider, &QSlider::setValue);
    connect(skullViewer, &Viewer::setRMSliderValue, rightPlaneSlider, &QSlider::setValue);

    // Add the slider (fibula)
    QSlider *fibulaSlider = new QSlider(Qt::Horizontal);
    contentLayoutFibula->addRow("Fibula position slider", fibulaSlider);
    fibulaSlider->setMaximum(fibulaOffsetMax);

    // Connect the fibula slider
    connect(fibulaSlider, static_cast<void (QSlider::*)(int)>(&QSlider::sliderMoved), fibulaViewer, &ViewerFibula::movePlanes);
    connect(fibulaSlider, &QSlider::sliderReleased, fibulaViewer, &ViewerFibula::planesMoved);
    connect(fibulaViewer, &ViewerFibula::setPlaneSliderValue, fibulaSlider, &QSlider::setValue);

    // Connect the two views
    connect(skullViewer, &Viewer::leftPosChanged, fibulaViewer, &ViewerFibula::movePlaneDistance);
    connect(skullViewer, &Viewer::rightPosChanged, fibulaViewer, &ViewerFibula::movePlaneDistance);

    connect(skullViewer, &Viewer::ghostPlanesAdded, fibulaViewer, &ViewerFibula::ghostPlanesRecieved);
    connect(skullViewer, &Viewer::ghostPlanesTranslated, fibulaViewer, &ViewerFibula::middlePlaneMoved);

    connect(fibulaViewer, &ViewerFibula::sendToManible, skullViewer, &Viewer::recieveFromFibulaMesh);

    /*connect(skullViewer, &Viewer::haltMeshUpdate, fibulaViewer, &ViewerFibula::toHaltMeshUpdate);
    connect(skullViewer, &Viewer::continueMeshUpdate, fibulaViewer, &ViewerFibula::toContinueMeshUpdate);*/

    connect(skullViewer, &Viewer::noGhostPlanesToSend, fibulaViewer, &ViewerFibula::noGhostPlanesToRecieve);
    connect(skullViewer, &Viewer::preparingToCut, fibulaViewer, &ViewerFibula::uncutMesh);

    contentsMand->setLayout(contentLayoutMand);
    contentsFibula->setLayout(contentLayoutFibula);

    layout->addWidget(contentsMand);
    layout->addWidget(contentsFibula);

    QWidget* controlWidget = new QWidget();
    controlWidget->setLayout(layout);

    skullDockWidget->setWidget(controlWidget);

    this->addDockWidget(Qt::BottomDockWidgetArea, skullDockWidget);
}

void MainWindow::initFileActions(){
    fileActionGroup = new QActionGroup(this);

    QAction *openFileSkullAction = new QAction("Open skull mesh", this);
    connect(openFileSkullAction, &QAction::triggered, this, &MainWindow::openSkullMesh);

    QAction *openFileFibulaAction = new QAction("Open fibula mesh", this);
    connect(openFileFibulaAction, &QAction::triggered, this, &MainWindow::openFibulaMesh);

    QAction *cutMeshAction = new QAction("Cut", this);
    connect(cutMeshAction, &QAction::triggered, skullViewer, &Viewer::cutMesh);
    connect(skullViewer, &Viewer::okToCut, fibulaViewer, &ViewerFibula::cutMesh);

    QAction *unCutMeshAction = new QAction("Undo cut", this);
    connect(unCutMeshAction, &QAction::triggered, skullViewer, &Viewer::uncutMesh);
    connect(unCutMeshAction, &QAction::triggered, fibulaViewer, &ViewerFibula::uncutMesh);

    QAction *drawMeshAction = new QAction("Draw Mesh", this);
    connect(drawMeshAction, &QAction::triggered, skullViewer, &Viewer::drawMesh);
    connect(drawMeshAction, &QAction::triggered, fibulaViewer, &ViewerFibula::drawMesh);

    fileActionGroup->addAction(openFileSkullAction);
    fileActionGroup->addAction(openFileFibulaAction);
    fileActionGroup->addAction(unCutMeshAction);
    fileActionGroup->addAction(cutMeshAction);
    fileActionGroup->addAction(drawMeshAction);

}

void MainWindow::initFileMenu(){
    initFileActions();

    QMenu *fileMenu = menuBar()->addMenu(tr("File"));
    fileMenu->addActions(fileActionGroup->actions());
}

void MainWindow::initToolBars () {
    QToolBar *fileToolBar = new QToolBar(this);
    fileToolBar->addActions(fileActionGroup->actions());
    addToolBar(fileToolBar);
}

void MainWindow::openSkullMesh(){
    QString openFileNameLabel, selectedFilter;

    QString fileFilter = "OFF (*.off)";

    QString fileName = QFileDialog::getOpenFileName(this, tr("Select a skull mesh"), openFileNameLabel, fileFilter, &selectedFilter);

    if(fileName.isEmpty()) return;

    skullViewer->openOFF(fileName);
}

void MainWindow::openFibulaMesh(){
    QString openFileNameLabel, selectedFilter;

    QString fileFilter = "OFF (*.off)";

    QString fileName = QFileDialog::getOpenFileName(this, tr("Select a fibula mesh"), openFileNameLabel, fileFilter, &selectedFilter);

    if(fileName.isEmpty()) return;

    fibulaViewer->openOFF(fileName);
}

