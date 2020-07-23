#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QActionGroup>
#include "Display/viewerfibula.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    // Main viewers
    Viewer *skullViewer;
    ViewerFibula *fibulaViewer;
    QDockWidget *skullDockWidget;
    QDockWidget *editMenuWidget;
    QDockWidget *editFragmentDockWidget;
    void initDisplayDockWidgets();

    // File menu
    QActionGroup *fileActionGroup;
    void initFileMenu();
    void initToolBars();
    void initFileActions();
    void initEditMenu();
    void initEditFragmentsMenu();

    void readJSON(const QJsonObject &json, Viewer *v);
    bool openJSON(Viewer* v);
    void filesOpened();

private Q_SLOTS:
    void openMandJSON();
    void openFibJSON();
    void setFragRadios();
    void enableFragmentEditing();
    void disableFragmentEditing();
    void displayEditMenu();
    void displayEditFragmentMenu();
    void displayFragmentMenuButton();
    void hideFragmentMenuButton();
    void editPlane(unsigned int);
    void editBoxCentre(unsigned int);
    void editBoxStart(unsigned int);
    void editBoxEnd(unsigned int);
    void toEditPlane(bool);
    void toEditBoxCentre(bool);
    void toEditBoxStart(bool);
    void toEditBoxEnd(bool);

private:
    int sliderMax = 100;
    int fibulaOffsetMax;
    QRadioButton *radioFrag1, *radioFrag2, *radioFrag3, *radioFragPlanes;
    QGroupBox *groupRadioBox;
    QDockWidget *loadedMeshes;
    QPushButton *editMenuButton, *editFragmentMenuButton, *toggleDrawMeshButton;
    bool isOpenMand = false;
    bool isOpenFib = false;
    unsigned int currentBox, currentPlane;
};

#endif // MAINWINDOW_H
