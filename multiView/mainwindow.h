#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "viewerfibula.h"

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
    void initDisplayDockWidgets();

    // File menu
    QActionGroup *fileActionGroup;
    void initFileMenu();
    void initToolBars();
    void initFileActions();

private Q_SLOTS:
    void openSkullMesh();
    void openFibulaMesh();

private:
    int sliderMax;
    int fibulaOffsetMax;
};

#endif // MAINWINDOW_H
