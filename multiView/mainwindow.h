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

    // Reading
    void readJSON(const QJsonObject &json, Viewer *v);
    void openJSON(Viewer* v);

private Q_SLOTS:
    void openSkullMesh();
    void openFibulaMesh();
    void openMandJSON();
    void openFibJSON();

private:
    int sliderMax;
    int fibulaOffsetMax;
};

#endif // MAINWINDOW_H
