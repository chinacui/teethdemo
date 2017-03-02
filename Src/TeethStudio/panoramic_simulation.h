#ifndef PANORAMIC_SIMULATION_H
#define PANORAMIC_SIMULATION_H

#include <QWidget>
#include "ui_panoramic_simulation.h"

class PanoramicSimulation : public QWidget
{
	Q_OBJECT

public:
	PanoramicSimulation(QWidget *parent = 0);
	~PanoramicSimulation();

private:
	Ui::PanoramicSimulation ui;
};

#endif // PANORAMIC_SIMULATION_H
