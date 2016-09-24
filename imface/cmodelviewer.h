#ifndef MODEL_VIEWER_H
#define MODEL_VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <vector>


class CModelViewer : public QGLViewer
{
protected:
	virtual void draw();
	virtual void init();
public:
	CModelViewer(QWidget *parent);
	~CModelViewer();


};

#endif // MODEL_VIEWER_H
