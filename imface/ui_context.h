#ifndef UI_CONTEXT_H
#define UI_CONTEXT_H
class CScene;
class CImfaceWindow;
class CUIContext
{
protected:
	static CScene* scene_;
	static CImfaceWindow* main_window_;
public:
	static CScene* GetScene();
	static void Init();

};
#endif