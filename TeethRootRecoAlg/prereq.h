#pragma once
#pragma warning(disable: 4251)

#ifdef TEETHROOTRECOALG_EXPORTS
#define TEETHROOTRECOALG_API extern "C" __declspec(dllexport)
#define TEETHROOTRECOALG_CLASS __declspec(dllexport)
#else
#define TEETHROOTRECOALG_API extern "C" __declspec(dllimport)
#define TEETHROOTRECOALG_CLASS __declspec(dllimport)
#endif