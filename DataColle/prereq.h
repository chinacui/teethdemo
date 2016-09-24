#pragma once
#pragma warning(disable: 4251)
#ifdef DATACOLLE_EXPORTS
#define DATACOLLE_API extern "C" __declspec(dllexport)
#define DATACOLLE_CLASS __declspec(dllexport)
#else
#define DATACOLLE_API extern "C" __declspec(dllimport)
#define DATACOLLE_CLASS __declspec(dllimport)
#endif