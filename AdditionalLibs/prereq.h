#pragma once
#pragma warning(disable: 4251)

#ifdef ADDITIONALLIBS_EXPORTS
#define ADDITIONALLIBS_API extern "C" __declspec(dllexport)
#define ADDITIONALLIBS_CLASS __declspec(dllexport)
#else
#define ADDITIONALLIBS_API extern "C" __declspec(dllimport)
#define ADDITIONALLIBS_CLASS __declspec(dllimport)
#endif