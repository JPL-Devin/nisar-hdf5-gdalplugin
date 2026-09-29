/******************************************************************************
 *
 * Project:  Hierarchical Data Format Release 5 (HDF5)
 * Authors:  Denis Nadeau <denis.nadeau@gmail.com>
 * Sam Gillingham <gillingham.sam@gmail.com>
 *
 ******************************************************************************
 * Copyright (c) 2008-2018, Even Rouault <even.rouault at spatialys.com>
 *
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

// This file contains the Virtual File Layer implementation that calls through
// to the VSI functions and should be included by HDF5 based drivers that wish
// to use the VFL for /vsi file system support.

// hdf5vfl.h
#ifndef HDF5VFL_H_INCLUDED_
#define HDF5VFL_H_INCLUDED_

#include "cpl_port.h"
#include <hdf5.h>
#include <cstdint>

// --------------------------------------------------------------------------
// CRITICAL: Namespace isolation for out-of-tree plugins.
// Prevents symbol collision with GDAL's internal HDF5/NetCDF drivers.
// --------------------------------------------------------------------------
namespace NisarVFL {
    // We only expose the two functions the rest of your plugin actually needs to see.
    // All the heavy lifting (structs, static functions) stays hidden inside hdf5vfl.cpp
    hid_t HDF5VFLGetFileDriver();
    void HDF5VFLUnloadFileDriver();

    // Process-wide counters of read callbacks served by this VFL (all files).
    void HDF5VFLGetReadStats(uint64_t *pnReadCalls, uint64_t *pnBytesRead);
}

#endif /* HDF5VFL_H_INCLUDED_ */
