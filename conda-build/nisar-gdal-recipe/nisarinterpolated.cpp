#include <cmath>
#include <cctype>
#include <algorithm>
#include <limits>

#include "nisarinterpolated.h"
#include "nisarinterpolatedrasterband.h"
#include "nisardataset.h"
#include "gdalwarper.h"
#include "vrtdataset.h"

// Split "NISAR:<file>[:<hdf5-path>]" into its parts (same rules as NisarDataset::Open):
// the last colon separates file and path unless it is followed by "//" (URL scheme).
static bool NisarSplitConnectionString(const char* pszInput, std::string& sFile, std::string& sPath)
{
    std::string s(pszInput);
    if (STARTS_WITH_CI(s.c_str(), "NISAR:")) s = s.substr(6);
    if (s.empty()) return false;

    sPath.clear();
    const size_t nColon = s.rfind(':');
    if (nColon != std::string::npos && s.compare(nColon, 3, "://") != 0) {
        sPath = s.substr(nColon + 1);
        sFile = s.substr(0, nColon);
    } else {
        sFile = s;
    }
    if (sFile.size() >= 2 && sFile.front() == '"' && sFile.back() == '"') {
        sFile = sFile.substr(1, sFile.size() - 2);
    }
    return !sFile.empty();
}

static std::string NisarUpper(const char* psz)
{
    std::string s(psz);
    for (auto& c : s) c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
    return s;
}

// Open a NISAR dataset with the driver forced; returns the concrete NisarDataset.
static NisarDataset* NisarOpenInternal(const std::string& sConn, char** papszOptions)
{
    const char* const apszAllowedDrivers[] = { "NISAR", nullptr };
    GDALDataset* poDS = static_cast<GDALDataset*>(GDALOpenEx(
        sConn.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY | GDAL_OF_INTERNAL,
        apszAllowedDrivers, papszOptions, nullptr));
    if (poDS == nullptr) return nullptr;
    NisarDataset* poNisar = dynamic_cast<NisarDataset*>(poDS);
    if (poNisar == nullptr) GDALClose(poDS);
    return poNisar;
}


// ====================================================================
// NisarInterpolatedDataset Implementation
// ====================================================================

NisarInterpolatedDataset::NisarInterpolatedDataset()
{
}

NisarInterpolatedDataset::~NisarInterpolatedDataset()
{
    // Clean up the DEM datasets when this dataset is closed
    if (m_poAlignedDEM) GDALClose(m_poAlignedDEM);
    if (m_poRawDEM) GDALClose(m_poRawDEM);
}

#if GDAL_VERSION_MAJOR < 3 || (GDAL_VERSION_MAJOR == 3 && GDAL_VERSION_MINOR < 12)
CPLErr NisarInterpolatedDataset::GetGeoTransform(double* padfTransform)
{
    memcpy(padfTransform, m_adfTargetGeoTransform, 6 * sizeof(double));
    return CE_None;
}
#else
CPLErr NisarInterpolatedDataset::GetGeoTransform(GDALGeoTransform& gt) const
{
    for (int i = 0; i < 6; ++i) {
        gt[i] = m_adfTargetGeoTransform[i];
    }
    return CE_None;
}
#endif

const OGRSpatialReference* NisarInterpolatedDataset::GetSpatialRef() const
{
    return m_oSRS.IsEmpty() ? nullptr : &m_oSRS;
}
GDALDataset* NisarInterpolatedDataset::Open(GDALOpenInfo* poOpenInfo)
{
    char** papszOO = poOpenInfo->papszOpenOptions;
    const char* pszDemFile = CSLFetchNameValue(papszOO, "DEM_FILE");
    const char* pszQuantity = CSLFetchNameValue(papszOO, "QUANTITY");
    const char* pszInstOpt = CSLFetchNameValue(papszOO, "INST");
    const char* pszFreqOpt = CSLFetchNameValue(papszOO, "FREQ");
    const char* pszPolOpt = CSLFetchNameValue(papszOO, "POL");

    std::string sFile, sCubePath;
    if (!NisarSplitConnectionString(poOpenInfo->pszFilename, sFile, sCubePath)) {
        CPLError(CE_Failure, CPLE_OpenFailed, "Interpolation: empty filename.");
        return nullptr;
    }
    if (sFile.size() < 3 || !EQUAL(sFile.c_str() + sFile.size() - 3, ".h5")) {
        return nullptr;  // not a NISAR granule
    }
    const std::string sQuotedFile = "NISAR:\"" + sFile + "\"";

    // 1. Identify the product; without an explicit cube path, resolve it from QUANTITY.
    std::string sInst, sProduct;
    bool bIsLevel1 = false;
    NisarDataset* poCoarseCubeDS = nullptr;

    if (sCubePath.empty()) {
        NisarDataset* poContainer = NisarOpenInternal(sQuotedFile, nullptr);
        if (poContainer == nullptr) {
            CPLError(CE_Failure, CPLE_OpenFailed, "Interpolation: failed to open %s", sFile.c_str());
            return nullptr;
        }
        sInst = poContainer->GetInstrument();
        sProduct = poContainer->GetProductType();
        bIsLevel1 = poContainer->IsLevel1();
        GDALClose(poContainer);

        if (pszInstOpt) sInst = NisarUpper(pszInstOpt);
        if (sInst.empty() || sProduct.empty()) {
            CPLError(CE_Failure, CPLE_OpenFailed,
                     "Interpolation: could not identify NISAR product in %s", sFile.c_str());
            return nullptr;
        }

        // L1 quantities live in metadata/geolocationGrid, L2/L3 in metadata/radarGrid.
        const char* pszCubeGroup = bIsLevel1 ? "geolocationGrid" : "radarGrid";
        sCubePath = "/science/" + sInst + "/" + sProduct + "/metadata/" + pszCubeGroup + "/" + pszQuantity;
        CPLDebug("NISAR_DRIVER", "Interpolation: QUANTITY=%s resolved to %s", pszQuantity, sCubePath.c_str());
    }

    CPLDebug("NISAR_DRIVER", "Interpolation: Opening coarse cube at %s", sCubePath.c_str());
    poCoarseCubeDS = NisarOpenInternal(sQuotedFile + ":" + sCubePath, nullptr);
    if (poCoarseCubeDS == nullptr || poCoarseCubeDS->GetRasterCount() < 2) {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to open valid 3D coarse metadata cube at %s", sCubePath.c_str());
        if (poCoarseCubeDS) GDALClose(poCoarseCubeDS);
        return nullptr;
    }

    if (sProduct.empty()) {
        sInst = pszInstOpt ? NisarUpper(pszInstOpt) : poCoarseCubeDS->GetInstrument();
        sProduct = poCoarseCubeDS->GetProductType();
        bIsLevel1 = poCoarseCubeDS->IsLevel1();
    }

    // 2. Open the reference (imaging) grid whose shape/georeferencing the output inherits.
    const std::string sFreq = pszFreqOpt ? NisarUpper(pszFreqOpt) : std::string("A");
    std::string sPol = pszPolOpt ? NisarUpper(pszPolOpt) : std::string();
    std::string sRefDesc;
    GDALDataset* poTargetGridDS = nullptr;

    if (bIsLevel1) {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Interpolation: Level-1 product %s is not supported yet (radar coordinates).",
                 sProduct.c_str());
        GDALClose(poCoarseCubeDS);
        return nullptr;
    }
    else if (EQUAL(sProduct.c_str(), "GCOV") || EQUAL(sProduct.c_str(), "GSLC")) {
        // grids/frequency<F>/<POL>: let NisarDataset build and validate the path
        char** papszRefOO = nullptr;
        papszRefOO = CSLSetNameValue(papszRefOO, "INST", sInst.c_str());
        papszRefOO = CSLSetNameValue(papszRefOO, "FREQ", sFreq.c_str());
        if (!sPol.empty()) papszRefOO = CSLSetNameValue(papszRefOO, "POL", sPol.c_str());
        sRefDesc = CPLSPrintf("%s grids/frequency%s/%s", sProduct.c_str(), sFreq.c_str(),
                              sPol.empty() ? "<default POL>" : sPol.c_str());
        poTargetGridDS = NisarOpenInternal(sQuotedFile, papszRefOO);
        CSLDestroy(papszRefOO);
    }
    else if (EQUAL(sProduct.c_str(), "GUNW")) {
        if (sPol.empty()) sPol = "HH";
        sRefDesc = "/science/" + sInst + "/GUNW/grids/frequency" + sFreq +
                   "/unwrappedInterferogram/" + sPol + "/unwrappedPhase";
        poTargetGridDS = NisarOpenInternal(sQuotedFile + ":" + sRefDesc, nullptr);
    }
    else {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Interpolation: unsupported product type '%s'. Cannot determine reference grid.",
                 sProduct.c_str());
        GDALClose(poCoarseCubeDS);
        return nullptr;
    }

    if (!poTargetGridDS) {
        CPLError(CE_Failure, CPLE_OpenFailed,
                 "Interpolation: failed to open reference grid (%s).", sRefDesc.c_str());
        GDALClose(poCoarseCubeDS);
        return nullptr;
    }
    CPLDebug("NISAR_DRIVER", "Interpolation: reference grid %s (%d x %d)", sRefDesc.c_str(),
             poTargetGridDS->GetRasterXSize(), poTargetGridDS->GetRasterYSize());

    // Instantiate our custom dataset and copy spatial info
    NisarInterpolatedDataset* poDS = new NisarInterpolatedDataset();

    poDS->nRasterXSize = poTargetGridDS->GetRasterXSize();
    poDS->nRasterYSize = poTargetGridDS->GetRasterYSize();
    GDALGetGeoTransform(poTargetGridDS, poDS->m_adfTargetGeoTransform);

    const OGRSpatialReference* poTargetSRS = poTargetGridDS->GetSpatialRef();
    if (poTargetSRS) {
        poDS->m_oSRS = *poTargetSRS;
    }

    poDS->SetDescription(poOpenInfo->pszFilename);
    poDS->SetMetadataItem("NISAR_PRODUCT_TYPE", sProduct.c_str());
    poDS->SetMetadataItem("NISAR_CUBE_PATH", sCubePath.c_str());
    const char* pszRefPath = poTargetGridDS->GetMetadataItem("HDF5_PATH");
    poDS->SetMetadataItem("NISAR_REFERENCE_GRID",
                          pszRefPath ? pszRefPath : poTargetGridDS->GetDescription());
    if (pszQuantity) poDS->SetMetadataItem("NISAR_QUANTITY", pszQuantity);
    for (const char* pszKey : { "description", "units", "long_name" }) {
        const char* pszVal = poCoarseCubeDS->GetMetadataItem(pszKey);
        if (pszVal) poDS->SetMetadataItem(pszKey, pszVal);
    }
    GDALClose(poTargetGridDS); // Done with the target grid reference

    // Open and Warp the DEM
    CPLDebug("NISAR_DRIVER", "Interpolation: Warping DEM from %s", pszDemFile);
    poDS->m_poRawDEM = (GDALDataset*)GDALOpenEx(pszDemFile, GDAL_OF_RASTER | GDAL_OF_READONLY, nullptr, nullptr, nullptr);
    if (!poDS->m_poRawDEM) {
        CPLError(CE_Failure, CPLE_OpenFailed, "Failed to open DEM file: %s", pszDemFile);
        delete poDS; GDALClose(poCoarseCubeDS); return nullptr;
    }

    char* pszTargetWKT = nullptr;
    poDS->m_oSRS.exportToWkt(&pszTargetWKT);

    //  Fetch the user's requested resampling method (Defaulting to CUBICSPLINE)
    const char* pszResampling = CSLFetchNameValueDef(poOpenInfo->papszOpenOptions, "DEM_RESAMPLING", "CUBICSPLINE");
    GDALResampleAlg eResampleAlg = GRA_CubicSpline; // Default

    if (EQUAL(pszResampling, "NEAREST")) {
        eResampleAlg = GRA_NearestNeighbour;
    } else if (EQUAL(pszResampling, "BILINEAR")) {
        eResampleAlg = GRA_Bilinear;
    } else if (EQUAL(pszResampling, "CUBIC")) {
        eResampleAlg = GRA_Cubic;
    } else if (EQUAL(pszResampling, "CUBICSPLINE")) {
        eResampleAlg = GRA_CubicSpline;
    } else {
        CPLDebug("NISAR_DRIVER", "Unrecognized DEM_RESAMPLING '%s', falling back to CUBICSPLINE.", pszResampling);
    }

    // Lazily warped DEM on the target grid: blocks are resampled on demand, so
    // memory stays bounded even for GSLC-sized grids. Uncovered pixels read as 0.
    GDALWarpOptions* psWarpOptions = GDALCreateWarpOptions();
    psWarpOptions->hSrcDS = poDS->m_poRawDEM;
    psWarpOptions->eResampleAlg = eResampleAlg;
    GDALWarpInitDefaultBandMapping(psWarpOptions, 1);
    psWarpOptions->papszWarpOptions =
        CSLSetNameValue(psWarpOptions->papszWarpOptions, "INIT_DEST", "0");

    char** papszTO = CSLSetNameValue(nullptr, "DST_SRS", pszTargetWKT);
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
    psWarpOptions->pTransformerArg =
        GDALCreateGenImgProjTransformer2(poDS->m_poRawDEM, nullptr, papszTO);
    CSLDestroy(papszTO);
    if (psWarpOptions->pTransformerArg == nullptr) {
        CPLError(CE_Failure, CPLE_AppDefined, "Failed to create DEM -> target grid transformer.");
        GDALDestroyWarpOptions(psWarpOptions);
        CPLFree(pszTargetWKT); delete poDS; GDALClose(poCoarseCubeDS); return nullptr;
    }
    GDALSetGenImgProjTransformerDstGeoTransform(psWarpOptions->pTransformerArg,
                                                poDS->m_adfTargetGeoTransform);

    // The warped VRT takes ownership of the transformer and a reference to the raw DEM.
    poDS->m_poAlignedDEM = GDALDataset::FromHandle(GDALCreateWarpedVRT(
        poDS->m_poRawDEM, poDS->nRasterXSize, poDS->nRasterYSize,
        poDS->m_adfTargetGeoTransform, psWarpOptions));
    if (poDS->m_poAlignedDEM == nullptr) {
        CPLError(CE_Failure, CPLE_AppDefined, "Failed to create warped DEM on the target grid.");
        GDALDestroyTransformer(psWarpOptions->pTransformerArg);
        GDALDestroyWarpOptions(psWarpOptions);
        CPLFree(pszTargetWKT); delete poDS; GDALClose(poCoarseCubeDS); return nullptr;
    }
    GDALDestroyWarpOptions(psWarpOptions);
    poDS->m_poAlignedDEM->SetProjection(pszTargetWKT);
    CPLFree(pszTargetWKT);

    // Load Coarse Cube Data into RAM
    poDS->m_nCubeXSize = poCoarseCubeDS->GetRasterXSize();
    poDS->m_nCubeYSize = poCoarseCubeDS->GetRasterYSize();
    poDS->m_nCubeZSize = poCoarseCubeDS->GetRasterCount();
    GDALGetGeoTransform(poCoarseCubeDS, poDS->m_adfCubeGeoTransform);

    // Calculate the inverse transform ONCE during setup
    if (!GDALInvGeoTransform(poDS->m_adfCubeGeoTransform, poDS->m_adfCubeInvGeoTransform)) 
    {
        CPLError(CE_Failure, CPLE_AppDefined, 
                 "Cannot invert coarse cube GeoTransform during dataset initialization.");
        delete poDS; 
        GDALClose(poCoarseCubeDS);
        return nullptr; 
    }

    // Read Z-axis (heightAboveEllipsoid) from metadata
    poDS->m_zVect.reserve(poDS->m_nCubeZSize);
    for (int i = 1; i <= poDS->m_nCubeZSize; ++i) {
        GDALRasterBand* pCBand = poCoarseCubeDS->GetRasterBand(i);
        const char* pszHeight = pCBand->GetMetadataItem("Z_VALUE");
        poDS->m_zVect.push_back(pszHeight ? CPLAtof(pszHeight) : static_cast<double>(i));
    }

    // Read the entire 3D cube into a flat 1D std::vector
    int nTotalCubePixels = poDS->m_nCubeXSize * poDS->m_nCubeYSize * poDS->m_nCubeZSize;
    poDS->m_cubeData.resize(nTotalCubePixels);

    for (int z = 0; z < poDS->m_nCubeZSize; ++z) {
        GDALRasterBand* pCBand = poCoarseCubeDS->GetRasterBand(z + 1);
        float* pDst = poDS->m_cubeData.data() + (z * poDS->m_nCubeXSize * poDS->m_nCubeYSize);
        if (pCBand->RasterIO(GF_Read, 0, 0, poDS->m_nCubeXSize, poDS->m_nCubeYSize,
                             pDst, poDS->m_nCubeXSize, poDS->m_nCubeYSize,
                             GDT_Float32, 0, 0, nullptr) != CE_None)
        {
            CPLError(CE_Failure, CPLE_AppDefined, "Failed to read coarse cube data via RasterIO.");
            delete poDS;
            GDALClose(poCoarseCubeDS);
            return nullptr;
        }
    }

    GDALClose(poCoarseCubeDS); // Done with the coarse cube
    CPLDebug("NISAR_DRIVER", "Interpolation: Successfully loaded 3D cube into RAM.");

    // Create custom raster band
    poDS->SetBand(1, new NisarInterpolatedRasterBand(poDS, 1));

    return poDS;
}
