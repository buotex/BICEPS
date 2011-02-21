//
// RawFileValues.cpp
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2005 Louis Warschaw Prostate Cancer Center
//   Cedars Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//

#define RAWFILE_SOURCE

#include "RawFileValues.h"

using namespace XRawfile;
using namespace std;


namespace pwiz {
namespace raw {
namespace RawFileValues {


ValueDescriptor<ValueID_Long> ValueData<ValueID_Long>::descriptors_[] =
{
    {VersionNumber, &IXRawfile::GetVersionNumber, "VersionNumber"},
    {IsError, &IXRawfile::IsError, "IsError"},
    {IsNewFile, &IXRawfile::IsNewFile, "IsNewFile"},
    {ErrorCode, &IXRawfile::GetErrorCode, "ErrorCode"},
    {SeqRowNumber, &IXRawfile::GetSeqRowNumber, "SeqRowNumber"},
    {SeqRowSampleType, &IXRawfile::GetSeqRowSampleType, "SeqRowSampleType"},
    {InAcquisition, &IXRawfile::InAcquisition, "InAcquisition"},
    {NumberOfControllers, &IXRawfile::GetNumberOfControllers, "NumberOfControllers"},
    {NumSpectra, &IXRawfile::GetNumSpectra, "NumSpectra"},
    {NumStatusLog, &IXRawfile::GetNumStatusLog, "NumStatusLog"},
    {NumErrorLog, &IXRawfile::GetNumErrorLog, "NumErrorLog"},
    {NumTuneData, &IXRawfile::GetNumTuneData, "NumTuneData"},
    {NumTrailerExtra, &IXRawfile::GetNumTrailerExtra, "NumTrailerExtra"},
    {MaxIntensity, &IXRawfile::GetMaxIntensity, "MaxIntensity"},
    {FirstSpectrumNumber, &IXRawfile::GetFirstSpectrumNumber, "FirstSpectrumNumber"},
    {LastSpectrumNumber, &IXRawfile::GetLastSpectrumNumber, "LastSpectrumNumber"},
    {InstrumentID, &IXRawfile::GetInstrumentID, "InstrumentID"},
    {InletID, &IXRawfile::GetInletID, "InletID"},
    {ErrorFlag, &IXRawfile::GetErrorFlag, "ErrorFlag"},
    {VialNumber, &IXRawfile::GetVialNumber, "VialNumber"},
    {NumInstMethods, &IXRawfile::GetNumInstMethods, "NumInstMethods"},
    {InstNumChannelLabels, &IXRawfile::GetInstNumChannelLabels, "InstNumChannelLabels"},
    {IsThereMSData, &IXRawfile::IsThereMSData, "IsThereMSData"},
    {HasExpMethod, &IXRawfile::HasExpMethod, "HasExpMethod"},
    {FilterMassPrecision, &IXRawfile::GetFilterMassPrecision, "FilterMassPrecision"},
    {ValueID_Long_Count, 0, 0}
};


ValueDescriptor<ValueID_Double> ValueData<ValueID_Double>::descriptors_[] =
{
    {SeqRowInjectionVolume, &IXRawfile::GetSeqRowInjectionVolume, "SeqRowInjectionVolume"},
    {SeqRowSampleWeight, &IXRawfile::GetSeqRowSampleWeight, "SeqRowSampleWeight"},
    {SeqRowSampleVolume, &IXRawfile::GetSeqRowSampleVolume, "SeqRowSampleVolume"},
    {SeqRowISTDAmount, &IXRawfile::GetSeqRowISTDAmount, "SeqRowISTDAmount"},
    {SeqRowDilutionFactor, &IXRawfile::GetSeqRowDilutionFactor, "SeqRowDilutionFactor"},
    {MassResolution, &IXRawfile::GetMassResolution, "MassResolution"},
    {ExpectedRunTime, &IXRawfile::GetExpectedRunTime, "ExpectedRunTime"},
    {LowMass, &IXRawfile::GetLowMass, "LowMass"},
    {HighMass, &IXRawfile::GetHighMass, "HighMass"},
    {StartTime, &IXRawfile::GetStartTime, "StartTime"},
    {EndTime, &IXRawfile::GetEndTime, "EndTime"},
    {MaxIntegratedIntensity, &IXRawfile::GetMaxIntegratedIntensity, "MaxIntegratedIntensity"},
    {SampleVolume, &IXRawfile::GetSampleVolume, "SampleVolume"},
    {SampleWeight, &IXRawfile::GetSampleWeight, "SampleWeight"},
    {InjectionVolume, &IXRawfile::GetInjectionVolume, "InjectionVolume"},
    {ValueID_Double_Count, 0, 0}
};


ValueDescriptor<ValueID_String> ValueData<ValueID_String>::descriptors_[] =
{
    {FileName, &IXRawfile::GetFileName, "FileName"},
    {CreatorID, &IXRawfile::GetCreatorID, "CreatorID"},
    {ErrorMessage, &IXRawfile::GetErrorMessage, "ErrorMessage"},
    {WarningMessage, &IXRawfile::GetWarningMessage, "WarningMessage"},
    {SeqRowDataPath, &IXRawfile::GetSeqRowDataPath, "SeqRowDataPath"},
    {SeqRowRawFileName, &IXRawfile::GetSeqRowRawFileName, "SeqRowRawFileName"},
    {SeqRowSampleName, &IXRawfile::GetSeqRowSampleName, "SeqRowSampleName"},
    {SeqRowSampleID, &IXRawfile::GetSeqRowSampleID, "SeqRowSampleID"},
    {SeqRowComment, &IXRawfile::GetSeqRowComment, "SeqRowComment"},
    {SeqRowLevelName, &IXRawfile::GetSeqRowLevelName, "SeqRowLevelName"},
    {SeqRowInstrumentMethod, &IXRawfile::GetSeqRowInstrumentMethod, "SeqRowInstrumentMethod"},
    {SeqRowProcessingMethod, &IXRawfile::GetSeqRowProcessingMethod, "SeqRowProcessingMethod"},
    {SeqRowCalibrationFile, &IXRawfile::GetSeqRowCalibrationFile, "SeqRowCalibrationFile"},
    {SeqRowVial, &IXRawfile::GetSeqRowVial, "SeqRowVial"},
    {Flags, &IXRawfile::GetFlags, "Flags"},
    {AcquisitionFileName, &IXRawfile::GetAcquisitionFileName, "AcquisitionFileName"},
    {InstrumentDescription, &IXRawfile::GetInstrumentDescription, "InstrumentDescription"},
    {AcquisitionDate, &IXRawfile::GetAcquisitionDate, "AcquisitionDate"},
    {Operator, &IXRawfile::GetOperator, "Operator"},
    {Comment1, &IXRawfile::GetComment1, "Comment1"},
    {Comment2, &IXRawfile::GetComment2, "Comment2"},
    {SampleAmountUnits, &IXRawfile::GetSampleAmountUnits, "SampleAmountUnits"},
    {InjectionAmountUnits, &IXRawfile::GetInjectionAmountUnits, "InjectionAmountUnits"},
    {SampleVolumeUnits, &IXRawfile::GetSampleVolumeUnits, "SampleVolumeUnits"},
    {InstName, &IXRawfile::GetInstName, "InstName"},
    {InstModel, &IXRawfile::GetInstModel, "InstModel"},
    {InstSerialNumber, &IXRawfile::GetInstSerialNumber, "InstSerialNumber"},
    {InstSoftwareVersion, &IXRawfile::GetInstSoftwareVersion, "InstSoftwareVersion"},
    {InstHardwareVersion, &IXRawfile::GetInstHardwareVersion, "InstHardwareVersion"},
    {InstFlags, &IXRawfile::GetInstFlags, "InstFlags"},
    {ValueID_String_Count, 0, 0}
};


ValueData<ValueID_Long>::map_type ValueData<ValueID_Long>::descriptorMap_;
ValueData<ValueID_Double>::map_type ValueData<ValueID_Double>::descriptorMap_;
ValueData<ValueID_String>::map_type ValueData<ValueID_String>::descriptorMap_;


template<typename id_type>
void initializeMap()
{
    ValueDescriptor<id_type>* vd = ValueData<id_type>::descriptors_;
    for (;!(vd->function == 0 && vd->name == 0); vd++)
        ValueData<id_type>::descriptorMap_[vd->id] = vd;
}


void initializeMaps()
{
    initializeMap<ValueID_Long>();
    initializeMap<ValueID_Double>();
    initializeMap<ValueID_String>();
}


} // namespace RawFileValues
} // namespace raw
} // namespace pwiz
