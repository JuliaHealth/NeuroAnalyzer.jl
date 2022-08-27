#############################################################################
#
# Copyright (c) 2020 Teunis van Beelen
# All rights reserved.
#
# Email: teuniz@protonmail.com
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the copyright holder nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES
# LOSS OF USE, DATA, OR PROFITS OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#############################################################################

import sys
import io
import os
import string
import array
from collections import namedtuple
import numpy as np
from datetime import datetime

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
  print("Must be using Python version >= 3.5.0")
  sys.exit()

if np.__version__ < "1.17.0":
  print("Must be using NumPy version >= 1.17.0")
  sys.exit()

################################################################################
# START class EDFwriter
################################################################################

class EDFwriter:
  """A writer for EDF+ and BDF+ files.

  EDF header:

  offset (hex, dec) length
  ---------------------------------------------------------------------
  0x00      0     8 ascii : version of this data format (0)
  0x08      8    80 ascii : local patient identification
  0x58     88    80 ascii : local recording identification
  0xA8    168     8 ascii : startdate of recording (dd.mm.yy)
  0xB0    176     8 ascii : starttime of recording (hh.mm.ss)
  0xB8    184     8 ascii : number of bytes in header record
  0xC0    192    44 ascii : reserved
  0xEC    236     8 ascii : number of data records (-1 if unknown)
  0xF4    244     8 ascii : duration of a data record, in seconds
  0xFC    252     4 ascii : number of signals

        0x00           0     ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp)
   ns * 0x10    ns *  16     ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
   ns * 0x60    ns *  96     ns *  8 ascii : ns * physical dimension (e.g. uV or degreeC)
   ns * 0x68    ns * 104     ns *  8 ascii : ns * physical minimum (e.g. -500 or 34)
   ns * 0x70    ns * 112     ns *  8 ascii : ns * physical maximum (e.g. 500 or 40)
   ns * 0x78    ns * 120     ns *  8 ascii : ns * digital minimum (e.g. -2048)
   ns * 0x80    ns * 128     ns *  8 ascii : ns * digital maximum (e.g. 2047)
   ns * 0x88    ns * 136     ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz N:60)
   ns * 0xD8    ns * 216     ns *  8 ascii : ns * nr of samples in each data record
   ns * 0xE0    ns * 224     ns * 32 ascii : ns * reserved

  ns: number of signals

  All fields are left aligned and filled up with spaces, no NULL's.

  Only printable ASCII characters are allowed.

  Decimal separator (if any) must be a dot. No grouping characters in numbers.

  For more info about the EDF and EDF+ format, visit: https://edfplus.info/specs/

  For more info about the BDF and BDF+ format, visit: https://www.teuniz.net/edfbrowser/bdfplus%20format%20description.html

  note: In EDF, the sensitivity (e.g. uV/bit) and offset are stored using four parameters:
  digital maximum and minimum, and physical maximum and minimum.
  Here, digital means the raw data coming from a sensor or ADC. Physical means the units like uV.
  The sensitivity in units/bit is calculated as follows:

  units per bit = (physical max - physical min) / (digital max - digital min)

  The digital offset is calculated as follows:

  offset = (physical max / units per bit) - digital max

  For a better explanation about the relation between digital data and physical data,
  read the document "Coding Schemes Used with Data Converters" (PDF):

  https://www.ti.com/general/docs/lit/getliterature.tsp?baseLiteratureNumber=sbaa042

  note: An EDF file usually contains multiple so-called datarecords. One datarecord usually has a duration of one second (this is the default but it is not mandatory!).
  In that case a file with a duration of five minutes contains 300 datarecords. The duration of a datarecord can be freely choosen but, if possible, use values from
  0.1 to 1 second for easier handling. Just make sure that the total size of one datarecord, expressed in bytes, does not exceed 10MByte (15MBytes for BDF(+)).

  The RECOMMENDATION of a maximum datarecordsize of 61440 bytes in the EDF and EDF+ specification was usefull in the time people were still using DOS as their main operating system.
  Using DOS and fast (near) pointers (16-bit pointers), the maximum allocatable block of memory was 64KByte.
  This is not a concern anymore so the maximum datarecord size now is limited to 10MByte for EDF(+) and 15MByte for BDF(+). This helps to accommodate for higher samplingrates
  used by modern Analog to Digital Converters.

  EDF header character encoding: The EDF specification says that only (printable) ASCII characters are allowed.
  When writing the header info, EDFlib will assume you are using Latin1 encoding and it will automatically convert
  characters with accents, umlauts, tilde, etc. to their "normal" equivalent without the accent/umlaut/tilde/etc.
  in order to create a valid EDF file.

  The description/name of an EDF+ annotation on the other hand, is encoded in UTF-8.

  author: Teunis van Beelen
  """

  EDFLIB_TIME_DIMENSION     = 10000000
  EDFLIB_MAXSIGNALS               = 640
  EDFLIB_MAX_ANNOTATION_LEN       = 512

  EDFSEEK_SET = 0
  EDFSEEK_CUR = 1
  EDFSEEK_END = 2

  EDFLIB_FILETYPE_EDF                 = 0
  EDFLIB_FILETYPE_EDFPLUS             = 1
  EDFLIB_FILETYPE_BDF                 = 2
  EDFLIB_FILETYPE_BDFPLUS             = 3
  EDFLIB_MALLOC_ERROR                = -1
  EDFLIB_NO_SUCH_FILE_OR_DIRECTORY   = -2
  EDFLIB_FILE_CONTAINS_FORMAT_ERRORS = -3

  EDFLIB_MAXFILES_REACHED            = -4
  EDFLIB_FILE_READ_ERROR             = -5
  EDFLIB_FILE_ALREADY_OPENED         = -6
  EDFLIB_FILETYPE_ERROR              = -7
  EDFLIB_FILE_WRITE_ERROR            = -8
  EDFLIB_NUMBER_OF_SIGNALS_INVALID   = -9
  EDFLIB_FILE_IS_DISCONTINUOUS      = -10
  EDFLIB_INVALID_READ_ANNOTS_VALUE  = -11
  EDFLIB_INVALID_ARGUMENT           = -12
  EDFLIB_FILE_CLOSED                = -13

  EDFLIB_DO_NOT_READ_ANNOTATIONS = 0
  EDFLIB_READ_ANNOTATIONS        = 1
  EDFLIB_READ_ALL_ANNOTATIONS    = 2

  EDFLIB_NO_SIGNALS                 = -20
  EDFLIB_TOO_MANY_SIGNALS           = -21
  EDFLIB_NO_SAMPLES_IN_RECORD       = -22
  EDFLIB_DIGMIN_IS_DIGMAX           = -23
  EDFLIB_DIGMAX_LOWER_THAN_DIGMIN   = -24
  EDFLIB_PHYSMIN_IS_PHYSMAX         = -25
  EDFLIB_DATARECORD_SIZE_TOO_BIG    = -26

  EDFLIB_VERSION = 106

# max size of annotationtext
  __EDFLIB_WRITE_MAX_ANNOTATION_LEN = 40

# bytes in datarecord for EDF annotations, must be an integer multiple of three and two
  __EDFLIB_ANNOTATION_BYTES = 114

# for writing only
  __EDFLIB_MAX_ANNOTATION_CHANNELS = 64

  __EDFLIB_ANNOT_MEMBLOCKSZ = 1000

  __EDFAnnotationStruct = namedtuple("annotation", ["onset", "duration", "description"])

  def __init__(self, p_path: str, f_file_type: int, number_of_signals: int):
    """Creates an instance of an EDF writer.

    Path is the path to the EDF file.
    Filetype must be EDFLIB_FILETYPE_EDFPLUS or EDFLIB_FILETYPE_BDFPLUS.
    Number of signals: the number of signals that you are planning to write into the file.
    """
    self.__path = p_path
    self.__filetype = f_file_type
    self.__edfsignals = number_of_signals
    self.__status_ok = 0
    self.__edf = 0
    self.__bdf = 0
    self.__plus_patientcode = ""
    self.__plus_gender = ""
    self.__plus_birthdate = ""
    self.__plus_patient_name = ""
    self.__plus_patient_additional = ""
    self.__plus_startdate = ""
    self.__plus_admincode = ""
    self.__plus_technician = ""
    self.__plus_equipment = ""
    self.__plus_recording_additional = ""
    self.__annotationslist = []
    self.__nr_annot_chns = 1
    self.__long_data_record_duration = self.EDFLIB_TIME_DIMENSION
    self.__annotlist_sz = 0
    self.__annots_in_file = 0
    self.__plus_gender = 2
    self.__datarecords = 0
    self.__recordsize = 0
    self.__signal_write_sequence_pos = 0
    self.__total_annot_bytes = 0
    self.__startdate_year = 0
    self.__startdate_month = 0
    self.__startdate_day = 0
    self.__starttime_hour = 0
    self.__starttime_minute = 0
    self.__starttime_second = 0
    self.__starttime_offset = 0
    self.__plus_birthdate_year = 0
    self.__plus_birthdate_month = 0
    self.__plus_birthdate_day = 0

    if sys.version_info[0] != 3 or sys.version_info[1] < 5:
      raise EDFexception("Must be using Python version >= 3.5.0")

    if (self.__edfsignals < 1) or (self.__edfsignals > self.EDFLIB_MAXSIGNALS):
      raise EDFexception("Invalid number of signals.")

    if (self.__filetype != self.EDFLIB_FILETYPE_EDFPLUS) and (self.__filetype != self.EDFLIB_FILETYPE_BDFPLUS):
      raise EDFexception("Invalid filetype.")

    try:
      self.__file_out = open(self.__path, "wb")
    except OSError as e:
      raise EDFexception("Can not open file for writing: %s" %(e.strerror))

    if self.__filetype == self.EDFLIB_FILETYPE_EDFPLUS:
      self.__edf = 1
    else:
      self.__bdf = 1

    self.__param_label = [""] * self.__edfsignals
    self.__param_transducer = [""] * self.__edfsignals
    self.__param_physdimension = [""] * self.__edfsignals
    self.__param_phys_min = [0.0] * self.__edfsignals
    self.__param_phys_max = [0.0] * self.__edfsignals
    self.__param_dig_min = [0] * self.__edfsignals
    self.__param_dig_max = [0] * self.__edfsignals
    self.__param_prefilter = [""] * self.__edfsignals
    self.__param_smp_per_record = [0] * self.__edfsignals
    self.__param_offset = [0.0] * self.__edfsignals
    self.__param_buf_offset = [0] * self.__edfsignals
    self.__param_bitvalue = [0.0] * self.__edfsignals

    self.__status_ok = 1

  def close(self) -> int:
    """Finalizes and closes the file.

    This function is required after writing. Failing to do so will cause a corrupted and incomplete file.
    Returns 0 on success, otherwise -1.
    """
    if self.__status_ok:
      if self.__datarecords < 100000000:
        self.__file_out.seek(236, io.SEEK_SET)
        if self.__fprint_int_number_nonlocalized(self.__file_out, self.__datarecords, 0, 0) < 2:
          self.__file_out.write(bytes(" ", encoding="ascii"))
      self.__write_annotations()
      self.__file_out.close()
      self.__status_ok = 0
      return 0
    else:
      return -1

  def version(self) -> int:
    """If version is 1.00 then it will return 100."""
    return self.EDFLIB_VERSION

  def setSampleFrequency(self, s: int, sf: int) -> int:
    """Sets the samplefrequency of signal s.

    (In reallity, it sets the number of samples in a datarecord.)
    The samplefrequency of a signal is determined as: sf = number of samples in a datarecord / datarecord duration.
    The samplefrequency equals the number of samples in a datarecord only when the datarecord duration is set to the default of one second.
    This function is required for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    sf is the samplefrequency.
    Returns 0 on success, otherwise -1.
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0) or (sf < 1):
      return -1
    self.__param_smp_per_record[s] = sf
    return 0

  def setPhysicalMaximum(self, s: int, phys_max: float) -> int:
    """Sets the maximum physical value of signal s.

    This is the value of the input of the ADC when the output equals the value of "digital maximum".
    It is the highest value that the equipment is able to record. It does not necessarily mean the signal recorded reaches this level.
    Must be un-equal to physical minimum.
    This function is required for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    phys_max is the maximum input value.

    Returns 0 on success, otherwise -1.

    note: In EDF, the sensitivity (e.g. uV/bit) and offset are stored using four parameters:
    digital maximum and minimum, and physical maximum and minimum.
    Here, digital means the raw data coming from a sensor or ADC. Physical means the units like uV.
    Usually they are the extreme input and output values of the ADC.
    The sensitivity in units/bit is calculated as follows:

    units per bit = (physical max - physical min) / (digital max - digital min)

    The digital offset is calculated as follows:

    offset = (physical max / units per bit) - digital max
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_phys_max[s] = phys_max
    return 0

  def setPhysicalMinimum(self, s: int, phys_min: float) -> int:
    """Sets the minimum physical value of signal s.

    This is the value of the input of the ADC when the output equals the value of "digital minimum".
    It is the lowest value that the equipment is able to record. It does not necessarily mean the signal recorded reaches this level.
    Must be un-equal to physical maximum.
    This function is required for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    phys_min is the minimum input value.

    Returns 0 on success, otherwise -1.

    note: In EDF, the sensitivity (e.g. uV/bit) and offset are stored using four parameters:
    digital maximum and minimum, and physical maximum and minimum.
    Here, digital means the raw data coming from a sensor or ADC. Physical means the units like uV.
    Usually they are the extreme input and output values of the ADC.
    The sensitivity in units/bit is calculated as follows:

    units per bit = (physical max - physical min) / (digital max - digital min)

    The digital offset is calculated as follows:

    offset = (physical max / units per bit) - digital max
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_phys_min[s] = phys_min
    return 0

  def setDigitalMaximum(self, s: int, dig_max: int) -> int:
    """Sets the maximum digital value of signal s.

    Usually it's the extreme output of the ADC. The maximum value is 32767 for EDF and 8388607 for BDF.
    It is the highest value that the equipment is able to record. It does not necessarily mean the signal recorded reaches this level.
    Must be higher than digital minimum.
    This function is required for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    dig_max is the maximum output value (<= 32767 for EDF and <= 8388607 for BDF).

    Returns 0 on success, otherwise -1.

    note: In EDF, the sensitivity (e.g. uV/bit) and offset are stored using four parameters:
    digital maximum and minimum, and physical maximum and minimum.
    Here, digital means the raw data coming from a sensor or ADC. Physical means the units like uV.
    Usually they are the extreme input and output values of the ADC.
    The sensitivity in units/bit is calculated as follows:

    units per bit = (physical max - physical min) / (digital max - digital min)

    The digital offset is calculated as follows:

    offset = (physical max / units per bit) - digital max
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    if self.__edf != 0:
      if dig_max > 32767:
        return -1
    else:
      if dig_max > 8388607:
        return -1
    self.__param_dig_max[s] = dig_max
    return 0

  def setDigitalMinimum(self, s: int, dig_min: int) -> int:
    """Sets the minimum digital value of signal s.

    Usually it's the extreme output of the ADC. The minimum value is -32768 for EDF and -8388608 for BDF.
    It is the lowest value that the equipment is able to record. It does not necessarily mean the signal recorded reaches this level.
    Must be lower than digital maximum.
    This function is required for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    dig_min is the minimum output value (>= -32768 for EDF and >= -8388608 for BDF).

    Returns 0 on success, otherwise -1.

    note: In EDF, the sensitivity (e.g. uV/bit) and offset are stored using four parameters:
    digital maximum and minimum, and physical maximum and minimum.
    Here, digital means the raw data coming from a sensor or ADC. Physical means the units like uV.
    Usually they are the extreme input and output values of the ADC.
    The sensitivity in units/bit is calculated as follows:

    units per bit = (physical max - physical min) / (digital max - digital min)

    The digital offset is calculated as follows:

    offset = (physical max / units per bit) - digital max
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    if self.__edf != 0:
      if dig_min < -32768:
        return -1
    else:
      if dig_min < -8388608:
        return -1
    self.__param_dig_min[s] = dig_min
    return 0

  def setSignalLabel(self, s: int, label: str) -> int:
    """Sets the label (name) of signal s.

    (e.g. "FP1", "SaO2", etc.) String must contain printable ASCII only.
    This function is recommended for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    label is the signallabel.
    Returns 0 on success, otherwise -1.
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_label[s] = label
    return 0

  def setPreFilter(self, s: int, prefilter: str) -> int:
    """Sets the prefilter description of signal s.

    (e.g. "HP:0.05Hz", "LP:250Hz", "N:60Hz", etc.) String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.

    s is the signal number (zero-based).
    prefilter is the prefilter description.
    Returns 0 on success, otherwise -1.
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_prefilter[s] = prefilter
    return 0

  def setTransducer(self, s: int, transducer: str) -> int:
    """Sets the transducer description of signal s.

    ("AgAgCl cup electrodes", etc.) String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.

    s is the signal number (zero-based).
    transducer is the transducer description.
    Returns 0 on success, otherwise -1.
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_transducer[s] = transducer
    return 0

  def setPhysicalDimension(self, s: int, physical_dimension: str) -> int:
    """Sets the physical_dimension (unit) of signal s.

    ("uV", "BPM", "mA", "Degr.", etc.) String must contain printable ASCII only.
    This function recommended for every signal and can be called only before the first sample write action.

    s is the signal number (zero-based).
    physical_dimension is the physical dimension description.
    Returns 0 on success, otherwise -1.
    """
    if (s < 0) or (s >= self.__edfsignals) or (self.__datarecords != 0):
      return -1
    self.__param_physdimension[s] = physical_dimension
    return 0

  def setStartDateTime(self, year: int, month: int, day: int, hour: int, minute: int, second: int, subsecond: int) -> int:
    """
    Sets the startdate and starttime.

    If not called, the system date and time at runtime will be used.
    This function is optional and can be called only before the first sample write action.
    If subsecond precision is not needed or not applicable, leave it at zero.
    Note: for anonymization purposes, the consensus is to use 1985-01-01 00:00:00 for the startdate and starttime.

    year: 1985 - 2084
    month: 1 - 12
    day: 1 - 31
    hour: 0 - 23
    minute: 0 - 59
    second: 0 - 59
    subsecond: 0 - 9999  expressed in units of 100 microSeconds

    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    if (year < 1985)   or (year > 2084) or \
       (month < 1)     or (month > 12)  or \
       (day < 1)       or (day > 31)    or \
       (hour < 0)      or (hour > 23)   or \
       (minute < 0)    or (minute > 59) or \
       (second < 0)    or (second > 59) or \
       (subsecond < 0) or (subsecond > 9999):
      return -1
    self.__startdate_year = year
    self.__startdate_month = month
    self.__startdate_day = day
    self.__starttime_hour = hour
    self.__starttime_minute = minute
    self.__starttime_second = second
    self.__starttime_offset = subsecond * 1000
    return 0

  def setPatientName(self, name: str) -> int:
    """Sets the patientname.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_patient_name = name
    return 0

  def setPatientCode(self, code: str) -> int:
    """Sets the patientcode.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_patientcode = code
    return 0

  def setPatientGender(self, gender: int) -> int:
    """Sets the patient's gender.

    gender: 0 = female, 1 = male, 2 = unknown or not applicable (default)
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    if (gender < 0) or (gender > 2):
      return -1
    self.__plus_gender = gender
    return 0

  def setPatientBirthDate(self, year: int, month: int, day: int) -> int:
    """Sets the patients' birthdate.

    This function is optional and can be called only before the first sample write action.

    year: 1800 - 3000
    month: 1 - 12
    day: 1 - 31

    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    if (year < 1800)   or (year > 3000) or \
       (month < 1)     or (month > 12)  or \
       (day < 1)       or (day > 31):
      return -1
    self.__plus_birthdate_year = year
    self.__plus_birthdate_month = month
    self.__plus_birthdate_day = day
    return 0

  def setAdditionalPatientInfo(self, additional: str) -> int:
    """Sets the additional information related to the patient.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_patient_additional = additional
    return 0

  def setAdministrationCode(self, admin_code: str) -> int:
    """Sets the administration code.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_admincode = admin_code
    return 0

  def setTechnician(self, technician: str) -> int:
    """Sets the name or id of the technician who performed the recording.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_technician = technician
    return 0

  def setEquipment(self, equipment: str) -> int:
    """Sets the description of the equipment used for the recording.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_equipment = equipment
    return 0

  def setAdditionalRecordingInfo(self, additional: str) -> int:
    """Sets the additional info related to the recording.

    String must contain printable ASCII only.
    This function is optional and can be called only before the first sample write action.
    Returns 0 on success, otherwise -1.
    """
    if self.__datarecords != 0:
      return -1
    self.__plus_recording_additional = additional
    return 0

  def writeSamples(self, buf: np.array) -> int:
    """Write samples.

    Writes sf samples into the file.
    Buf must be a one-dimensional numpy array containing samples of one signal of datatype int32, float_ or float64.
    For EDF, dataype int16 can also be used.
    If buf is of type integer, the samples are written into the file without any conversion.
    If buf is of type float, the physical samples will be converted to digital samples using the
    values of physical maximum, physical minimum, digital maximum and digital minimum.
    The number of samples written is equal to the samplefrequency of the signal.
    (actually, it's the value that is set with setSampleFrequency()).
    Size of buf should be equal to or bigger than the samplefrequency.
    Call this function for every signal in the file. The order is important!
    When there are 4 signals in the file,  the order of calling this function
    must be: signal 0, signal 1, signal 2, signal 3, signal 0, signal 1, signal 2, etc.
    The end of a recording must always be at the end of a complete cycle.

    buf is a one-dimensional numpy array of datatype int32, float_ or float64. For EDF, dataype int16 can also be used.
    Returns 0 on success, otherwise -1.
    """
    if self.__status_ok == 0:
      return -1

    if buf.ndim != 1:
      return -1

    if (buf.dtype != np.int16) and (buf.dtype != np.int32) and (buf.dtype != np.float_) and (buf.dtype != np.float64):
      return -1

    if (buf.dtype == np.int16) and (self.__bdf != 0):
      return -1

    edfsignal = self.__signal_write_sequence_pos

    if self.__datarecords == 0:
      if edfsignal == 0:
        error = self.__write_edf_header()
        if error != 0:
          return error

    sf = self.__param_smp_per_record[edfsignal]

    digmax = self.__param_dig_max[edfsignal]

    digmin = self.__param_dig_min[edfsignal]

    if sf > buf.size:
      return -1

    if self.__edf != 0:
      if (buf.dtype == np.int16) or (buf.dtype == np.int32):
        for i in range(0, sf):
          if buf[i] > digmax:
            buf[i] = digmax
          if buf[i] < digmin:
            buf[i] = digmin
          self.__file_out.write(buf[i].astype("int16").tobytes(order="C"))
      else:
        for i in range(0, sf):
          value = int((buf[i] / self.__param_bitvalue[edfsignal]) - self.__param_offset[edfsignal])
          if value > digmax:
            value = digmax
          if value < digmin:
            value = digmin
          self.__file_out.write(value.to_bytes(2, byteorder="little", signed=True))
    else:
      if buf.dtype == np.int32:
        for i in range(0, sf):
          value = int(buf[i])
          if value > digmax:
            value = digmax
          if value < digmin:
            value = digmin
          self.__file_out.write(value.to_bytes(3, byteorder="little", signed=True))
      else:
        for i in range(0, sf):
          value = int((buf[i] / self.__param_bitvalue[edfsignal]) - self.__param_offset[edfsignal])
          if value > digmax:
            value = digmax
          if value < digmin:
            value = digmin
          self.__file_out.write(value.to_bytes(3, byteorder="little", signed=True))

    self.__signal_write_sequence_pos += 1

    if self.__signal_write_sequence_pos == self.__edfsignals:
      self.__signal_write_sequence_pos = 0

      if self.__write_tal(self.__file_out) != 0:
        return -1

      self.__datarecords += 1

    return 0

  def setDataRecordDuration(self, duration: int) -> int:
    """Sets the datarecord duration.

    This function is optional, normally you don't need to change the default value of one second.
    This function is NOT REQUIRED but can be called only before the first sample write action.

    This function can be used when you want to use a non-integer samplerate.
    For example, if you want to use a samplerate of 0.5 Hz, set the samplefrequency to 5 Hz and
    the datarecord duration to 10 seconds, or alternatively, set the samplefrequency to 1 Hz and
    the datarecord duration to 2 seconds.
    This function can also be used when you want to use a very high samplerate.
    For example, if you want to use a samplerate of 5 GHz,
    set the samplefrequency to 5000 Hz and the datarecord duration to 1 microSecond.
    Do not use this function if not necessary.

    duration is expressed in microSeconds, range: 1 - 60000000  (1uSec. - 60 sec.)
    Returns 0 on success, otherwise -1.
    """
    if (duration < 1) or (duration > 60000000) or (self.__datarecords != 0):
      return -1

    self.__long_data_record_duration = duration * 10

    return 0

  def setNumberOfAnnotationSignals(self, annot_signals: int) -> int:
    """Sets the number of annotation signals.

    The default value is 1.
    This function is optional and, if used, must be called before the first sample write action.
    Normally you don't need to change the default value. Only when the number of annotations
    you want to write is higher than the number of datarecords in the recording, you can use
    this function to increase the storage space for annotations.
    """
    if (annot_signals < 1) or (annot_signals >= self.__EDFLIB_MAX_ANNOTATION_CHANNELS) or (self.__datarecords != 0):
      return -1

    self.__nr_annot_chns = annot_signals

    return 0

  def writeAnnotation(self, onset: int, duration: int, description: str) -> int:
    """Writes an annotation/event to the file.

    onset is relative to the starttime of the recording and must be >= 0.
    onset and duration are in units of 100 microSeconds. Resolution is 0.0001 second.
    E.g. 34.071 seconds must be written as 340710.
    If duration is unknown or not applicable: set a negative number (-1).
    Description is a string containing the text that describes the event.
    This function is optional.
    """
    if (self.__status_ok == 0) or (onset < 0):
      return -1

    self.__annotationslist.append(self.__EDFAnnotationStruct(onset = onset, duration = duration, description = description))

    self.__annots_in_file += 1

    return 0

################################################################################
# from here only internal utils
################################################################################

# writes the EDF header
  def __write_edf_header(self):
    if self.__status_ok == 0:
      return -1

    eq_sf = 1

    self.__recordsize = 0

    str_ = bytearray(256)

    self.__total_annot_bytes = self.__EDFLIB_ANNOTATION_BYTES * self.__nr_annot_chns

    for i in range(0, self.__edfsignals):
      if self.__param_smp_per_record[i] < 1:
        return self.EDFLIB_NO_SAMPLES_IN_RECORD

      if self.__param_dig_max[i] == self.__param_dig_min[i]:
        return self.EDFLIB_DIGMIN_IS_DIGMAX

      if self.__param_dig_max[i] < self.__param_dig_min[i]:
        return self.EDFLIB_DIGMAX_LOWER_THAN_DIGMIN

      if self.__param_phys_max[i] == self.__param_phys_min[i]:
        return self.EDFLIB_PHYSMIN_IS_PHYSMAX

      self.__recordsize += self.__param_smp_per_record[i]

      if i > 0:
        if self.__param_smp_per_record[i] != self.__param_smp_per_record[i-1]:
          eq_sf = 0

    if self.__edf != 0:
      self.__recordsize *= 2
      self.__recordsize += self.__total_annot_bytes
      if self.__recordsize > (10 * 1024 * 1024):  # datarecord size should not exceed 10MB for EDF
        return self.EDFLIB_DATARECORD_SIZE_TOO_BIG
          # if your application gets hit by this limitation, lower the value for the datarecord duration
          # using the function edf_set_datarecord_duration()
    else:
      self.__recordsize *= 3
      self.__recordsize += self.__total_annot_bytes
      if self.__recordsize > (15 * 1024 * 1024):  #datarecord size should not exceed 15MB for BDF
        return self.EDFLIB_DATARECORD_SIZE_TOO_BIG
          # if your application gets hit by this limitation, lower the value for the datarecord duration
          # using the function edf_set_datarecord_duration()

    for i in range(0, self.__edfsignals):
      self.__param_bitvalue[i] = (self.__param_phys_max[i] - self.__param_phys_min[i]) / (self.__param_dig_max[i] - self.__param_dig_min[i])
      self.__param_offset[i] = self.__param_phys_max[i] / self.__param_bitvalue[i] - self.__param_dig_max[i]

    total_signals = self.__edfsignals + self.__nr_annot_chns

    hdr_sz = (total_signals + 1) * 256
    header = bytearray(hdr_sz)
    for i in range(0, hdr_sz):
      header[i] = 32

    self.__file_out.seek(0, io.SEEK_SET)

    if self.__edf != 0:
      header[0 : 8] = bytes("0       ", encoding="ascii")
    else:
      header[0 : 8] = bytes("\xffBIOSEMI", encoding="latin_1")

    p = 0
    if self.__plus_birthdate_year == 0:
      rest = 73
    else:
      rest = 63

    self.__plus_patientcode.lstrip(" ")
    self.__plus_patientcode.rstrip(" ")
    l = len(self.__plus_patientcode)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_patientcode, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
        rest = 0
      else:
        rest -= l
      self.__latin1_to_ascii(str_, l)
      for i in range(0 , l):
        if str_[i] == 32:
          str_[i] = 95
      header[8 + p : 8 + p + l] = str_[0 : l]
      p += l
      header[8 + p] = 32
      p += 1
    else:
      header[8 + p : 8 + p + 2] = bytes("X ", encoding="ascii")
      p += 2
      rest -= 1

    if self.__plus_gender == 1:
      header[8 + p : 8 + p + 2] = bytes("M ", encoding="ascii")
    else:
      if self.__plus_gender == 0:
        header[8 + p : 8 + p + 2] = bytes("F ", encoding="ascii")
      else:
        header[8 + p : 8 + p + 2] = bytes("X ", encoding="ascii")
    p += 2

    if self.__plus_birthdate_year == 0:
      header[8 + p : 8 + p + 2] = bytes("X ", encoding="ascii")
      p += 2
    else:
      month_str = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
      header[8 + p : 8 + p + 12] = bytes("%02d-%s-%04d " \
      %(self.__plus_birthdate_day, month_str[self.__plus_birthdate_month - 1], self.__plus_birthdate_year), encoding="ascii")
      p += 12

    self.__plus_patient_name.lstrip(" ")
    self.__plus_patient_name.rstrip(" ")
    l = len(self.__plus_patient_name)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_patient_name, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
        rest = 0
      else:
        rest -= l
      self.__latin1_to_ascii(str_, l)
      for i in range(0 , l):
        if str_[i] == 32:
          str_[i] = 95
      header[8 + p : 8 + p + l] = str_[0 : l]
      p += l
      header[8 + p] = 32
      p += 1
    else:
      header[8 + p : 8 + p + 2] = bytes("X ", encoding="ascii")
      p += 2
      rest -= 1

    self.__plus_patient_additional.lstrip(" ")
    self.__plus_patient_additional.rstrip(" ")
    l = len(self.__plus_patient_additional)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_patient_additional, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
      self.__latin1_to_ascii(str_, l)
      header[8 + p : 8 + p + l] = str_[0 : l]
      p += l

    for j in range(p, 80):
      header[8 + j] = 32

    if self.__startdate_year == 0:
      dt = datetime.now()
      self.__startdate_year = dt.year
      self.__startdate_month = dt.month
      self.__startdate_day = dt.day
      self.__starttime_hour = dt.hour
      self.__starttime_minute = dt.minute
      self.__starttime_second = dt.second

    month_str = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    header[88 : 110] = bytes("Startdate %02d-%s-%04d " \
    %(self.__startdate_day, month_str[self.__startdate_month - 1], self.__startdate_year), encoding="ascii")

    p = 22
    rest = 50

    self.__plus_admincode.lstrip(" ")
    self.__plus_admincode.rstrip(" ")
    l = len(self.__plus_admincode)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_admincode, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
        rest = 0
      else:
        rest -= l
      self.__latin1_to_ascii(str_, l)
      for i in range(0 , l):
        if str_[i] == 32:
          str_[i] = 95
      header[88 + p : 88 + p + l] = str_[0 : l]
      p += l
      header[88 + p] = 32
      p += 1
    else:
      header[88 + p : 88 + p + 2] =  bytes("X ", encoding="ascii")
      p += 2
      rest -= 1

    self.__plus_technician.lstrip(" ")
    self.__plus_technician.rstrip(" ")
    l = len(self.__plus_technician)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_technician, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
        rest = 0
      else:
        rest -= l
      self.__latin1_to_ascii(str_, l)
      for i in range(0 , l):
        if str_[i] == 32:
          str_[i] = 95
      header[88 + p : 88 + p + l] = str_[0 : l]
      p += l
      header[88 + p] = 32
      p += 1
    else:
      header[88 + p : 88 + p + 2] =  bytes("X ", encoding="ascii")
      p += 2
      rest -= 1

    self.__plus_equipment.lstrip(" ")
    self.__plus_equipment.rstrip(" ")
    l = len(self.__plus_equipment)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_equipment, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
        rest = 0
      else:
        rest -= l
      self.__latin1_to_ascii(str_, l)
      for i in range(0 , l):
        if str_[i] == 32:
          str_[i] = 95
      header[88 + p : 88 + p + l] = str_[0 : l]
      p += l
      header[88 + p] = 32
      p += 1
    else:
      header[88 + p : 88 + p + 2] =  bytes("X ", encoding="ascii")
      p += 2
      rest -= 1

    self.__plus_recording_additional.lstrip(" ")
    self.__plus_recording_additional.rstrip(" ")
    l = len(self.__plus_recording_additional)
    if (l != 0) and (rest > 0):
      str_ = bytearray(self.__plus_recording_additional, encoding="latin_1")
      l = len(str_)
      if l > rest:
        l = rest
      self.__latin1_to_ascii(str_, l)
      header[88 + p : 88 + p + l] = str_[0 : l]
      p += l

    header[168 : 168 + 16] = bytes("%02d.%02d.%02d%02d.%02d.%02d" \
    %(self.__startdate_day, self.__startdate_month, self.__startdate_year % 100, self.__starttime_hour, self.__starttime_minute, self.__starttime_second), encoding="ascii")

    str_ = bytearray(256)
    l = self.__sprint_int_number_nonlocalized(str_, (total_signals + 1) * 256, 0, 0)
    if l > 8:
      l = 8
    header[184 : 184 + l] = str_[0 : l]
    if self.__edf != 0:
      header[192 : 192 + 5] = bytes("EDF+C", encoding="ascii")
    else:
      header[192 : 192 + 5] = bytes("BDF+C", encoding="ascii")
    header[236 : 236 + 8] = bytes("-1      ", encoding="ascii")
    if self.__long_data_record_duration == self.EDFLIB_TIME_DIMENSION:
      header[244 : 244 + 8] = bytes("1       ", encoding="ascii")
    else:
      l = self.__sprint_number_nonlocalized(str_, self.__long_data_record_duration / self.EDFLIB_TIME_DIMENSION)
      if l > 8:
        l = 8
      header[244 : 244 + l] = str_[0 : l]
    l = self.__sprint_int_number_nonlocalized(str_, total_signals, 0, 0)
    if l > 4:
      l = 4
    header[252 : 252 + l] = str_[0 : l]

    for i in range(0, self.__edfsignals):
      l = len(self.__param_label[i])
      if l != 0:
        str_ = bytearray(self.__param_label[i], encoding="latin_1")[0 : 16]
        l = len(str_)
        if l > 16:
          l = 16
        self.__latin1_to_ascii(str_, l)
        header[256 + (i * 16) : 256 + (i * 16) + l] = str_[0 : l]
    for i in range(self.__edfsignals, total_signals):
      if(self.__edf != 0):
        header[256 + (i * 16) : 256 + (i * 16) + 16] = bytes("EDF Annotations ", encoding="ascii")
      else:
        header[256 + (i * 16) : 256 + (i * 16) + 16] = bytes("BDF Annotations ", encoding="ascii")

    for i in range(0, self.__edfsignals):
      l = len(self.__param_transducer[i])
      if l != 0:
        str_ = bytearray(self.__param_transducer[i], encoding="latin_1")[0 : 80]
        l = len(str_)
        if l > 80:
          l = 80
        self.__latin1_to_ascii(str_, l)
        header[256 + (total_signals * 16) + (i * 80) : 256 + (total_signals * 16) + (i * 80) + l] = str_[0 : l]

    for i in range(0, self.__edfsignals):
      l = len(self.__param_physdimension[i])
      if l != 0:
        str_ = bytearray(self.__param_physdimension[i], encoding="latin_1")[0 : 8]
        l = len(str_)
        if l > 8:
          l = 8
        self.__latin1_to_ascii(str_, l)
        header[256 + (total_signals * 96) + (i * 8) : 256 + (total_signals * 96) + (i * 8) + l] = str_[0 : l]

    str_ = bytearray(256)
    for i in range(0, self.__edfsignals):
      l = self.__sprint_number_nonlocalized(str_, self.__param_phys_min[i])
      if l > 8:
        l = 8
      header[256 + (total_signals * 104) + (i * 8) : 256 + (total_signals * 104) + (i * 8) + l] = str_[0 : l]

    for i in range(self.__edfsignals, total_signals):
        header[256 + (total_signals * 104) + (i * 8) : 256 + (total_signals * 104) + (i * 8) + 2] = bytes("-1", encoding="ascii")

    for i in range(0, self.__edfsignals):
      l = self.__sprint_number_nonlocalized(str_, self.__param_phys_max[i])
      if l > 8:
        l = 8
      header[256 + (total_signals * 112) + (i * 8) : 256 + (total_signals * 112) + (i * 8) + l] = str_[0 : l]

    for i in range(self.__edfsignals, total_signals):
        header[256 + (total_signals * 112) + (i * 8) : 256 + (total_signals * 112) + (i * 8) + 1] = bytes("1", encoding="ascii")

    for i in range(0, self.__edfsignals):
      l = self.__sprint_int_number_nonlocalized(str_, self.__param_dig_min[i], 0, 0)
      if l > 8:
        l = 8
      header[256 + (total_signals * 120) + (i * 8) : 256 + (total_signals * 120) + (i * 8) + l] = str_[0 : l]

    for i in range(self.__edfsignals, total_signals):
      if self.__edf != 0:
        header[256 + (total_signals * 120) + (i * 8) : 256 + (total_signals * 120) + (i * 8) + 6] = bytes("-32768", encoding="ascii")
      else:
        header[256 + (total_signals * 120) + (i * 8) : 256 + (total_signals * 120) + (i * 8) + 8] = bytes("-8388608", encoding="ascii")

    for i in range(0, self.__edfsignals):
      l = self.__sprint_int_number_nonlocalized(str_, self.__param_dig_max[i], 0, 0)
      if l > 8:
        l = 8
      header[256 + (total_signals * 128) + (i * 8) : 256 + (total_signals * 128) + (i * 8) + l] = str_[0 : l]

    for i in range(self.__edfsignals, total_signals):
      if self.__edf != 0:
        header[256 + (total_signals * 128) + (i * 8) : 256 + (total_signals * 128) + (i * 8) + 5] = bytes("32767", encoding="ascii")
      else:
        header[256 + (total_signals * 128) + (i * 8) : 256 + (total_signals * 128) + (i * 8) + 7] = bytes("8388607", encoding="ascii")

    for i in range(0, self.__edfsignals):
      l = len(self.__param_prefilter[i])
      if l != 0:
        str_ = bytearray(self.__param_prefilter[i], encoding="latin_1")[0 : 80]
        l = len(str_)
        if l > 80:
          l = 80
        self.__latin1_to_ascii(str_, l)
        header[256 + (total_signals * 136) + (i * 80) : 256 + (total_signals * 136) + (i * 80) + l] = str_[0 : l]

    str_ = bytearray(256)
    for i in range(0, self.__edfsignals):
      l = self.__sprint_int_number_nonlocalized(str_, self.__param_smp_per_record[i], 0, 0)
      if l > 8:
        l = 8
      header[256 + (total_signals * 216) + (i * 8) : 256 + (total_signals * 216) + (i * 8) + l] = str_[0 : l]

    for i in range(self.__edfsignals, total_signals):
      if self.__edf != 0:
        l = self.__sprint_int_number_nonlocalized(str_, self.__EDFLIB_ANNOTATION_BYTES // 2, 0, 0)
      else:
        l = self.__sprint_int_number_nonlocalized(str_, self.__EDFLIB_ANNOTATION_BYTES // 3, 0, 0)
      if l > 8:
        l = 8
      header[256 + (total_signals * 216) + (i * 8) : 256 + (total_signals * 216) + (i * 8) + l] = str_[0 : l]

    self.__file_out.write(header)

    return 0

# writes a TAL
  def __write_tal(self, f):
    scratchpad = bytearray(self.__total_annot_bytes)

    p = self.__snprint_ll_number_nonlocalized(scratchpad, 0, (self.__datarecords * self.__long_data_record_duration + self.__starttime_offset) / self.EDFLIB_TIME_DIMENSION, 0, 1)
    if ((self.__long_data_record_duration % self.EDFLIB_TIME_DIMENSION) != 0) or (self.__starttime_offset != 0):
      scratchpad[p] = 46
      p += 1
      p += self.__snprint_ll_number_nonlocalized(scratchpad, p, (self.__datarecords * self.__long_data_record_duration + self.__starttime_offset) % self.EDFLIB_TIME_DIMENSION, 7, 0)
    scratchpad[p] = 20
    p += 1
    scratchpad[p] = 20
    p += 1
    for i in range(p, self.__total_annot_bytes):
      scratchpad[i] = 0
    f.write(scratchpad)

    return 0

# writes the annotations to the file
  def __write_annotations(self):
    err = 0
    datrecs = 0

    str_ = bytearray(self.__EDFLIB_ANNOTATION_BYTES)

    offset = (self.__edfsignals + self.__nr_annot_chns + 1) * 256

    file_sz = offset + (self.__datarecords * self.__recordsize)

    datrecsize = self.__total_annot_bytes

    for i in range(0, self.__edfsignals):
      if self.__edf != 0:
        offset += self.__param_smp_per_record[i] * 2

        datrecsize += self.__param_smp_per_record[i] * 2
      else:
        offset += self.__param_smp_per_record[i] * 3

        datrecsize += self.__param_smp_per_record[i] * 3

    j = 0

    for k in range(0, self.__annots_in_file):
      annot2 = self.__annotationslist[k]

      onset = annot2.onset + (self.__starttime_offset // 1000)

      p = 0

      if j == 0:  # first annotation signal
        if (offset + self.__total_annot_bytes) > file_sz:
          break

        self.__file_out.seek(offset, io.SEEK_SET)

        p += self.__snprint_ll_number_nonlocalized(str_, 0, (datrecs * self.__long_data_record_duration + self.__starttime_offset) // self.EDFLIB_TIME_DIMENSION, 0, 1)

        if ((self.__long_data_record_duration % self.EDFLIB_TIME_DIMENSION) != 0) or (self.__starttime_offset != 0):
          str_[p] = 46
          p += 1
          n = self.__snprint_ll_number_nonlocalized(str_, p, (datrecs * self.__long_data_record_duration + self.__starttime_offset) % self.EDFLIB_TIME_DIMENSION, 7, 0)
          p += n
        str_[p] = 20
        p += 1
        str_[p] = 20
        p += 1
        str_[p] =  0
        p += 1

      n = self.__snprint_ll_number_nonlocalized(str_, p, onset // 10000, 0, 1)
      p += n
      if (onset % 10000) != 0:
        str_[p] = 46
        p += 1
        n = self.__snprint_ll_number_nonlocalized(str_, p, onset % 10000, 4, 0)
        p += n
      if annot2.duration >= 0:
        str_[p] = 21
        p += 1
        n = self.__snprint_ll_number_nonlocalized(str_, p, annot2.duration // 10000, 0, 0)
        p += n
        if (annot2.duration % 10000) != 0:
          str_[p] = 46
          p += 1
          n = self.__snprint_ll_number_nonlocalized(str_, p, annot2.duration % 10000, 4, 0)
          p += n
      str_[p] = 20
      p += 1
      annot2.description.lstrip(" ")
      annot2.description.rstrip(" ")
      ba_tmp = bytearray(annot2.description, encoding="utf-8")
      l = self.__strlen(ba_tmp)
      if l > self.__EDFLIB_WRITE_MAX_ANNOTATION_LEN:
        l = self.__EDFLIB_WRITE_MAX_ANNOTATION_LEN
      str_[p : p + l] = ba_tmp[0 : l]
      p += l
      str_[p] = 20
      p += 1

      for p in range(p, self.__EDFLIB_ANNOTATION_BYTES):
        str_[p] = 0

      self.__file_out.write(str_)

      j += 1
      if j >= self.__nr_annot_chns:
        j = 0
        offset += datrecsize
        datrecs += 1
        if datrecs >= self.__datarecords:
          break

    return 0

# minimum is the minimum digits that will be printed (minus sign not included), leading zero's will be added if necessary
# if sign is zero, only negative numbers will have the sign '-' character
# if sign is one, the sign '+' or '-' character will always be printed
# returns the number of characters printed
  def __sprint_int_number_nonlocalized(self, str_, q, minimum, sign):
    flag = 0
    z = 0
    i = 0
    j = 0
    base = 1000000000

    if minimum < 0:
      minimum = 0

    if minimum > 9:
      flag = 1

    if q < 0:
      str_[j] = 45
      j += 1
      base = -base
    else:
      if sign != 0:
        str_[j] = 43
        j += 1

    for i in range(10, 0, -1):
      if minimum == i:
        flag = 1

      z = q // base

      q = int(q % base)

      if (z != 0) or (flag != 0):
        str_[j] = 48 + z
        j += 1
        flag = 1

      base //= 10

    if flag == 0:
      str_[j] = 48
      j += 1

    str_[j] = 0

    return j

  def __sprint_number_nonlocalized(self, dest, val):
    flag = 0
    z = 0
    i = 0
    j = 0
    base = 1000000000

    sz = len(dest)

    if sz < 1:
      return 0

    q = int(val)

    var = val - q

    if val < 0.0:
      dest[j] = 45
      j += 1

      if q < 0:
        base = -base

    if j == sz:
      j -= 1
      dest[j] = 0
      return j

    for i in range(10, 0, -1):
      z = q // base

      q = int(q % base)

      if (z != 0) or (flag != 0):
        dest[j] = 48 + z
        j += 1

        if j == sz:
          j -= 1
          dest[j] = 0
          return j

        flag = 1

      base //= 10

    if flag == 0:
      dest[j] = 48
      j += 1

    if j == sz:
      j -= 1
      dest[j] = 0
      return j

    base = 100000000

    var *= (base * 10)

    q = int(var)

    if q < 0:
      base = -base

    if q == 0:
      dest[j] = 0
      return j

    dest[j] = 46
    j += 1

    if j == sz:
      j -= 1
      dest[j] = 0
      return j

    for i in range(9, 0, -1):
      z = q // base

      q = int(q % base)

      dest[j] = 48 + z
      j += 1

      if j == sz:
        j -= 1
        dest[j] = 0
        return j

      base //= 10

    dest[j] = 0

    j -= 1

    for j in range(j, 0, -1):
      if dest[j] == 48:
        dest[j] = 0
      else:
        j += 1
        break

    return j

# minimum is the minimum digits that will be printed (minus sign not included), leading zero's will be added if necessary
# if sign is zero, only negative numbers will have the sign '-' character
# if sign is one, the sign '+' or '-' character will always be printed
# returns the number of characters printed
  def __snprint_ll_number_nonlocalized(self, dest, offset, q, minimum, sign):

    flag = 0
    z = 0
    i = 0
    j = offset
    sz = 0
    base = 1000000000000000000

    sz = len(dest)

    if (sz - offset) < 1:
      return 0

    if minimum < 0:
      minimum = 0

    if minimum > 18:
      flag = 1

    if q < 0:
      dest[j] = 45
      j += 1

      base = -base
    else:
      if sign != 0:
        dest[j] = 43
        j += 1

    if j == sz:
      j -= 1
      dest[j] = 0
      return (j - offset)

    for i in range(19, 0, -1):
      if minimum == i:
        flag = 1

      z = q // base

      q = int(q % base)

      if (z != 0) or (flag != 0):
        dest[j] = 48 + z
        j += 1

        if j == sz:
          dest[j] = 0
          j -= 1
          return (j - offset)

        flag = 1

      base = base // 10

    if flag == 0:
      dest[j] = 48
      j += 1

    if j == sz:
      dest[j] = 0
      j -= 1
      return (j - offset)

    dest[j] = 0
    return (j - offset)

# minimum is the minimum digits that will be printed (minus sign not included), leading zero's will be added if necessary
# if sign is zero, only negative numbers will have the sign '-' character
# if sign is one, the sign '+' or '-' character will always be printed
# returns the amount of characters printed
  def __fprint_int_number_nonlocalized(self, f, q, minimum, sign):
    flag = 0
    z = 0
    i = 0
    j = 0
    base = 1000000000

    if minimum < 0:
      minimum = 0

    if minimum > 9:
      flag = 1

    if q < 0:
      f.write(bytes("-", encoding="ascii"))
      j += 1
      base = -base
    else:
      if sign != 0:
        f.write(bytes("+", encoding="ascii"))
        j += 1

    for i in range(10, 0, -1):
      if minimum == i:
        flag = 1

      z = q // base

      q = int(q % base)

      if (z != 0) or (flag != 0):
        f.write(bytes(chr(48 + z), encoding="ascii"))
        j += 1
        flag = 1

      base = base // 10

    if flag == 0:
      f.write(bytes("0", encoding="ascii"))
      j += 1

    return j

# get string length
  def __strlen(self, str_):
    l = len(str_)
    for i in range(0, l):
      if str_[i] == 0:
        return i
    return (i + 1)

# copy a string
  def __strcpy(self, dest, src):
    sz = len(dest) - 1

    srclen = self.__strlen(src)

    if srclen > sz:
      srclen = sz

    if srclen < 0:
      return 0

    for i in range(0, srclen):
      dest[i] = src[i]

    dest[srclen] = 0

    return srclen

# converts Latin-1 to ASCII
  def __latin1_to_ascii(self, str_, l):

    i = 0
    value = 0

    if l > len(str_):
      l = len(str_)

    conv_table = bytearray(" E ,F\".++^mS<E Z  `\'\"\".--~ s>e zY           <        \'u     >   ?AAAAAAECEEEEIIIIDNOOOOOxOUUUUYIsaaaaaaeceeeeiiiidnooooo-0uuuuyty", encoding="ascii")

    for i in range(0, l):
      value = str_[i]

      if (value > 31) and (value < 127):
        continue
      if value < 0:
        value += 256
      if value < 32:
        str_[i] = 32
        continue
      str_[i] = conv_table[value - 127]

################################################################################
# END class EDFwriter
################################################################################

################################################################################
# START class EDFexception
################################################################################

class EDFexception(Exception):
  def __init__(self, message):
    self.message = message
    super().__init__(self.message)

################################################################################
# END class EDFexception
################################################################################




















