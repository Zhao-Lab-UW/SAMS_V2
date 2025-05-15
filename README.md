# SAMS (Semi-Automated MEA Spike sorting )

## Features
SAMS is a graphical user interface application built using MATLAB App Designer that provides:
- Automated spike sorting capabilities
- Support for single electrode and network burst analysis
- Output visualization files 

SAMS_ManualSpikeSorting is a complementary tool to the SAMS automated spike sorting pipeline. It allows users to:
- Visualize and review sorting results from SAMS
- Manually correct any mis-sorted electrodes
- Generate updated sorting results

## Important Notice

This application utilizes file loading functions from the [Axion File Loader](https://github.com/axionbio/AxionFileLoader). Users should check the Axion BioSystems GitHub repository for the most up-to-date loaders to ensure proper data loading functionality.

## Dependencies

This application depends on:
- [Axion File Loader](https://github.com/axionbio/AxionFileLoader) for reading Axion BioSystems data files
- MATLAB Runtime R2023b

## Credits and Acknowledgments

- File loading functionality provided by [Axion BioSystems](https://github.com/axionbio/AxionFileLoader)
- Special thanks to Axion BioSystems for allowing the use of their file loader functions in this open-source tool

## License

This software is provided free for academic and research use. The Axion File Loader components are subject to their original license terms - please refer to [Axion BioSystems' GitHub repository](https://github.com/axionbio/AxionFileLoader) for the most current information.



## Installation Requirements

1. **MATLAB Runtime R2023b (Required)**
   - Download MATLAB Runtime R2023b (free) from: [MATLAB Runtime Download Page](https://www.mathworks.com/products/compiler/matlab-runtime.html)
   - Choose version R2023b (9.15)
   - Select your operating system (Windows 64-bit)
   - Install MATLAB Runtime before running SAMS

2. **System Requirements**
   - Windows 10 or later
   - Minimum 16 GB RAM recommended
   - 8 GB free disk space for MATLAB Runtime

## Installation Steps

1. **Install MATLAB Runtime first:**
   - Download MATLAB Runtime R2023b
   - Run the installer
   - Follow installation prompts
   - Wait for installation to complete (may take 10-15 minutes)

2. **Install SAMS:**
   - Download and extract`SAMS` folder
   - Place it in your preferred location
   - Double-click to run

## Prerequisites

- MATLAB Runtime R2023b installed
- To use SAMS_ManualSpikeSorting, Output files from SAMS pipeline is needed:
  - spike_sorting.xlsx
  - ***.spk files
  - burst_info_all.mat


## Usage Guide

## Usage

### SAMS GUI - Automated Sorting
![automatic sorting GUI](https://github.com/user-attachments/assets/aeb531ad-87dd-4539-baa6-e2efbf659345)

1. **File Selection**
   - Click the "file selection" button
   - Select folder containing .spk files

2. **Parameter Adjustment**
   - Set spike sorting parameters:
     ```
     Threshold for merge: 1.5
     Refractory period (ms): 1.5
     Cutoff frequency (Hz): 0.1
     Flag threshold: 0.2
     Standard deviation cutoff: 3
     ```
   - Configure burst parameters:
     ```
     Single Electrode Bursts:
     - Min spikes: 5
     - Max ISI (ms): 100

     Network Bursts:
     - Min spikes: 50
     - Max ISI (ms): 100
     - Min electrodes participating (%): 35
     ```
3. ** Well Selection (optional) **
   - Click on well selection button to select well(s) to be processed, if no well is selected, all avaliable wells would be processed
4. **Processing**
   - Click "Run"
   - Progress bar shows processing status

### Manual Sorting Interface
<img width="542" alt="manual sorting GUI" src="https://github.com/user-attachments/assets/5f2d3b7f-caff-48cb-8b15-2928e41123d4" />


1. **Load Files**
   - Select "spike_sorting.xlsx" using excel file selection button
   - Select .spk file using spk file selection button
   - Select "burst_info_all.mat" using mat file selection button
<img width="671" alt="output files " src="https://github.com/user-attachments/assets/6dac49bd-5ad6-4864-8b8e-d000081552bf" />

2. **Configure and Load**
   - Adjust parameters if needed
   - Click "apply" to load files

3. **Electrode Processing**
   - Select electrode from dropdown menu
   - Click "select waveforms" for manual sorting
   - Use "add unit" or "exclude unit" buttons as needed
<img width="907" alt="add unit steps" src="https://github.com/user-attachments/assets/8f7994b5-163e-4dc0-bde8-af5a6e9da485" />


4. **Export Results**
   - Click "export" when finished
   - Results saved in "adjusted sorting results" subfolder

## Interface Components

### Main Window
- **Top Panel**: File selection and parameter adjustment
- **Middle Panel**: Burst parameter configuration
- **Bottom Panel**: Progress monitoring

### Manual Sorting Window
- **Upper Plot**: Waveform visualization
- **Lower Plot**: PCA clustering display
- **Right Panel**: Control buttons for sorting operations

## Troubleshooting

1. **Application Won't Start**
   - Verify MATLAB Runtime R2023b is installed
   - Try reinstalling MATLAB Runtime
   - Run as administrator if needed
2. **Files Not Found
   - Verify "SAMS1.pptx" is in the same folder with SAMS app
