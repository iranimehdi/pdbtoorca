# PDBtoORCA Toolkit

## Overview

**PDBtoORCA** is a comprehensive Python-based toolkit designed for setting up multiscale quantum mechanics/molecular mechanics (QM/MM), QM1/QM2, and QM1/QM2/MM calculations using the **ORCA** program. It also provides a rich set of utilities to analyze and format **PDB** files.

As described in section 8.13 of the ORCA manual, the QM/MM setup divides the molecular system into two main regions: the **QM region** and the **MM region**. The MM region is further divided into **active** and **non-active** subregions. Active atoms in the MM region are optimized during geometry optimization, while non-active atoms are constrained but still included.

Our program automates the tedious and error-prone task of MM region setup. By specifying the QM atoms (via a text file) and a distance-based thickness for the active region, the tool can identify nearby atoms or entire residue groups. Users can then choose to define the regions either in the ORCA input file or in a PDB file using the occupancy and B-factor columns.

### MM Region Definition Options

The user is prompted to choose how to define the **non-active** region:

1. **Do not use an extended active region**
2. **Add only atoms bonded covalently to active atoms**
3. **(ORCA default)** Use a distance-based criterion
4. **Custom definition** using our toolkit

During setup, the user will also be prompted to enter:
- The level of theory for the QM part
- Topology file name
- Number of CPU cores
- Charge and multiplicity of the QM region

The program then generates either:
- An ORCA input file (`qmmm.inp` by default), or
- A PDB file with QM and MM region encoding in occupancy and B-factor fields.

Similar procedures are available for setting up **QM1/QM2** and **QM1/QM2/MM** calculations using the `qmqm2` and `qmqm2mm` commands.

## Features

- ✅ **PDB File Cleanup**
- ✅ **RMSD Calculation**
- ✅ **Total Charge Calculation**
- ✅ **Chain Extraction**
- ✅ **Geometric Tools** (distance, angle, torsion)
- ✅ **TER Record Insertion**
- ✅ **Protein Residue Renumbering**
- ✅ **Stereochemistry Tools** (Ile/Thr Cβ inversion)

## Installation

### On Linux

1. **Prerequisites:**
   - Python 3.10 is required. Verify with:
     ```bash
     python3.10 --version
     ```
   - If not installed, download from [Python Downloads](https://www.python.org/downloads/).

2. **Download the program:**
   - Download `pdbtoorca.tar.gz` from the [GitHub repository](https://github.com/iranimehdi/pdbtoorca).

3. **Extract the archive:**
   ```bash
   tar -xzvf pdbtoorca.tar.gz
   cd pdbtoorca
   ```

4. **Install to /usr/bin/:**
   ```bash
   sudo cp * /usr/bin/
   ```

5. **Run the tool:**
   ```bash
   pdbtoorca
   ```

### On Windows

1. **Download the EXE:**
   - Grab `pdbtoorca.exe` from the [GitHub releases page](https://github.com/iranimehdi/pdbtoorca/).

2. **Run the tool:**
   - Double-click the EXE or run it from Command Prompt.

> Ensure the EXE is in a directory where you have write permissions, as it will generate output files.

### On macOS

1. **Download the Mac executable:**
   - Available from the same GitHub repository.

2. **Run the program:**
   - Navigate to its directory in Terminal and run:
     ```bash
     ./pdbtoorca
     ```

---

**Citation:**  
If you use this tool in your research, please cite:
> Haji Dehabadi, M.; Saidi, H.; Zafari, F.; Irani, M. *Sci. Rep.* **2024**, *14*(1), 16791.  
> DOI: [10.1038/s41598-024-67468-x](https://doi.org/10.1038/s41598-024-67468-x)

**Contact:**  
Dr. Mehdi Irani — [m.irani@uok.ac.ir](mailto:m.irani@uok.ac.ir)  
Website: [https://prof.uok.ac.ir/m.irani/index.htm](https://prof.uok.ac.ir/m.irani/index.htm)
