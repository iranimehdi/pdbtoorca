This Python code sets up multiscale calculations for the ORCA program. As described in section 8.13 of the ORCA manual, in the QM/MM calculations the whole system is separated into two regions, the QM region and the MM region. The MM region is further partitioned into two regions, the active and non-active regions (cf. Figure 8.38 in the ORCA manual). The active atoms in the MM region are optimized during the geometry optimizations, while the non-active atoms are also included in the geometry optimization, but their positions are constrained. The MM region, which is usually a sphere around the QM region, can be defined either directly in the ORCA input file or via the B-factor column of a PDB file. The MM region sometimes includes hundreds of atoms and makes it tedious and error-prone to define by hand. Our program does this automatically by asking the user for his/her desired thickness of the active region around the QM system (the atoms of the QM region have to be specified in a text file). It also optionally includes either only atoms that are within a certain distance from the QM system, or the groups that have at least one atom within that distance. It also optionally specifies the QM and active regions either in the ORCA input file or in a PDB file. The program also asks the user for his/her desired option to specify the non-active atoms in the MM region: Option 1 - Do not use extended active region Option 2 - Add only atoms bonded covalently to active atoms Option 3 - (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) Option 4 - Manually (using our code) define the extended active region. During the setup steps, the user is also asked to specify the level of theory for the QM part of the QM/MM calculations, the name of the topology file, the number of CPU cores for parallel calculation, the charge of the QM system, and the multiplicity of the QM system. Finally, the program generates an ORCA input file (qmmm.inp by default or an extra PDB file that the QM and MM regions are defined in its occupancy and B-factor columns) to perform the QM/MM calculations. A similar process can also be done to set up QM1/QM2 and QM1/QM2/MM calculations by our program using the qmqm2 and qmqm2mm commands.

Thank you for providing the information about the program and the installation procedure for Linux systems. Let's proceed with adding the installation procedure for Windows users to the README file. Here's how we can update the README:

---

**Installation**

**On Linux:**

1. **Prerequisites:**
   - Make sure you have Python 3.10 installed on your system. You can verify this by running `python3.10 --version` in your terminal. If not installed, download and install it from the official Python website ([Python Downloads](https://www.python.org/downloads/)).

2. **Download the program:**
   - Download the `pdbtoorca.tar.gz` file from the [GitHub repository](https://github.com/iranimehdi/pdbtoorca).

3. **Extract the program files:**
   - Open your terminal.
   - Navigate to the directory containing the downloaded `pdbtoorca.tar.gz` file.
   - Use the following command to extract the contents:
     ```
     tar -xzvf pdbtoorca.tar.gz
     ```

4. **Copy the program files to /usr/bin/:**
   - Navigate into the extracted `pdbtoorca` directory:
     ```
     cd pdbtoorca
     ```
   - Use sudo to copy all contents of the directory to the `/usr/bin/` directory:
     ```
     sudo cp * /usr/bin/
     ```

5. **Verify installation:**
   - You can verify that the program has been installed correctly by running:
     ```
     pdbtoorca 
     ```
Follow the on-screen instructions to convert PDB files and set up multiscale calculations using the ORCA program.

**On Windows:**

1. **Download the program:**
   - Download the standalone Windows binary version (`pdbtoorca.exe`) of the program from the [GitHub releases page](https://github.com/iranimehdi/pdbtoorca/).

2. **Run the program:**
   - After downloading the `pdbtoorca.exe` file, simply double-click the file to run the program on your Windows system.
Follow the on-screen instructions to convert PDB files and set up multiscale calculations using the ORCA program.

**Note:** Ensure that the `pdbtoorca.exe` file is placed in a directory where you have write permissions, as it may generate output files during the conversion process. Please also note that the performance of this version may be slower than its Linux counterpart when run on the same hardware. The executable file can be easily integrated into your workflow on Windows systems, providing a convenient option for users who prefer this platform.
