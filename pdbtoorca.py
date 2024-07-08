## coding=utf8
# coding: iso-8859-1 -*-
## -*- encoding: utf-8 -*-
# #!/usr/bin/python3.10
# Starting Time
import os
from datetime import datetime
import math
#import numpy as np

start_time = datetime.now()
print("This Python code performs PDB structure analysis and prepares multiscale calculations for the ORCA program.")
print("It was written by Mehdi and Maryam.")
print("The program is free to use. However, if you use it in your project, please cite it as follows:")
print("DOI: ................")

# Open the PDB file
pdb_file_1 = input("Enter the name of your PDB file: ")
pdb_text_1 = open(pdb_file_1, "r")

# Read the contents of the PDB file into a string
data_pdb1 = pdb_text_1.read()

# Count the occurrences of ATOM and HETATM in the PDB file
atom_number_1 = 0
for line in data_pdb1.split('\n'):
    if line.startswith("ATOM") or line.startswith("HETATM"):
        atom_number_1 += 1
print('Number of atoms: ' + "%5s" % str(atom_number_1))

# Count the number of water molecules
WAT_Count = 0
for line in data_pdb1.split('\n'):
    if line[13:20] == "O   WAT" or line[13:20] == "O   HOH":
        WAT_Count += 1

print(f'Number of water molecules: ' + "%5s" % str(WAT_Count))

ca_count = 0
for line in data_pdb1.split('\n'):
    if line.startswith("ATOM") and line[13:15] == "CA":
        ca_count += 1
print(f'Number of amino acids: ' + "%5s" % str(ca_count))
# Commands
while True:
    print("Press Enter to view the available commands.")
    command = input("Enter a command: ").lower()
    if command == "":
        print("orcapdb  - Convert your PDB file to a format suitable for ORCA (with all values in the occupancy and B-factor columns set to 0.00)")
        print("qmmm     - Set up QM/MM calculations for ORCA")
        print("qmqm2    - Set up QM1/QM2 calculations for ORCA")
        print("qmqm2mm  - Set up QM1/QM2/MM calculations for ORCA")
        print("tidy     - Clean up the protein (keep only lines that start with ATOM, HETATM, or TER)")
        print("rmsd     - Calculate RMSD of two PDB files")
        print("charge   - Calculate the total charge of the protein (for PDB files that have partial atomic charges in the last column)")
        print("chain    - Extract a desired chain from the PDB file")
        print("dist     - Calculate the distance between two atoms")
        print("ang      - Calculate the angle between three atoms")
        print("tor      - Calculate the torsion angle between four atoms")
        print("ter      - Insert TER between groups that are not covalently bonded (e.g., solvent molecules)")
        print("terp     - Insert TER in appropriate places of the PDB file that represents a protein system (e.g., between protein chains, after a line with an OXT atom, if there is a missing residue inside a chain, and between hetero groups)")
        print("renumber - Renumber protein residues)")
        print("q        - Quit")
    elif command == "*" or command == "q":
        print("The program has been stopped at your request.")
        break
    elif command.lower() == "tidy" or command.lower() == "t":
        out_pdb = input("Enter the name of your output PDB file (overwrite): ")
        if out_pdb == "":
            out_pdb = pdb_file_1
        else:
            out_pdb = out_pdb

        with open(pdb_file_1, "r") as pdb_text_1:
            data_pdb1 = pdb_text_1.read()

        with open(out_pdb, "w") as output:
            for line in data_pdb1.split("\n"):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    output.write(line + "\n")

        print(f"Your tidied-up PDB file ({out_pdb}) has been generated.")
    elif command.lower() in ["chain", "ch"]:
        # Prompt for chain ID
        chain_id = input("Which chain of the structure do you want to extract? ").strip()
        if not chain_id:
            chain_id = "A"
        # Prompt for output PDB file name
        out_pdb = input("Enter the name of your output PDB file (overwrite): ").strip()
        if not out_pdb:
            out_pdb = pdb_file_1
        # Read the input PDB file and write the specified chain to the output PDB file
        with open(pdb_file_1, "r") as pdb_text_1, open(out_pdb, "w") as output:
            for line in pdb_text_1:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == chain_id:
                    output.write(line)
        print(f"Your PDB file with chain {chain_id} has been generated.")
    elif command.lower() == "renumber" or command.lower() == "re":
        i_number = int(input("Enter the residue number of the first residue in the input file (1): "))
        if i_number == "":
            i_number = int(1)
        else:
            i_number = i_number
        f_number = int(input("Enter the residue number of the first residue in the output file: "))
        out_pdb = input("Enter the name of your output PDB file (overwrite): ")
        if out_pdb == "":
            out_pdb = pdb_file_1
        else:
            out_pdb = out_pdb
        # pdb_text_1 = open(pdb_file_1, "r")
        # Read the original PDB file
        with open(pdb_file_1, 'r') as f:
            lines = f.readlines()
        # Modify residue numbers
        new_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                old_resid = int(line[22:26])
                new_resid = old_resid + int(f_number - i_number)  # Adjust residue numbers by adding 8
                new_line = f"{line[:22]}{new_resid:4}{line[26:]}"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        # Write modified PDB file
        with open(out_pdb, 'w') as f:
            f.writelines(new_lines)
        print(f"The {out_pdb} PDB file with the modified residue number is created.")
    elif command.lower() == "orcapdb" or command.lower() == "o":
        out_pdb = input("Enter the name of your output PDB file (overwrite): ")
        if out_pdb == "":
            out_pdb = pdb_file_1
        else:
            out_pdb = out_pdb

        with open(pdb_file_1, "r") as pdb_text_1:
            data_pdb1 = pdb_text_1.read()

        with open(out_pdb, "w") as output:
            for line in data_pdb1.split("\n"):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    modified_line = line[0:55] + "  0.00  0.00          " + line[76:] + "\n"
                    output.write(modified_line)
            print(f"Your ORCA-compatible PDB file ({out_pdb}) has been generated.")
    elif command.lower() == "rmsd":
        second_pdb = input("Enter the name of the second PDB file: ")
        second_pdb_text = open(second_pdb, "r")

        data_pdb2 = second_pdb_text.read()

        second_n_atoms = 0
        for line in data_pdb2.split('\n'):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                second_n_atoms += 1

        if atom_number_1 != second_n_atoms:
            print(f"The number of atoms in '{pdb_file_1}' and '{second_pdb}' PDB files does not match.")
            break
        else:
            print(f"The number of atoms in '{pdb_file_1}' and '{second_pdb}' PDB files is the same.")
            print("Now checking whether the atom names and symbols in the two PDB files match.")

            out_1_pdb = "first_tided_up_pdb.tmp"
            with open(pdb_file_1, "r") as pdb_text_1:
                data_pdb1 = pdb_text_1.read()
            with open(out_1_pdb, "w") as out_1_pdb_file:
                for line in data_pdb1.split("\n"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        out_1_pdb_file.write(line + "\n")
                print("The first PDB file has been tidied up and stored in a temporary file.")
                print("This file will be deleted later.")

            out_2_pdb = "second_tided_up_pdb.tmp"
            with open(second_pdb, "r") as second_pdb_text:
                data_pdb2 = second_pdb_text.read()

            with open(out_2_pdb, "w") as out_2_pdb_file:
                for line in data_pdb2.split("\n"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        out_2_pdb_file.write(line + "\n")
                print("The second PDB file has been tidied up and stored in a temporary file.")
                print("This file will be deleted later.")

            file_1 = "first_tided_up_pdb.tmp"
            with open(file_1, "r") as text_1:
                data_1 = text_1.read()

            file_2 = "second_tided_up_pdb.tmp"
            with open(file_2, "r") as text_2:
                data_2 = text_2.read()

            boolean = True
            for ind, line_a in enumerate(data_1.split("\n")):
                a = line_a[13:16]
                line_b = data_2.split("\n")[ind]
                b = line_b[13:16]
                if a != b:
                    print(f"Some atoms do not match in the {pdb_file_1} and {second_pdb} PDB files.")
                    boolean = False
                    break

            if boolean:
                print(f"All atoms match in the {pdb_file_1} and {second_pdb} PDB files.")
                print("Now calculating the RMSD value.")

                x_rsd, y_rsd, z_rsd = 0.00, 0.00, 0.00
                for ind, line_a in enumerate(data_1.split("\n")):
                    if line_a[30:38] == "":
                        xi, yi, zi = 0.00, 0.00, 0.00
                    else:
                        xi = float(line_a[30:38])
                        yi = float(line_a[38:46])
                        zi = float(line_a[46:54])
                    line_b = data_2.split("\n")[ind]
                    if line_b[30:38] == "":
                        xj, yj, zj = 0.00, 0.00, 0.00
                    else:
                        xj = float(line_b[30:38])
                        yj = float(line_b[38:46])
                        zj = float(line_b[46:54])
                    x_rsd += (xi - xj) ** 2
                    y_rsd += (yi - yj) ** 2
                    z_rsd += (zi - zj) ** 2
                final_rmsd = ((x_rsd + y_rsd + z_rsd) / atom_number_1) ** 0.5
                r_rmsd = round(final_rmsd, 2)
                print(f"The RMSD value between the {pdb_file_1} and {second_pdb} files is {r_rmsd} Å.")
            text_1.close()
            text_2.close()
            os.remove('first_tided_up_pdb.tmp')
            os.remove('second_tided_up_pdb.tmp')
    elif command.lower() == "dist" or command.lower() == "d":
        first_atom_n = input("Enter the atom number of the first atom: ")
        seco_atom_n = input("Enter the atom number of the second atom: ")

        xd1, yd1, zd1 = 0.00, 0.00, 0.00
        xd2, yd2, zd2 = 0.00, 0.00, 0.00
        distance = 0.00
        for line_d1 in data_pdb1.split('\n'):
            if line_d1.startswith("ATOM") or line_d1.startswith("HETATM"):
                if int(line_d1[6:11]) == int(first_atom_n):
                    xd1 = float(line_d1[30:38])
                    yd1 = float(line_d1[38:46])
                    zd1 = float(line_d1[46:54])
                if int(line_d1[6:11]) == int(seco_atom_n):
                    xd2 = float(line_d1[30:38])
                    yd2 = float(line_d1[38:46])
                    zd2 = float(line_d1[46:54])
                distance = round(((xd1 - xd2) ** 2 + (yd1 - yd2) ** 2 + (zd1 - zd2) ** 2) ** 0.5, 6)

        print(f"The distance between atom {first_atom_n} and atom {seco_atom_n} is {distance} angstrom.")
    elif command.lower() == "ang":
        atom_i = input("Enter the atom number of the first            atom: ")
        atom_j = input("Enter the atom number of the second (central) atom: ")
        atom_k = input("Enter the atom number of the third            atom: ")

        xi, yi, zi = 0.00, 0.00, 0.00
        xj, yj, zj = 0.00, 0.00, 0.00
        xk, yk, zk = 0.00, 0.00, 0.00
        abs_r_ji, abs_r_jk = 0.00, 0.00
        x_ji, y_ji, z_ji = 0.00, 0.00, 0.00
        x_jk, y_jk, z_jk = 0.00, 0.00, 0.00
        angle_ijk = 0.00

        for line_d1 in data_pdb1.split('\n'):
            if line_d1.startswith("ATOM") or line_d1.startswith("HETATM"):
                if int(line_d1[6:11]) == int(atom_i):
                    xi = float(line_d1[30:38])
                    yi = float(line_d1[38:46])
                    zi = float(line_d1[46:54])

        for line_d1 in data_pdb1.split('\n'):
            if line_d1.startswith("ATOM") or line_d1.startswith("HETATM"):
                if int(line_d1[6:11]) == int(atom_j):
                    xj = float(line_d1[30:38])
                    yj = float(line_d1[38:46])
                    zj = float(line_d1[46:54])

        for line_d1 in data_pdb1.split('\n'):
            if line_d1.startswith("ATOM") or line_d1.startswith("HETATM"):
                if int(line_d1[6:11]) == int(atom_k):
                    xk = float(line_d1[30:38])
                    yk = float(line_d1[38:46])
                    zk = float(line_d1[46:54])

        abs_r_ji = ((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2) ** 0.5
        abs_r_jk = ((xk - xj) ** 2 + (yk - yj) ** 2 + (zk - zj) ** 2) ** 0.5

        x_ji, y_ji, z_ji = xi - xj, yi - yj, zi - zj
        x_jk, y_jk, z_jk = xk - xj, yk - yj, zk - zj

        angle_ijk = round(
            math.acos(((1 / (abs_r_jk * abs_r_ji)) * (x_ji * x_jk + y_ji * y_jk + z_ji * z_jk))) * 180 / math.pi, 4)

        print(f"The angle between atoms {atom_i}, {atom_j}, {atom_k} is '{angle_ijk}' degrees.")
    elif command.lower() == "tor":
        print("We are still working on this command")
    ###########################################################################################
        # Preparation of ORCA input files for QM/MM calculations
    elif command.lower() == "qmmm" or command.lower() == "qm" or command.lower() == "qmm":
        pdb_or_inp = input(f"Do you want to specify the QM and active regions in an ORCA input file (o) or in a PDB file (p)? ")
        if pdb_or_inp.lower() == "o" or pdb_or_inp == "":
            with open("orcapdb.pdb", "w") as output:
                for line in data_pdb1.split("\n"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        rline = line[0:55] + "  0.00  0.00          " + line[76:] + "\n"
                        output.write(rline)
                #output.write("END   00000          00000       0.000  0.0000  0.0000   0.00  0.00")

            # Read contents of the temp.pdb file to a string
            pdb_file_orca = 'orcapdb.pdb'
            pdb_text_1 = open("orcapdb.pdb", "r")
            data_pdb1 = pdb_text_1.read()
            # Ask the user for the S1 file name, the file describing the QM region
            s1_file = input("Please provide the file name that contains the QM region description (s1): ")
            if s1_file == "":
                s1_file = "s1"
            else:
                s1_file == s1_file
            # Open the S1 file, read its lines and connect the lines
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                linear_string = str(",".join(s1_line.strip() for s1_line in s1_lines_No_Blank_lines))
                # Convert the numeric string ranges to a list (linear_qm_atoms_list). E.g., convert 5-7 to 5, 6, 7
                def f(linear_string):
                    result = []
                    for part in linear_string.split(','):
                        if '-' in part:
                            a, b = part.split('-')
                            a, b = int(a), int(b)
                            result.extend(range(a, b + 1))
                        else:
                            a = int(part)
                            result.append(a)
                    return result
                Linear_QM_Atoms_List_PDB_FORMAT = (f(linear_string))
                Linear_QM_Atoms_List_PDB_FORMAT_sorted = sorted(Linear_QM_Atoms_List_PDB_FORMAT)
            #print(f"QM atoms list in PDB format: {Linear_QM_Atoms_List_PDB_FORMAT_sorted}")
            ################################################
            # Convert the S1 file to the ORCA input file format e.g., {1 5:10 12}
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                qm_atoms_for_orca_input_file = ""
                for qm_atoms in s1_lines_No_Blank_lines:
                    qm_atoms_orca_format = qm_atoms.strip().split('-')
                    qm_atoms_orca_format = [int(element) - 1 for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = [str(element) for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = ":".join(qm_atoms_orca_format)
                    qm_atoms_for_orca_input_file += qm_atoms_orca_format + " "
            Qm_lst = ['QMatoms ', '{', qm_atoms_for_orca_input_file, '}', ' end']
            QM_atoms_list_ORCA_format = ''.join(map(str, Qm_lst))
            print(f"The following {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are included in the QM region (excluding junctions).")
            print(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
            print("ORCA identifies junctions itself and will include them in the QM region.")
            #########################################################################
            print(
                f"The PDB file containing the QM atoms (QM.pdb) is generated.")            
            pdb_text_11 = open(pdb_file_1, "r")
            data_pdb11 = pdb_text_11.readlines()
            with open ("QM.pdb", "w") as output:
                for Pline in data_pdb11:
                    if int(Pline[6:11].strip()) in Linear_QM_Atoms_List_PDB_FORMAT_sorted:
                        output.write(Pline)           
            #########################################################################    
            # Ask the user for the thickness of the active region
            active_region_thickness = input("Enter the thickness of the active region around the QM region in angstrom (6.0): ")
            if active_region_thickness == "":
                active_region_thickness = float(6.0)
            else:
                active_region_thickness = float(active_region_thickness)
            # Define the variables
            x_qm = 0.00
            y_qm = 0.00
            z_qm = 0.00
            x_mm = 0.00
            y_mm = 0.00
            z_mm = 0.00
            SD = 0.00
            # Define the new list
            atoms_or_groups = input("Do you wish the active region to consist only of atoms that are within the specified distance from the QM region (a), or to consist of groups that have at least one atom within that distance (g)? ")
            if atoms_or_groups.lower() == "g" or atoms_or_groups.lower() == "":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ##
                list_resi_num = []
                list_all_mm = []
                for line in data_pdb1.split('\n')[:-1]:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if int(line[6:11].strip()) in mm_list_sorted:
                            atm_num = line[20:26].strip()
                            list_resi_num.append(atm_num)
                            for line in data_pdb1.split('\n')[:-1]:
                                if line[20:26].strip() == atm_num:
                                    list_all_mm.append(line[6:11].strip())
                                else:
                                    continue
                        else:
                            continue
                    else:
                        continue
                new_list_all_mm = []
                for one_mm in list_all_mm:
                    if int(one_mm) not in new_list_all_mm:
                        new_list_all_mm.append(int(one_mm))
                    else:
                        continue
                ## simplify the list
                # mm_list_sortedf = mm_list_sorted[:-1]
                # mmlistoc = []
                # for i in range(0,len(mm_list_sortedf)):
                # if i <= 0:
                # for j in range (0,len(mm_list_sortedf)):
                # if mm_list_sortedf [i] in range(mm_list_sortedf[i] and mm_list_sortedf[j]+1) and mm_list_sortedf[j]+1 not in mm_list_sortedf and mm_list_sortedf[i]-1 not in mm_list_sortedf:
                # if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                # if mm_list_sortedf[i] != mm_list_sortedf[j]:
                # mmlistoc.append(str(mm_list_sortedf[i]-1)+":"+str(mm_list_sortedf[j]-1))
                # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                # i = j+1
                # else:
                # mmlistoc.append(str(mm_list_sortedf[i]-1))
                # print(str(mm_list_sortedf[i]))
                # i = j+1
                # lsf=str(mmlistoc).replace('[','')
                # lsff=str(lsf).replace(']','')
                # lsfff=str(lsff).replace("'","")
                # lsffff=str(lsfff).replace(",","")
                # mm_lsOrca = 'ActiveAtoms'+' {'+ lsffff+'} '+ 'end'
                # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                ## simplify the new-list_mm
                oc_list_all_mm = new_list_all_mm
                mmlistoc_new = []
                for i in range(0, len(oc_list_all_mm)):
                    if i <= 0:
                        for j in range(0, len(oc_list_all_mm)):
                            if oc_list_all_mm[i] in range(oc_list_all_mm[i] and oc_list_all_mm[j] + 1) and oc_list_all_mm[
                                j] + 1 not in oc_list_all_mm and oc_list_all_mm[i] - 1 not in oc_list_all_mm:
                                if oc_list_all_mm[i] != min(oc_list_all_mm) or oc_list_all_mm[j] != max(oc_list_all_mm):
                                    if oc_list_all_mm[i] != oc_list_all_mm[j]:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1) + ":" + str(oc_list_all_mm[j] - 1))
                                        # print(str(oc_list_all_mm[i])+":"+str(oc_list_all_mm[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1))
                                        # print(str(oc_list_all_mm[i]))
                                        i = j + 1
                        lsf_new = str(mmlistoc_new).replace('[', '')
                        lsff_new = str(lsf_new).replace(']', '')
                        lsfff_new = str(lsff_new).replace("'", "")
                        lsffff_new = str(lsfff_new).replace(",", "")
                        mm_lsOrca_new = 'ActiveAtoms' + ' {' + lsffff_new + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff_new +'} '+ 'end')
                # print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                # print(mm_list_sorted)
                print(f"The following {len(new_list_all_mm)} atoms are included in the active region.")
                print(new_list_all_mm)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            elif atoms_or_groups.lower() == "a":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ## simplify the list
                mm_list_sortedf = mm_list_sorted
                mmlistoc = []
                for i in range(0, len(mm_list_sortedf)):
                    if i <= 0:
                        for j in range(0, len(mm_list_sortedf)):
                            if mm_list_sortedf[i] in range(mm_list_sortedf[i] and mm_list_sortedf[j] + 1) and \
                                    mm_list_sortedf[j] + 1 not in mm_list_sortedf and mm_list_sortedf[
                                i] - 1 not in mm_list_sortedf:
                                if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                                    if mm_list_sortedf[i] != mm_list_sortedf[j]:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1) + ":" + str(mm_list_sortedf[j] - 1))
                                        # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1))
                                        # print(str(mm_list_sortedf[i]))
                                        i = j + 1
                        lsf = str(mmlistoc).replace('[', '')
                        lsff = str(lsf).replace(']', '')
                        lsfff = str(lsff).replace("'", "")
                        lsffff = str(lsfff).replace(",", "")
                        mm_lsOrca = 'ActiveAtoms' + ' {' + lsffff + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                print(mm_list_sorted)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            else:
                print("type 'a' or 'g'")
                break
            #########################################################################
            print(
                f"The PDB file containing the QM and MM active regions (QM-MMactive.pdb) is generated.")            
            pdb_text_12 = open(pdb_file_1, "r")
            data_pdb12 = pdb_text_12.readlines()
            if atoms_or_groups.lower() == "a":
                with open ("QM-MMactive.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in mm_list_sorted:
                            output.write(Pline)
            if atoms_or_groups.lower() == "g":
                with open ("QM-MMactive.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in new_list_all_mm:
                            output.write(Pline)           
            #########################################################################
            print("Now, we are assigning the extension shell. Be patient, this might take a while.")
            print("ORCA optimizes only the positions of active atoms in its QM/MM geometry optimization algorithm."
                  "Nevertheless, the forces exerted on these active atoms are influenced by their interactions with the"
                  " non-active atoms around them. To obtain a smooth convergence of quasi-Newton algorithms in internal"
                  " coordinates, the Hessian values between the active atoms and the directly surrounding non-active atoms"
                  " must be available. For that reason, the active atoms are extended by a shell of surrounding non-active"
                  " atoms which are also included in the geometry optimization, but whose positions are constrained. This"
                  " shell of non-active atoms can be automatically chosen by ORCA (see section 8.13.1 in the manual). This"
                  " code includes two options for choosing the extension shell. The first is the default option of ORCA"
                  " which includes those non-active atoms in the extension shell that have a distance of less than the sum"
                  " of their VDW radii plus Dist AtomsAroundOpt. The second option is to create a shell of atoms in the"
                  " non-active region at a specific distance from the active region. The user will be asked to provide the"
                  " distance.")
            # Ask the user to choose the atoms belonging to the extension shell
            Ex_sh_options =  input("Select an option to choose the atoms belonging to optRegionExt: " + "\n"
                                   "Option 1       - Do not use extended active region" + "\n"
                                   "Option 2       - Add only atoms bonded covalently to active atoms" + "\n"
                                   "Option 3       - (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1))"+ "\n"
                                   "Option 4       - Manually (using our code) define the extended active region (Default: empty list)" + "\n"
                                   "Enter your option: ")
            # option 1
            if Ex_sh_options == "1":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB false \n")
                #
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB false \n")
                qmmm_file.write(f"# Option 1: no extension shell will be included in the QM/MM calculations \n")
                qmmm_file.write("ExtendActiveRegion No \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_orca} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 2
            if Ex_sh_options == "2":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB false \n")
                #
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB false \n")
                qmmm_file.write(f"# Option 2: add only atoms bonded covalently to active atoms \n")
                qmmm_file.write("ExtendActiveRegion cov_bonds \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_orca} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 3
            if Ex_sh_options == "3":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB false \n")
                #
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB false \n")
                qmmm_file.write(f"# Option 3: (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) \n")
                qmmm_file.write("ExtendActiveRegion distance \n")
                distance_criterion = input("Enter the Dist_AtomsAroundOpt value (Default; 1): ")
                if distance_criterion == "":
                    qmmm_file.write("Dist_AtomsAroundOpt 1.\n")
                else:
                    qmmm_file.write(f"Dist_AtomsAroundOpt {distance_criterion}\n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_orca} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 4
            if Ex_sh_options == "4":
                if atoms_or_groups.lower() == "a":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in mm_list_sorted:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    #list_resi_num_es = []
                    #list_all_es = []
                    #for line in data_pdb1.split('\n')[:-1]:
                        #if line.startswith("ATOM") or line.startswith("HETATM"):
                            #if int(line[6:11].strip()) in es_list:
                                #atm_num = line[22:26].strip()
                                #list_resi_num_es.append(atm_num)
                                #for line in data_pdb1.split('\n')[:-1]:
                                    #if line[22:26].strip() == atm_num:
                                        #list_all_es.append(line[6:11].strip())
                                    #else:
                                        #continue
                            #else:
                                #continue
                        #else:
                            #continue
    
                    #new_list_all_es = []
                    #for one_mm in list_all_es:
                        #if int(one_mm) not in new_list_all_es:
                            #new_list_all_es.append(int(one_mm))
                        #else:
                            #continue
                    es_list_new = set(es_list) - set(mm_list_sorted)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                elif atoms_or_groups.lower() == "g":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in new_list_all_mm:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    list_resi_num_es = []
                    list_all_es = []
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(line[6:11].strip()) in es_list:
                                atm_num = line[20:26].strip()
                                list_resi_num_es.append(atm_num)
                                for line in data_pdb1.split('\n')[:-1]:
                                    if line[20:26].strip() == atm_num:
                                        list_all_es.append(line[6:11].strip())
                                    else:
                                        continue
                            else:
                                continue
                        else:
                            continue
    
                    new_list_all_es = []
                    for one_mm in list_all_es:
                        if int(one_mm) not in new_list_all_es:
                            new_list_all_es.append(int(one_mm))
                        else:
                            continue
                    #es_list_new = set(new_list_all_es)                
                    es_list_new = set(new_list_all_es) - set(new_list_all_mm)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                # Assign the output PDB file for QM/MM calculations by ORCA
                #out_pdb_for_ORCA_QMMM = input(
                    #"Enter the name of your output PDB file, to assign the QM atoms and active region in its "
                    #"occupancy and B-factor columns, respectively. (out_pdb_ORCA_QMMM.pdb): ")
                #if out_pdb_for_ORCA_QMMM == "":
                    #out_pdb_for_ORCA_QMMM = "out_pdb_ORCA_QMMM.pdb"
                #else:
                    #out_pdb_for_ORCA_QMMM == out_pdb_for_ORCA_QMMM
                # We expect from the following lines to read the input PDB file (defined in the beginig of the code), line by line.
                # Read the mm_list_sorted list item by item.
                # If indexes 6–11 in a line of the input pdb file are equal to the item, replace indexes 64–67 of that line with 1.00.
                #lower = 0
                #upper = len(mm_list_sorted)
                #low_pdb = 0
                #line_in_input_pdb = data_pdb1.split("\n")[:-1]
                #up_pdb = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
                #with open(out_pdb_for_ORCA_QMMM, "w") as output:
                    #for mm_atom in range(lower, upper):
                        #if mm_atom <= lower:
                            #mm_list_sorted.append("")
                            #for qm_atom in range(low_qm,up_qm):
                                #if qm_atom <= low_qm:
                                    #Linear_QM_Atoms_List_PDB_FORMAT_sorted.append("")
                                    #for line_pdb in range (low_pdb, up_pdb):
                                        #if str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip()and str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + "  1.00"+ line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + line_in_input_pdb[line_pdb][61:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:63] + "1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                        #else:
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb] + "\n"
                                        #output.write(line_in_output_pdb)
                    #output.write("END                                                                ")
                #########################################################################
                print(
                    f"The PDB file containing the QM, MM active,and -fixed  regions (QM-MMactive-fixed.pdb) is generated.")            
                pdb_text_13 = open(pdb_file_1, "r")
                data_pdb13 = pdb_text_13.readlines()
                if atoms_or_groups.lower() == "a":
                    totpdb = sorted(mm_list_sorted + es_list_sorted_PDB_format)
                    with open ("QM-MMactive-fixed.pdb", "w") as output:
                        for Pline in data_pdb13:
                            if int(Pline[6:11].strip()) in totpdb:
                                output.write(Pline)
                if atoms_or_groups.lower() == "g":
                    totpdb = sorted(new_list_all_mm + es_list_sorted_PDB_format)
                    with open ("QM-MMactive-fixed.pdb", "w") as output:
                        for Pline in data_pdb13:
                            if int(Pline[6:11].strip()) in totpdb:
                                output.write(Pline)           
                #########################################################################

                # Write the information into the ORCA input file (qmmm.inp).
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB false \n")
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB false \n")
                #elif pdb_or_inp == "p" or pdb_or_inp == "pdb":
                    #qmmm_file.write(
                        #f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region."
                        #f" They are assigned with '1.00' in the occupancy column of the {out_pdb_for_ORCA_QMMM} file. \n")
                    #qmmm_file.write("Use_QM_InfoFromPDB true \n")
                    #qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                    #f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                    #f"region. They are assigned with '1.00' in the B-factor column of the {out_pdb_for_ORCA_QMMM} file. \n")
                    #qmmm_file.write("Use_Active_InfoFromPDB true \n")
                #E_Shell_List = input(f"Do you want to include Extension Shell in ORCA input (no)  or (yes)? ")
                #if E_Shell_List.lower() == "no" or E_Shell_List.lower() == "n" or E_Shell_List.lower() == "":
                    #print("Non-active atoms that have a distance of less than the sum of their VDW radii plus Dist"
                          #" AtomsAroundOpt will be included in the extension shell. This the default option of ORCA.")
                    #qmmm_file.write(
                    #    f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                    #    f"included in the active region ({len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms in the QM region, "
                    #    f"excluding junctions). \n")
                    #qmmm_file.write(f"# Non-active atoms that have a distance of less than the sum of their VDW radii plus Dist"
                          #" AtomsAroundOpt will be included in the extension shell. This the default option of ORCA. \n")
                #elif E_Shell_List.lower() == "yes" or E_Shell_List.lower() == "y":
                qmmm_file.write(f"# An extension shell of {extension_shell_thickness}-angstrom thickness will surround the "
                                    f"active region. \n")
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(
                        f"# There are {atom_number_1} atoms in your PDB file, from which {len(new_list_all_mm)} atoms are "
                        f"included in the active region and {len(es_list_sorted_PDB_format)} atoms are included in the"
                        f" extension shell. Hence, {atom_number_1 - len(es_list_new) - len(new_list_all_mm)}"
                        f" atoms will be ignored in the calculations. \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(
                        f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                        f"included in the active region and {len(es_list_sorted_PDB_format)} atoms are included in the"
                        f" extension shell. Hence, {atom_number_1 - len(es_list_new) - len(mm_list_sorted)}"
                        f" atoms will be ignored in the calculations. \n")                    
                qmmm_file.write(f"{es_lsOrca_new} \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_orca} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
        elif pdb_or_inp.lower() == "p":
            with open("orcapdb.pdb", "w") as output:
                for line in data_pdb1.split("\n"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        rline = line[0:55] + "  0.00  0.00          " + line[76:] + "\n"
                        output.write(rline)
                #output.write("END   00000          00000       0.000  0.0000  0.0000   0.00  0.00")

            # Read contents of the temp.pdb file to a string
            pdb_file_orca = 'orcapdb.pdb'
            pdb_text_1 = open("orcapdb.pdb", "r")
            data_pdb1 = pdb_text_1.read()
            
            # Ask the user for the S1 file name, the file describing the QM region
            s1_file = input("Please provide the file name that contains the QM region description (s1): ")
            if s1_file == "":
                s1_file = "s1"
            else:
                s1_file == s1_file
            # Open the S1 file, read its lines and connect the lines
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                linear_string = str(",".join(s1_line.strip() for s1_line in s1_lines_No_Blank_lines))
                # Convert the numeric string ranges to a list (linear_qm_atoms_list). E.g., convert 5-7 to 5, 6, 7
                def f(linear_string):
                    result = []
                    for part in linear_string.split(','):
                        if '-' in part:
                            a, b = part.split('-')
                            a, b = int(a), int(b)
                            result.extend(range(a, b + 1))
                        else:
                            a = int(part)
                            result.append(a)
                    return result
                Linear_QM_Atoms_List_PDB_FORMAT = (f(linear_string))
                Linear_QM_Atoms_List_PDB_FORMAT_sorted = sorted(Linear_QM_Atoms_List_PDB_FORMAT)
            #print(f"QM atoms list in PDB format: {Linear_QM_Atoms_List_PDB_FORMAT_sorted}")
            ################################################
            # Convert the S1 file to the ORCA input file format e.g., {1 5:10 12}
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                qm_atoms_for_orca_input_file = ""
                for qm_atoms in s1_lines_No_Blank_lines:
                    qm_atoms_orca_format = qm_atoms.strip().split('-')
                    qm_atoms_orca_format = [int(element) - 1 for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = [str(element) for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = ":".join(qm_atoms_orca_format)
                    qm_atoms_for_orca_input_file += qm_atoms_orca_format + " "
            Qm_lst = ['QMatoms ', '{', qm_atoms_for_orca_input_file, '}', ' end']
            QM_atoms_list_ORCA_format = ''.join(map(str, Qm_lst))
            print(f"The following {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are included in the QM region (excluding junctions).")
            print(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
            print("ORCA identifies junctions itself and will include them in the QM region.")
    
            # Ask the user for the thickness of the active region
            active_region_thickness = input("Enter the thickness of the active region around the QM region in angstrom (6.0): ")
            if active_region_thickness == "":
                active_region_thickness = float(6.0)
            else:
                active_region_thickness = float(active_region_thickness)
            # Define the variables
            x_qm = 0.00
            y_qm = 0.00
            z_qm = 0.00
            x_mm = 0.00
            y_mm = 0.00
            z_mm = 0.00
            SD = 0.00
            # Define the new list
            atoms_or_groups = input("How do you want to define the active region? Do you want to include only atoms that are within a certain distance from the QM region (a), or the groups that have at least one atom within that distance (g)? ")
            if atoms_or_groups.lower() == "g":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ##
                list_resi_num = []
                list_all_mm = []
                for line in data_pdb1.split('\n')[:-1]:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if int(line[6:11].strip()) in mm_list_sorted:
                            atm_num = line[20:26].strip()
                            list_resi_num.append(atm_num)
                            for line in data_pdb1.split('\n')[:-1]:
                                if line[20:26].strip() == atm_num:
                                    list_all_mm.append(line[6:11].strip())
                                else:
                                    continue
                        else:
                            continue
                    else:
                        continue
                new_list_all_mm = []
                for one_mm in list_all_mm:
                    if int(one_mm) not in new_list_all_mm:
                        new_list_all_mm.append(int(one_mm))
                    else:
                        continue
                ## simplify the list
                # mm_list_sortedf = mm_list_sorted[:-1]
                # mmlistoc = []
                # for i in range(0,len(mm_list_sortedf)):
                # if i <= 0:
                # for j in range (0,len(mm_list_sortedf)):
                # if mm_list_sortedf [i] in range(mm_list_sortedf[i] and mm_list_sortedf[j]+1) and mm_list_sortedf[j]+1 not in mm_list_sortedf and mm_list_sortedf[i]-1 not in mm_list_sortedf:
                # if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                # if mm_list_sortedf[i] != mm_list_sortedf[j]:
                # mmlistoc.append(str(mm_list_sortedf[i]-1)+":"+str(mm_list_sortedf[j]-1))
                # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                # i = j+1
                # else:
                # mmlistoc.append(str(mm_list_sortedf[i]-1))
                # print(str(mm_list_sortedf[i]))
                # i = j+1
                # lsf=str(mmlistoc).replace('[','')
                # lsff=str(lsf).replace(']','')
                # lsfff=str(lsff).replace("'","")
                # lsffff=str(lsfff).replace(",","")
                # mm_lsOrca = 'ActiveAtoms'+' {'+ lsffff+'} '+ 'end'
                # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                ## simplify the new-list_mm
                oc_list_all_mm = new_list_all_mm
                mmlistoc_new = []
                for i in range(0, len(oc_list_all_mm)):
                    if i <= 0:
                        for j in range(0, len(oc_list_all_mm)):
                            if oc_list_all_mm[i] in range(oc_list_all_mm[i] and oc_list_all_mm[j] + 1) and oc_list_all_mm[
                                j] + 1 not in oc_list_all_mm and oc_list_all_mm[i] - 1 not in oc_list_all_mm:
                                if oc_list_all_mm[i] != min(oc_list_all_mm) or oc_list_all_mm[j] != max(oc_list_all_mm):
                                    if oc_list_all_mm[i] != oc_list_all_mm[j]:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1) + ":" + str(oc_list_all_mm[j] - 1))
                                        # print(str(oc_list_all_mm[i])+":"+str(oc_list_all_mm[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1))
                                        # print(str(oc_list_all_mm[i]))
                                        i = j + 1
                        lsf_new = str(mmlistoc_new).replace('[', '')
                        lsff_new = str(lsf_new).replace(']', '')
                        lsfff_new = str(lsff_new).replace("'", "")
                        lsffff_new = str(lsfff_new).replace(",", "")
                        mm_lsOrca_new = 'ActiveAtoms' + ' {' + lsffff_new + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff_new +'} '+ 'end')
                # print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                # print(mm_list_sorted)
                print(f"The following {len(new_list_all_mm)} atoms are included in the active region.")
                print(new_list_all_mm)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            elif atoms_or_groups.lower() == "a":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ## simplify the list
                mm_list_sortedf = mm_list_sorted[:-1]
                mmlistoc = []
                for i in range(0, len(mm_list_sortedf)):
                    if i <= 0:
                        for j in range(0, len(mm_list_sortedf)):
                            if mm_list_sortedf[i] in range(mm_list_sortedf[i] and mm_list_sortedf[j] + 1) and \
                                    mm_list_sortedf[j] + 1 not in mm_list_sortedf and mm_list_sortedf[
                                i] - 1 not in mm_list_sortedf:
                                if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                                    if mm_list_sortedf[i] != mm_list_sortedf[j]:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1) + ":" + str(mm_list_sortedf[j] - 1))
                                        # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1))
                                        # print(str(mm_list_sortedf[i]))
                                        i = j + 1
                        lsf = str(mmlistoc).replace('[', '')
                        lsff = str(lsf).replace(']', '')
                        lsfff = str(lsff).replace("'", "")
                        lsffff = str(lsfff).replace(",", "")
                        mm_lsOrca = 'ActiveAtoms' + ' {' + lsffff + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                print(mm_list_sorted)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            else:
                print("type 'a' or 'g'")
                break

            # Assign the output PDB file for QM/MM calculations by ORCA
            if atoms_or_groups.lower() == "a":
                out_pdb_for_ORCA_QMMM = input(
                    "Enter the name of your output PDB file, to assign the QM atoms and active region in its "
                    "occupancy and B-factor columns, respectively. (out_pdb_ORCA_QMMM.pdb): ")
                if out_pdb_for_ORCA_QMMM == "":
                    out_pdb_for_ORCA_QMMM = "out_pdb_ORCA_QMMM.pdb"
                else:
                    out_pdb_for_ORCA_QMMM == out_pdb_for_ORCA_QMMM
                # We expect from the following lines to read the input PDB file (defined in the beginig of the code), line by line.
                # Read the mm_list_sorted list item by item.
                # If indexes 6–11 in a line of the input pdb file are equal to the item, replace indexes 64–67 of that line with 1.00.
                lower = 0
                upper = len(mm_list_sorted)
                low_pdb = 0
                line_in_input_pdb = data_pdb1.split("\n")[:-1]
                up_pdb = len(line_in_input_pdb)
                low_qm = 0
                up_qm = len(line_in_input_pdb)
                low_qm = 0
                up_qm = len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
                with open(out_pdb_for_ORCA_QMMM, "w") as output:
                    for mm_atom in range(lower, upper):
                        if mm_atom <= lower:
                            mm_list_sorted.append("")
                            for qm_atom in range(low_qm,up_qm):
                                if qm_atom <= low_qm:
                                    Linear_QM_Atoms_List_PDB_FORMAT_sorted.append("")
                                    for line_pdb in range (low_pdb, up_pdb):
                                        if line_in_input_pdb[line_pdb].startswith("TER") or line_in_input_pdb[line_pdb].startswith("END"):
                                            line_in_output_pdb = line_in_input_pdb[line_pdb] + "\n" 
                                        elif str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip()and str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "1.00" + "  1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                            qm_atom += 1
                                        elif str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "1.00" + "  0.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                            qm_atom += 1
                                        elif str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "0.00" + "  1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                        else:
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "0.00" + "  0.00"+ line_in_input_pdb[line_pdb][67:] + "\n"
                                        output.write(line_in_output_pdb)
                mm_list_sorted = mm_list_sorted[:-1]
             
            # Assign the output PDB file for QM/MM calculations by ORCA
            elif atoms_or_groups.lower() == "g":
                out_pdb_for_ORCA_QMMM = input(
                    "Enter the name of your output PDB file, to assign the QM atoms and active region in its "
                    "occupancy and B-factor columns, respectively. (out_pdb_ORCA_QMMM.pdb): ")
                if out_pdb_for_ORCA_QMMM == "":
                    out_pdb_for_ORCA_QMMM = "out_pdb_ORCA_QMMM.pdb"
                else:
                    out_pdb_for_ORCA_QMMM == out_pdb_for_ORCA_QMMM
                # We expect from the following lines to read the input PDB file (defined in the beginig of the code), line by line.
                # Read the mm_list_sorted list item by item.
                # If indexes 6–11 in a line of the input pdb file are equal to the item, replace indexes 64–67 of that line with 1.00.
                lower = 0
                upper = len(new_list_all_mm)
                low_pdb = 0
                line_in_input_pdb = data_pdb1.split("\n")[:-1]
                up_pdb = len(line_in_input_pdb)
                low_qm = 0
                up_qm = len(line_in_input_pdb)
                low_qm = 0
                up_qm = len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
                with open(out_pdb_for_ORCA_QMMM, "w") as output:
                    for mm_atom in range(lower, upper):
                        if mm_atom <= lower:
                            new_list_all_mm.append("")
                            for qm_atom in range(low_qm,up_qm):
                                if qm_atom <= low_qm:
                                    Linear_QM_Atoms_List_PDB_FORMAT_sorted.append("")
                                    for line_pdb in range (low_pdb, up_pdb):
                                        if line_in_input_pdb[line_pdb].startswith("TER") or line_in_input_pdb[line_pdb].startswith("END"):
                                            line_in_output_pdb = line_in_input_pdb[line_pdb] + "\n" 
                                        elif str(new_list_all_mm[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip()and str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "1.00" + "  1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                            qm_atom += 1
                                        elif str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "1.00" + "  0.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                            qm_atom += 1
                                        elif str(new_list_all_mm[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "0.00" + "  1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            mm_atom += 1
                                        else:
                                            line_in_output_pdb = line_in_input_pdb[line_pdb][0:56] + "0.00" + "  0.00"+ line_in_input_pdb[line_pdb][67:] + "\n"
                                        output.write(line_in_output_pdb)
        
                
                new_list_all_mm = new_list_all_mm[:-1]
            print("Now, we are assigning the extension shell. Be patient, this might take a while.")
            print("ORCA optimizes only the positions of active atoms in its QM/MM geometry optimization algorithm."
                  "Nevertheless, the forces exerted on these active atoms are influenced by their interactions with the"
                  " non-active atoms around them. To obtain a smooth convergence of quasi-Newton algorithms in internal"
                  " coordinates, the Hessian values between the active atoms and the directly surrounding non-active atoms"
                  " must be available. For that reason, the active atoms are extended by a shell of surrounding non-active"
                  " atoms which are also included in the geometry optimization, but whose positions are constrained. This"
                  " shell of non-active atoms can be automatically chosen by ORCA (see section 8.13.1 in the manual). This"
                  " code includes two options for choosing the extension shell. The first is the default option of ORCA"
                  " which includes those non-active atoms in the extension shell that have a distance of less than the sum"
                  " of their VDW radii plus Dist AtomsAroundOpt. The second option is to create a shell of atoms in the"
                  " non-active region at a specific distance from the active region. The user will be asked to provide the"
                  " distance.")
            # Ask the user to choose the atoms belonging to the extension shell
            Ex_sh_options =  input("Select an option to choose the atoms belonging to optRegionExt: " + "\n"
                                   "Option 1       - Do not use extended active region" + "\n"
                                   "Option 2       - Add only atoms bonded covalently to active atoms" + "\n"
                                   "Option 3       - (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1))"+ "\n"
                                   "Option 4       - Manually (using our code) define the extended active region (Default: empty list)" + "\n"
                                   "Enter your option: ")
            # option 1
            if Ex_sh_options == "1":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                #qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB true \n")
                #

                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    #qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    #qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB true \n")
                qmmm_file.write(f"# Option 1: no extension shell will be included in the QM/MM calculations \n")
                qmmm_file.write("ExtendActiveRegion No \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {out_pdb_for_ORCA_QMMM} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 2
            if Ex_sh_options == "2":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                #qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB true \n")
                #
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    #qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    #qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB true \n")
                qmmm_file.write(f"# Option 2: add only atoms bonded covalently to active atoms \n")
                qmmm_file.write("ExtendActiveRegion cov_bonds \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {out_pdb_for_ORCA_QMMM} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 3
            if Ex_sh_options == "3":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                #qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB true \n")
                #
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    #qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    #qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB true \n")
                qmmm_file.write(f"# Option 3: (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) \n")
                qmmm_file.write("ExtendActiveRegion distance \n")
                distance_criterion = input("Enter the Dist_AtomsAroundOpt value (Default; 1): ")
                if distance_criterion == "":
                    qmmm_file.write("Dist_AtomsAroundOpt 1.\n")
                else:
                    qmmm_file.write(f"Dist_AtomsAroundOpt {distance_criterion}\n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")                  
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {out_pdb_for_ORCA_QMMM} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()
                #print(f"We have created the {qmmm_file} file for your QM/MM calculations")
            # option 4
            if Ex_sh_options == "4":
                if atoms_or_groups.lower() == "a":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in mm_list_sorted:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    #list_resi_num_es = []
                    #list_all_es = []
                    #for line in data_pdb1.split('\n')[:-1]:
                        #if line.startswith("ATOM") or line.startswith("HETATM"):
                            #if int(line[6:11].strip()) in es_list:
                                #atm_num = line[22:26].strip()
                                #list_resi_num_es.append(atm_num)
                                #for line in data_pdb1.split('\n')[:-1]:
                                    #if line[22:26].strip() == atm_num:
                                        #list_all_es.append(line[6:11].strip())
                                    #else:
                                        #continue
                            #else:
                                #continue
                        #else:
                            #continue
    
                    #new_list_all_es = []
                    #for one_mm in list_all_es:
                        #if int(one_mm) not in new_list_all_es:
                            #new_list_all_es.append(int(one_mm))
                        #else:
                            #continue
                    es_list_new = set(es_list) - set(mm_list_sorted)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                elif atoms_or_groups.lower() == "g":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in new_list_all_mm:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    list_resi_num_es = []
                    list_all_es = []
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(line[6:11].strip()) in es_list:
                                atm_num = line[20:26].strip()
                                list_resi_num_es.append(atm_num)
                                for line in data_pdb1.split('\n')[:-1]:
                                    if line[20:26].strip() == atm_num:
                                        list_all_es.append(line[6:11].strip())
                                    else:
                                        continue
                            else:
                                continue
                        else:
                            continue
    
                    new_list_all_es = []
                    for one_mm in list_all_es:
                        if int(one_mm) not in new_list_all_es:
                            new_list_all_es.append(int(one_mm))
                        else:
                            continue
                    #es_list_new = set(new_list_all_es)                
                    es_list_new = set(new_list_all_es) - set(new_list_all_mm)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                # Assign the output PDB file for QM/MM calculations by ORCA
                #out_pdb_for_ORCA_QMMM = input(
                    #"Enter the name of your output PDB file, to assign the QM atoms and active region in its "
                    #"occupancy and B-factor columns, respectively. (out_pdb_ORCA_QMMM.pdb): ")
                #if out_pdb_for_ORCA_QMMM == "":
                    #out_pdb_for_ORCA_QMMM = "out_pdb_ORCA_QMMM.pdb"
                #else:
                    #out_pdb_for_ORCA_QMMM == out_pdb_for_ORCA_QMMM
                # We expect from the following lines to read the input PDB file (defined in the beginig of the code), line by line.
                # Read the mm_list_sorted list item by item.
                # If indexes 6–11 in a line of the input pdb file are equal to the item, replace indexes 64–67 of that line with 1.00.
                #lower = 0
                #upper = len(mm_list_sorted)
                #low_pdb = 0
                #line_in_input_pdb = data_pdb1.split("\n")[:-1]
                #up_pdb = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
                #with open(out_pdb_for_ORCA_QMMM, "w") as output:
                    #for mm_atom in range(lower, upper):
                        #if mm_atom <= lower:
                            #mm_list_sorted.append("")
                            #for qm_atom in range(low_qm,up_qm):
                                #if qm_atom <= low_qm:
                                    #Linear_QM_Atoms_List_PDB_FORMAT_sorted.append("")
                                    #for line_pdb in range (low_pdb, up_pdb):
                                        #if str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip()and str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + "  1.00"+ line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + line_in_input_pdb[line_pdb][61:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:63] + "1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                        #else:
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb] + "\n"
                                        #output.write(line_in_output_pdb)
                    #output.write("END                                                                ")
                # Write the information into the ORCA input file (qmmm.inp).
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmmm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmmm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    qmmm_file.write("! RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                qmmm_file.write("!QMMM \n")
                qmmm_file.write("%qmmm \n")
                topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                if topology_file == "":
                    topology_file = "prmtop.ORCAFF.prms"
                else:
                    topology_file == topology_file
                qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')
                #qmmm_file.write(f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                #                f"included in the active region and from which {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are"
                #                f" included in the QM region (excluding junctions). \n")
                #qmmm_file.write(
                #    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                #qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the active region. \n")
                #qmmm_file.write(
                #    f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                #pdb_or_inp = input(f"Read the active region from the ORCA input (o) or the PDB (p) file? ")
                #if pdb_or_inp == "o" or pdb_or_inp == "":
                qmmm_file.write(
                    f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region, as follows. \n")
                #qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                qmmm_file.write("Use_QM_InfoFromPDB true \n")
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(new_list_all_mm)} atoms are included in the active "
                                f"region, as follows. \n")
                    #qmmm_file.write(f"{mm_lsOrca_new} \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms "
                                f"(from the QM region). A total of {len(mm_list_sorted)} atoms are included in the active "
                                f"region, as follows. \n")                    
                    #qmmm_file.write(f"{mm_lsOrca} \n")
                qmmm_file.write("Use_Active_InfoFromPDB true \n")
                #elif E_Shell_List.lower() == "yes" or E_Shell_List.lower() == "y":
                qmmm_file.write(f"# An extension shell of {extension_shell_thickness}-angstrom thickness will surround the "
                                    f"active region. \n")
                if atoms_or_groups.lower() == "g":
                    qmmm_file.write(
                        f"# There are {atom_number_1} atoms in your PDB file, from which {len(new_list_all_mm)} atoms are "
                        f"included in the active region and {len(es_list_sorted_PDB_format)} atoms are included in the"
                        f" extension shell. Hence, {atom_number_1 - len(es_list_new) - len(new_list_all_mm)}"
                        f" atoms will be ignored in the calculations. \n")
                elif atoms_or_groups.lower() == "a":
                    qmmm_file.write(
                        f"# There are {atom_number_1} atoms in your PDB file, from which {len(mm_list_sorted)} atoms are "
                        f"included in the active region and {len(es_list_sorted_PDB_format)} atoms are included in the"
                        f" extension shell. Hence, {atom_number_1 - len(es_list_new) - len(mm_list_sorted)}"
                        f" atoms will be ignored in the calculations. \n")                    
                qmmm_file.write(f"{es_lsOrca_new} \n")
                qmmm_file.write("end \n")
                # Ask geom scan
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")                
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {out_pdb_for_ORCA_QMMM} \n")
                qmmm_file.write(" \n")
                qmmm_file.close()                
        else:
            break
        print("All necessary files for multi-scale calculations have been generated in your working directory. You can now execute another command or exit the program.")
    ######################################################################
    # The TER command. TER lines between molecular structures of organic solvents.
    elif command.lower() == "ter":
        print("This command inserts TER lines between all residues.")
        print("It is proper for PDB files containing small molecules solvated with organic solvents.")
        print("These PDB files are generated by Packmol.")
        print("To inset TER between protein chains use the 'terp' command.")
        pdb_text_1 = open(pdb_file_1, "r")
        data_pdb1 = pdb_text_1.readlines()
        out_pdb_ter = input("Enter the name of your output PDB file with TER (overwrite): ")
        if out_pdb_ter == "":
            out_pdb_ter = pdb_file_1
        else:
            out_pdb_ter = out_pdb_ter
        #lower = 0
        #upper = len(data_pdb1)-1
        def append_ter_line(out_pdb_ter, text_to_append):
            with open(out_pdb_ter, "w") as file_object:
                for line_pdb in range(len(data_pdb1)):
                    if line_pdb == len(data_pdb1)-1:
                        line_in_output_pdb = data_pdb1[line_pdb]
                    elif str(data_pdb1[line_pdb][17:20]) != str(data_pdb1[line_pdb + 1][17:20]):
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    elif str(data_pdb1[line_pdb][17:26]) != str(data_pdb1[line_pdb + 1][17:26]):
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    else:
                        line_in_output_pdb = data_pdb1[line_pdb]
                    file_object.write(line_in_output_pdb)
        append_ter_line(out_pdb_ter, "TER")
        print(f"Your PDB file ({out_pdb_ter}) with TERs is generated.")
    ######################################################################
    elif command.lower() == "terp":
        out_pdb = "tided_up_pdb.tmp"
        pdb_text_1 = open(pdb_file_1, "r")
        data_pdb1 = pdb_text_1.read()
        with open(out_pdb, "w") as output:
            for line in data_pdb1.split("\n"):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    output.write(line + "\n")
        #print(f"Your tided up PDB file ({out_pdb}) is generated.")
        pdb_text_1 = open(out_pdb, "r")
        data_pdb1 = pdb_text_1.readlines()
        out_pdb_ter = input("Enter the name of your output PDB file with TER (overwrite): ")
        if out_pdb_ter == "":
            out_pdb_ter = pdb_file_1
        else:
            out_pdb_ter = out_pdb_ter
            lower = 0
            upper = len(data_pdb1) - 1
        def append_terp_line(out_pdb_ter, text_to_append):
            with open(out_pdb_ter, "w") as file_object:
                for line_pdb in range(lower, upper):
                    # chain
                    if str(data_pdb1[line_pdb][21:22]) != str(data_pdb1[line_pdb + 1][21:22]):
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    # missing residue inside a chain
                    elif int(data_pdb1[line_pdb][22:26].strip())+1 != int(data_pdb1[line_pdb+1][22:26].strip()) and data_pdb1[line_pdb][17:20] != data_pdb1[line_pdb+1][17:20]: 
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    # OXT atom
                    elif str(data_pdb1[line_pdb][13:16]) == "OXT":
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"                          
                    # between ATOM and HETATM
                    elif str(data_pdb1[line_pdb][0:6]) != str(data_pdb1[line_pdb+1][0:6]):
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    # HETATM
                    elif data_pdb1[line_pdb].startswith("HETATM") and data_pdb1[line_pdb][17:22] != data_pdb1[line_pdb+1][17:22]:
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    elif data_pdb1[line_pdb].startswith("HETATM") \
                            and data_pdb1[line_pdb][17:22] == data_pdb1[line_pdb+1][17:22] \
                            and data_pdb1[line_pdb][23:26] != data_pdb1[line_pdb+1][23:26]:
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    elif data_pdb1[line_pdb].startswith("ATOM") \
                            and data_pdb1[line_pdb][17:20] == "WAT"\
                            and data_pdb1[line_pdb][17:22] == data_pdb1[line_pdb+1][17:22] \
                            and data_pdb1[line_pdb][23:26] != data_pdb1[line_pdb+1][23:26]:
                        line_in_output_pdb = str(data_pdb1[line_pdb]) + text_to_append + "\n"
                    else:
                        line_in_output_pdb = data_pdb1[line_pdb]
                    file_object.write(line_in_output_pdb)
        append_terp_line(out_pdb_ter, "TER")
        pdb_text_1.close()
        #out_pdb.close()
        os.remove('tided_up_pdb.tmp')
    elif command.lower() == "charge" or command.lower() == "cc":
        inp_pdb_text = open(pdb_file_1, "r")
        lines_inp_pdb = inp_pdb_text.readlines()
        # Read charges of each atoms in PDB file
        total_charge = 0
        for line in lines_inp_pdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                total_charge += float(line[54:64].strip())
        print("Total charge of the system is :", "{:.6f}".format(total_charge))
        ######################################################################
    elif command.lower() == "qmqm2" or command.lower() == "qm2":
        #pdb_or_inp = input(f"Read the QM region and active region from an ORCA input (o) or a PDB (p) file? ")
        #if pdb_or_inp.lower() == "o" or pdb_or_inp == "":
            #with open("orcapdb.pdb", "w") as output:
                #for line in data_pdb1.split("\n"):
                    #if line.startswith("ATOM") or line.startswith("HETATM"):
                        #rline = line[0:55] + "  0.00  0.00          " + line[76:] + "\n"
                        #output.write(rline)
                #output.write("END   00000          00000       0.000  0.0000  0.0000   0.00  0.00")

            # Read contents of the temp.pdb file to a string
            #pdb_file_orca = 'orcapdb.pdb'
            #pdb_text_1 = open("orcapdb.pdb", "r")
            #data_pdb1 = pdb_text_1.read()
            # Ask the user for the S1 file name, the file describing the QM region
            s1_file = input("Please provide the file name that contains the QM region description (s1): ")
            if s1_file == "":
                s1_file = "s1"
            else:
                s1_file == s1_file
            # Open the S1 file, read its lines and connect the lines
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                linear_string = str(",".join(s1_line.strip() for s1_line in s1_lines_No_Blank_lines))
                # Convert the numeric string ranges to a list (linear_qm_atoms_list). E.g., convert 5-7 to 5, 6, 7
                def f(linear_string):
                    result = []
                    for part in linear_string.split(','):
                        if '-' in part:
                            a, b = part.split('-')
                            a, b = int(a), int(b)
                            result.extend(range(a, b + 1))
                        else:
                            a = int(part)
                            result.append(a)
                    return result
                Linear_QM_Atoms_List_PDB_FORMAT = (f(linear_string))
                Linear_QM_Atoms_List_PDB_FORMAT_sorted = sorted(Linear_QM_Atoms_List_PDB_FORMAT)
            #print(f"QM atoms list in PDB format: {Linear_QM_Atoms_List_PDB_FORMAT_sorted}")
            ################################################
            # Convert the S1 file to the ORCA input file format e.g., {1 5:10 12}
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                qm_atoms_for_orca_input_file = ""
                for qm_atoms in s1_lines_No_Blank_lines:
                    qm_atoms_orca_format = qm_atoms.strip().split('-')
                    qm_atoms_orca_format = [int(element) - 1 for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = [str(element) for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = ":".join(qm_atoms_orca_format)
                    qm_atoms_for_orca_input_file += qm_atoms_orca_format + " "
            Qm_lst = ['QMatoms ', '{', qm_atoms_for_orca_input_file, '}', ' end']
            QM_atoms_list_ORCA_format = ''.join(map(str, Qm_lst))
            print(f"The following {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are included in the QM region (excluding junctions).")
            print(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
            #########################################################################
            print(
                f"The PDB file containing the QM atoms (QM.pdb) is generated.")            
            pdb_text_11 = open(pdb_file_1, "r")
            data_pdb11 = pdb_text_11.readlines()
            with open ("QM.pdb", "w") as output:
                for Pline in data_pdb11:
                    if int(Pline[6:11].strip()) in Linear_QM_Atoms_List_PDB_FORMAT_sorted:
                        output.write(Pline)           
            #########################################################################            
            print("ORCA identifies junctions itself and will include them in the QM region.")            
            ##################################################################################    
            # Ask the user for the thickness of the active region
            active_region_thickness = input("Enter the thickness of the QM2 region around the original QM region (default is 3.0 Angstrom): ")
            if active_region_thickness == "":
                active_region_thickness = float(3.0)
            else:
                active_region_thickness = float(active_region_thickness)
            # Define the variables
            x_qm = 0.00
            y_qm = 0.00
            z_qm = 0.00
            x_mm = 0.00
            y_mm = 0.00
            z_mm = 0.00
            SD = 0.00
            # Define the new list
            atoms_or_groups = input("Do you want to include all atoms within this distance from the QM region "
                                    "in the QM2 region (a), or the groups that touch or cross this distance "
                                    "from the QM region (g)? ")
            if atoms_or_groups.lower() == "g":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ##
                list_resi_num = []
                list_all_mm = []
                for line in data_pdb1.split('\n')[:-1]:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if int(line[6:11].strip()) in mm_list_sorted:
                            atm_num = line[20:26].strip()
                            list_resi_num.append(atm_num)
                            for line in data_pdb1.split('\n')[:-1]:
                                if line[20:26].strip() == atm_num:
                                    list_all_mm.append(line[6:11].strip())
                                else:
                                    continue
                        else:
                            continue
                    else:
                        continue
                new_list_all_mm = []
                for one_mm in list_all_mm:
                    if int(one_mm) not in new_list_all_mm:
                        new_list_all_mm.append(int(one_mm))
                    else:
                        continue
                # Remove duplicates from the QM and MM atoms lists
                qm2mm_list = set(new_list_all_mm) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                qm2mm_sortedlist = sorted(list(set(qm2mm_list)))
                ## simplify the list
                # mm_list_sortedf = mm_list_sorted[:-1]
                # mmlistoc = []
                # for i in range(0,len(mm_list_sortedf)):
                # if i <= 0:
                # for j in range (0,len(mm_list_sortedf)):
                # if mm_list_sortedf [i] in range(mm_list_sortedf[i] and mm_list_sortedf[j]+1) and mm_list_sortedf[j]+1 not in mm_list_sortedf and mm_list_sortedf[i]-1 not in mm_list_sortedf:
                # if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                # if mm_list_sortedf[i] != mm_list_sortedf[j]:
                # mmlistoc.append(str(mm_list_sortedf[i]-1)+":"+str(mm_list_sortedf[j]-1))
                # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                # i = j+1
                # else:
                # mmlistoc.append(str(mm_list_sortedf[i]-1))
                # print(str(mm_list_sortedf[i]))
                # i = j+1
                # lsf=str(mmlistoc).replace('[','')
                # lsff=str(lsf).replace(']','')
                # lsfff=str(lsff).replace("'","")
                # lsffff=str(lsfff).replace(",","")
                # mm_lsOrca = 'ActiveAtoms'+' {'+ lsffff+'} '+ 'end'
                # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                ## simplify the new-list_mm
                oc_list_all_mm = new_list_all_mm
                mmlistoc_new = []
                for i in range(0, len(oc_list_all_mm)):
                    if i <= 0:
                        for j in range(0, len(oc_list_all_mm)):
                            if oc_list_all_mm[i] in range(oc_list_all_mm[i] and oc_list_all_mm[j] + 1) and oc_list_all_mm[
                                j] + 1 not in oc_list_all_mm and oc_list_all_mm[i] - 1 not in oc_list_all_mm:
                                if oc_list_all_mm[i] != min(oc_list_all_mm) or oc_list_all_mm[j] != max(oc_list_all_mm):
                                    if oc_list_all_mm[i] != oc_list_all_mm[j]:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1) + ":" + str(oc_list_all_mm[j] - 1))
                                        # print(str(oc_list_all_mm[i])+":"+str(oc_list_all_mm[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1))
                                        # print(str(oc_list_all_mm[i]))
                                        i = j + 1
                        lsf_new = str(mmlistoc_new).replace('[', '')
                        lsff_new = str(lsf_new).replace(']', '')
                        lsfff_new = str(lsff_new).replace("'", "")
                        lsffff_new = str(lsfff_new).replace(",", "")
                        mm_lsOrca_new = 'ActiveAtoms' + ' {' + lsffff_new + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff_new +'} '+ 'end')
                # print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                # print(mm_list_sorted)
                print(f"The following {len(new_list_all_mm)} atoms are included in the QM2 region.")
                print(new_list_all_mm)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")

            elif atoms_or_groups.lower() == "a":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                qm2mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                qm2mm_sortedlist = sorted(list(set(qm2mm_list)))
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['ActiveAtoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ## simplify the list
                mm_list_sortedf = mm_list_sorted
                mmlistoc = []
                for i in range(0, len(mm_list_sortedf)):
                    if i <= 0:
                        for j in range(0, len(mm_list_sortedf)):
                            if mm_list_sortedf[i] in range(mm_list_sortedf[i] and mm_list_sortedf[j] + 1) and \
                                    mm_list_sortedf[j] + 1 not in mm_list_sortedf and mm_list_sortedf[
                                i] - 1 not in mm_list_sortedf:
                                if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                                    if mm_list_sortedf[i] != mm_list_sortedf[j]:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1) + ":" + str(mm_list_sortedf[j] - 1))
                                        # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1))
                                        # print(str(mm_list_sortedf[i]))
                                        i = j + 1
                        lsf = str(mmlistoc).replace('[', '')
                        lsff = str(lsf).replace(']', '')
                        lsfff = str(lsff).replace("'", "")
                        lsffff = str(lsfff).replace(",", "")
                        mm_lsOrca = 'ActiveAtoms' + ' {' + lsffff + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                print(f"The following {len(mm_list_sorted)} atoms are included in the QM2 region.")
                print(mm_list_sorted)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            else:
                print("type 'a' or 'g'")
                break

            #########################################################################
            print(
                f"The PDB file containing the QM and QM2 atoms (inp_qmqm2.pdb) is generated.")            
            pdb_text_12 = open(pdb_file_1, "r")
            data_pdb12 = pdb_text_12.readlines()
            if atoms_or_groups.lower() == "a":
                with open ("inp_qmqm2.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in mm_list_sorted:
                            output.write(Pline)
            if atoms_or_groups.lower() == "g":
                with open ("inp_qmqm2.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in new_list_all_mm:
                            output.write(Pline)           
            #########################################################################
            qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmqm2.inp): ")
            if qmmm_file == "":
                qmmm_file = "qmqm2.inp"
            else:
                qmmm_file == qmmm_file
            qmmm_file = open(qmmm_file, "w+")
            qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
            qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
            qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
            qmmm_file.write("# DOI: ................ \n")
            # The working directory and the SetUp time
            import os
            WorDir = os.getcwd()
            qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
            from datetime import date,time, datetime
            today = date.today()
            now = datetime. now()
            current_time = time(now.hour, now.minute, now.second)
            qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
            print("We are setting up the level of of theory.")
            print("The default level of theory is:")
            print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
            theory = input("Do you accept the default level of theory (yes/no): ").lower()
            if theory == "" or theory == "yes":
                functional = "TPSSh"
                basis = "def2-SVP"
                qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
            elif theory == "no":
                print("You have not accepted the default level of theory")
                functional = input("Enter the name of your desired functional: ")
                basis = input("Enter the name of your desired basis set: ")
                qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
            print("The default QM2 method is XTB")
            method = input("Do you accept the QM2 method (yes/no): ").lower()
            if method == "" or method == "yes":
                qmmm_file.write("!QM/XTB \n")
                qmmm_file.write("%qmmm \n")
                qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                if atoms_or_groups.lower() == "a":
                    qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the QM1+QM2 region. \n")
                elif atoms_or_groups.lower() == "g":
                    qmmm_file.write(f"# A total of {len(new_list_all_mm)} atoms are included in the QM1+QM2 region. \n")
                qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")          
                #qmmm_file.write(f"{mm_lsOrca} \n")
  
                # Ask the user for the charge of the medium (QM1+QM2) system
                charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                if charge_of_the_medium_region == "":
                    charge_of_the_medium_region = int(0)
                else:
                    charge_of_the_medium_region = int(charge_of_the_medium_region)
                # Ask the user for the multiplicity of the medium (QM1+QM2) system
                multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                if multiplicity_of_the_medium_region == "":
                    multiplicity_of_the_medium_region = int(1)
                else:
                    multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")                        
                qmmm_file.write("end \n")
                geom = input("Do you want to do geom scan (yes/no): ").lower()
                if geom == "y" or geom == "yes":
                    qmmm_file.write(f"% geom scan \n")  
                    atom1 = input("Enter the number of the first atom: ").lower()
                    orca_atom1 = str(int(atom1)-1)
                    atom2 = input("Enter the number of the sencond atom: ").lower()
                    orca_atom2 = str(int(atom2)-1)
                    dist1= input("Enter the initial distance: ").lower()
                    dist2 = input("Enter the final distance: ").lower()
                    num = input("Enter the number of scan proceses: ").lower()
                    qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                    qmmm_file.write("end \n")
                    qmmm_file.write("end \n")
                # Ask the user for the number of CPU cores for parallel calculations
                cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                if cpu_cores == "":
                    cpu_cores = int(28)
                else:
                    cpu_cores = int(cpu_cores)
                qmmm_file.write("% pal \n")
                qmmm_file.write("nprocs \n")
                qmmm_file.write(f"{cpu_cores} \n")
                qmmm_file.write("end \n")
                # Ask the user for the charge of the active region
                charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                if charge_of_the_QM_region == "":
                    charge_of_the_QM_region = int(0)
                else:
                    charge_of_the_QM_region = int(charge_of_the_QM_region)
                # Ask the user for the multiplicity of the active region
                multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                if multiplicity_of_the_QM_region == "":
                    multiplicity_of_the_QM_region = int(1)
                else:
                    multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                #pdb_file_orca = 'orcapdb.pdb'
                qmmm_file.write(
                    f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {'inp_qmqm2.pdb'}  # charge and mult. of the high level region\n")
                qmmm_file.write(" \n")
                qmmm_file.close()
            #######################    
            elif method == "no":
                print("You have not accepted the default QM2 method")
                Q2_method = input("Enter the name of your desired QM2 method: ")
                if Q2_method == "QM2" or Q2_method == "qm2":
                    print("You should determine your functional and basis set")
                    functional_q2 = input("Enter the name of your desired functional: ")
                    basis_q2 = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"!QM/QM2 \n")
                    qmmm_file.write("%qmmm \n")
                    qmmm_file.write(
                        f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the QM1+QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm)} atoms are included in the QM1+QM2 region. \n")
                    qmmm_file.write(f"# The thickness of the active region is {active_region_thickness} angstroms (from the QM region). \n")
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")          
                    qmmm_file.write(f'QM2CUSTOMMETHOD "{functional_q2}" \n')
                    qmmm_file.write(f'QM2CUSTOMBASIS  "{basis_q2} def2/J"\n')
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")                        
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {'inp_qmqm2.pdb'}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()

                    #####################
                else:
                    qmmm_file.write(f"!QM/{Q2_method} \n")
                    qmmm_file.write("%qmmm \n")
                    qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")                 
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted)} atoms are included in the QM1+QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm)} atoms are included in the QM1+QM2 region. \n")
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")          
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")                        
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {'inp_qmqm2.pdb'}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()
            print("All necessary files for multi-scale calculations have been generated in your working directory. You can now execute another command or exit the program.")
            ################################################################################################## 

    elif command.lower() == "qmqm2mm" or command.lower() == "qm2mm":
        #pdb_or_inp = input(f"Read the QM region and active region from an ORCA input (o) or a PDB (p) file? ")
        #if pdb_or_inp.lower() == "o" or pdb_or_inp == "":
            #with open("orcapdb.pdb", "w") as output:
                #for line in data_pdb1.split("\n"):
                    #if line.startswith("ATOM") or line.startswith("HETATM"):
                        #rline = line[0:55] + "  0.00  0.00          " + line[76:] + "\n"
                        #output.write(rline)
                #output.write("END   00000          00000       0.000  0.0000  0.0000   0.00  0.00")

            # Read contents of the temp.pdb file to a string
            #pdb_file_orca = 'orcapdb.pdb'
            #pdb_text_1 = open("orcapdb.pdb", "r")
            #data_pdb1 = pdb_text_1.read()
            # Ask the user for the S1 file name, the file describing the QM region
            s1_file = input("Please provide the file name that contains the QM region description (s1): ")
            if s1_file == "":
                s1_file = "s1"
            else:
                s1_file == s1_file
            # Open the S1 file, read its lines and connect the lines
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                linear_string = str(",".join(s1_line.strip() for s1_line in s1_lines_No_Blank_lines))
                # Convert the numeric string ranges to a list (linear_qm_atoms_list). E.g., convert 5-7 to 5, 6, 7
                def f(linear_string):
                    result = []
                    for part in linear_string.split(','):
                        if '-' in part:
                            a, b = part.split('-')
                            a, b = int(a), int(b)
                            result.extend(range(a, b + 1))
                        else:
                            a = int(part)
                            result.append(a)
                    return result
                Linear_QM_Atoms_List_PDB_FORMAT = (f(linear_string))
                Linear_QM_Atoms_List_PDB_FORMAT_sorted = sorted(Linear_QM_Atoms_List_PDB_FORMAT)
            #print(f"QM atoms list in PDB format: {Linear_QM_Atoms_List_PDB_FORMAT_sorted}")
            ################################################
            # Convert the S1 file to the ORCA input file format e.g., {1 5:10 12}
            with open(s1_file) as s1_text:
                s1_lines_No_Blank_lines = filter(None, (s1_line.rstrip() for s1_line in s1_text))
                qm_atoms_for_orca_input_file = ""
                for qm_atoms in s1_lines_No_Blank_lines:
                    qm_atoms_orca_format = qm_atoms.strip().split('-')
                    qm_atoms_orca_format = [int(element) - 1 for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = [str(element) for element in qm_atoms_orca_format]
                    qm_atoms_orca_format = ":".join(qm_atoms_orca_format)
                    qm_atoms_for_orca_input_file += qm_atoms_orca_format + " "
            Qm_lst = ['QMatoms ', '{', qm_atoms_for_orca_input_file, '}', ' end']
            QM_atoms_list_ORCA_format = ''.join(map(str, Qm_lst))
            print(f"The following {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are included in the QM region (excluding junctions).")
            print(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
            #########################################################################
            print(
                f"The PDB file containing the QM atoms (QM.pdb) is generated.")            
            pdb_text_11 = open(pdb_file_1, "r")
            data_pdb11 = pdb_text_11.readlines()
            with open ("QM.pdb", "w") as output:
                for Pline in data_pdb11:
                    if int(Pline[6:11].strip()) in Linear_QM_Atoms_List_PDB_FORMAT_sorted:
                        output.write(Pline)           
            #########################################################################            
            print("ORCA identifies junctions itself and will include them in the QM region.")
            ##################################################################################
    
            # Ask the user for the thickness of the active region
            active_region_thickness = input("Enter the thickness of the QM2 region around the original QM region (default is 3.0 Angstrom): ")
            if active_region_thickness == "":
                active_region_thickness = float(3.0)
            else:
                active_region_thickness = float(active_region_thickness)
            # Define the variables
            x_qm = 0.00
            y_qm = 0.00
            z_qm = 0.00
            x_mm = 0.00
            y_mm = 0.00
            z_mm = 0.00
            SD = 0.00
            # Define the new list
            atoms_or_groups = input("Do you want to include all atoms within this distance from the QM region "
                                    "in the QM2 region (a), or the groups that touch or cross this distance "
                                    "from the QM region (g)? ")
            if atoms_or_groups.lower() == "g":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['QM2Atoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ##
                list_resi_num = []
                list_all_mm = []
                for line in data_pdb1.split('\n')[:-1]:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if int(line[6:11].strip()) in mm_list_sorted:
                            atm_num = line[20:26].strip()
                            list_resi_num.append(atm_num)
                            for line in data_pdb1.split('\n')[:-1]:
                                if line[20:26].strip() == atm_num:
                                    list_all_mm.append(line[6:11].strip())
                                else:
                                    continue
                        else:
                            continue
                    else:
                        continue
                new_list_all_mm = []
                for one_mm in list_all_mm:
                    if int(one_mm) not in new_list_all_mm:
                        new_list_all_mm.append(int(one_mm))
                    else:
                        continue
                # Remove duplicates from the QM and MM atoms lists
                qm2mm_list = set(new_list_all_mm) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                qm2mm_sortedlist = sorted(list(set(qm2mm_list)))
                ## simplify the list
                # mm_list_sortedf = mm_list_sorted[:-1]
                # mmlistoc = []
                # for i in range(0,len(mm_list_sortedf)):
                # if i <= 0:
                # for j in range (0,len(mm_list_sortedf)):
                # if mm_list_sortedf [i] in range(mm_list_sortedf[i] and mm_list_sortedf[j]+1) and mm_list_sortedf[j]+1 not in mm_list_sortedf and mm_list_sortedf[i]-1 not in mm_list_sortedf:
                # if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                # if mm_list_sortedf[i] != mm_list_sortedf[j]:
                # mmlistoc.append(str(mm_list_sortedf[i]-1)+":"+str(mm_list_sortedf[j]-1))
                # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                # i = j+1
                # else:
                # mmlistoc.append(str(mm_list_sortedf[i]-1))
                # print(str(mm_list_sortedf[i]))
                # i = j+1
                # lsf=str(mmlistoc).replace('[','')
                # lsff=str(lsf).replace(']','')
                # lsfff=str(lsff).replace("'","")
                # lsffff=str(lsfff).replace(",","")
                # mm_lsOrca = 'ActiveAtoms'+' {'+ lsffff+'} '+ 'end'
                # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                ## simplify the new-list_mm
                #oc_list_all_mm = new_list_all_mm
                #oc_list_all_mm = mm_list_sorted
                oc_list_all_mm = qm2mm_sortedlist
                mmlistoc_new = []
                for i in range(0, len(oc_list_all_mm)):
                    if i <= 0:
                        for j in range(0, len(oc_list_all_mm)):
                            if oc_list_all_mm[i] in range(oc_list_all_mm[i] and oc_list_all_mm[j] + 1) and oc_list_all_mm[
                                j] + 1 not in oc_list_all_mm and oc_list_all_mm[i] - 1 not in oc_list_all_mm:
                                if oc_list_all_mm[i] != min(oc_list_all_mm) or oc_list_all_mm[j] != max(oc_list_all_mm):
                                    if oc_list_all_mm[i] != oc_list_all_mm[j]:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1) + ":" + str(oc_list_all_mm[j] - 1))
                                        # print(str(oc_list_all_mm[i])+":"+str(oc_list_all_mm[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc_new.append(str(oc_list_all_mm[i] - 1))
                                        # print(str(oc_list_all_mm[i]))
                                        i = j + 1
                        lsf_new = str(mmlistoc_new).replace('[', '')
                        lsff_new = str(lsf_new).replace(']', '')
                        lsfff_new = str(lsff_new).replace("'", "")
                        lsffff_new = str(lsfff_new).replace(",", "")
                        mm_lsOrca_new = 'QM2Atoms' + ' {' + lsffff_new + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff_new +'} '+ 'end')
                # print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                # print(mm_list_sorted)
                print(f"The following {len(new_list_all_mm)} atoms are included in the QM2 region.")
                print(new_list_all_mm)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")

            elif atoms_or_groups.lower() == "a":
                mm_list = []
                for item in Linear_QM_Atoms_List_PDB_FORMAT:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= active_region_thickness:
                                            mm = int(lines[6:11])
                                            mm_list.append(mm)
                # Remove duplicates from the MM list
                mm_list = list(dict.fromkeys(mm_list))
                # Remove duplicates from the QM and MM atoms lists
                qm2mm_list = set(mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                qm2mm_sortedlist = sorted(list(set(qm2mm_list)))
                # Sort the MM atoms list
                mm_list_sorted = sorted(list(set(mm_list)))
                MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in mm_list_sorted)
                mm_lst = ['QM2Atoms ', '{', MM_list_ORCA_FORMAT, '}', ' end']
                mm_lst2 = ''.join(map(str, mm_lst))
                ## simplify the list
                #mm_list_sortedf = mm_list_sorted
                mm_list_sortedf = qm2mm_sortedlist
                mmlistoc = []
                for i in range(0, len(mm_list_sortedf)):
                    if i <= 0:
                        for j in range(0, len(mm_list_sortedf)):
                            if mm_list_sortedf[i] in range(mm_list_sortedf[i] and mm_list_sortedf[j] + 1) and \
                                    mm_list_sortedf[j] + 1 not in mm_list_sortedf and mm_list_sortedf[
                                i] - 1 not in mm_list_sortedf:
                                if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                                    if mm_list_sortedf[i] != mm_list_sortedf[j]:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1) + ":" + str(mm_list_sortedf[j] - 1))
                                        # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                                        i = j + 1
                                    else:
                                        mmlistoc.append(str(mm_list_sortedf[i] - 1))
                                        # print(str(mm_list_sortedf[i]))
                                        i = j + 1
                        lsf = str(mmlistoc).replace('[', '')
                        lsff = str(lsf).replace(']', '')
                        lsfff = str(lsff).replace("'", "")
                        lsffff = str(lsfff).replace(",", "")
                        mm_lsOrca = 'QM2Atoms' + ' {' + lsffff + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                print(f"The following {len(mm_list_sorted)} atoms are included in the QM2 region.")
                print(mm_list_sorted)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} atoms are also in the QM region.")
            else:
                print("type 'a' or 'g'")
                break

            #########################################################################
            print(
                f"The PDB file containing the QM1 and QM2 atoms (QM1-QM2.pdb) is generated.")            
            pdb_text_12 = open(pdb_file_1, "r")
            data_pdb12 = pdb_text_12.readlines()
            if atoms_or_groups.lower() == "a":
                with open ("QM1-QM2.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in mm_list_sorted:
                            output.write(Pline)
            if atoms_or_groups.lower() == "g":
                with open ("QM1-QM2.pdb", "w") as output:
                    for Pline in data_pdb12:
                        if int(Pline[6:11].strip()) in new_list_all_mm:
                            output.write(Pline)           
            #########################################################################			
	    # ezafe shode 
	    ##########################################################################
            # Ask the user for the thickness of the active region
            mm_region_thickness = input("Enter the thickness of the active region (default is 6.0 Angstrom): ")
            if mm_region_thickness == "":
                mm_region_thickness = float(6.0)
            else:
                mm_region_thickness = float(mm_region_thickness)
            # Define the variables
            x_qm = 0.00
            y_qm = 0.00
            z_qm = 0.00
            x_mm = 0.00
            y_mm = 0.00
            z_mm = 0.00
            SD = 0.00
            # Define the new list
            mm_groups = input("Do you want to include all atoms within this distance from the QM2 system "
                                    "in the active region (a), or the groups that touch or cross this distance "
                                    "from the QM2 system (g)? ")
            if mm_groups.lower() == "g" or mm_groups.lower() == "":
                q2mm_list = []
                for item in qm2mm_sortedlist:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= mm_region_thickness:
                                            mm = int(lines[6:11])
                                            q2mm_list.append(mm)
                # Remove duplicates from the MM list
                q2mm_list = list(dict.fromkeys(q2mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # q2mm_list = set(q2mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                q2mm_list_sorted = sorted(list(set(q2mm_list)))
                q2MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in q2mm_list_sorted)
                q2mm_lst = ['ActiveAtoms ', '{', q2MM_list_ORCA_FORMAT, '}', ' end']
                q2mm_lst2 = ''.join(map(str, q2mm_lst))
                ##
                list_resi_numq2 = []
                list_all_mmq2 = []
                for line in data_pdb1.split('\n')[:-1]:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if int(line[6:11].strip()) in q2mm_list_sorted:
                            atm_num = line[20:26].strip()
                            list_resi_numq2.append(atm_num)
                            for line in data_pdb1.split('\n')[:-1]:
                                if line[20:26].strip() == atm_num:
                                    list_all_mmq2.append(line[6:11].strip())
                                else:
                                    continue
                        else:
                            continue
                    else:
                        continue
                new_list_all_mmq2 = []
                for one_mm in list_all_mmq2:
                    if int(one_mm) not in new_list_all_mmq2:
                        new_list_all_mmq2.append(int(one_mm))
                    else:
                        continue
                ## simplify the list
                # mm_list_sortedf = mm_list_sorted[:-1]
                # mmlistoc = []
                # for i in range(0,len(mm_list_sortedf)):
                # if i <= 0:
                # for j in range (0,len(mm_list_sortedf)):
                # if mm_list_sortedf [i] in range(mm_list_sortedf[i] and mm_list_sortedf[j]+1) and mm_list_sortedf[j]+1 not in mm_list_sortedf and mm_list_sortedf[i]-1 not in mm_list_sortedf:
                # if mm_list_sortedf[i] != min(mm_list_sortedf) or mm_list_sortedf[j] != max(mm_list_sortedf):
                # if mm_list_sortedf[i] != mm_list_sortedf[j]:
                # mmlistoc.append(str(mm_list_sortedf[i]-1)+":"+str(mm_list_sortedf[j]-1))
                # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                # i = j+1
                # else:
                # mmlistoc.append(str(mm_list_sortedf[i]-1))
                # print(str(mm_list_sortedf[i]))
                # i = j+1
                # lsf=str(mmlistoc).replace('[','')
                # lsff=str(lsf).replace(']','')
                # lsfff=str(lsff).replace("'","")
                # lsffff=str(lsfff).replace(",","")
                # mm_lsOrca = 'ActiveAtoms'+' {'+ lsffff+'} '+ 'end'
                # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                ## simplify the new-list_mm
                oc_list_all_mmq2 = new_list_all_mmq2
                q2mmlistoc_new = []
                for i in range(0, len(oc_list_all_mmq2)):
                    if i <= 0:
                        for j in range(0, len(oc_list_all_mmq2)):
                            if oc_list_all_mmq2[i] in range(oc_list_all_mmq2[i] and oc_list_all_mmq2[j] + 1) and oc_list_all_mmq2[
                                j] + 1 not in oc_list_all_mmq2 and oc_list_all_mmq2[i] - 1 not in oc_list_all_mmq2:
                                if oc_list_all_mmq2[i] != min(oc_list_all_mmq2) or oc_list_all_mmq2[j] != max(oc_list_all_mmq2):
                                    if oc_list_all_mmq2[i] != oc_list_all_mmq2[j]:
                                        q2mmlistoc_new.append(str(oc_list_all_mmq2[i] - 1) + ":" + str(oc_list_all_mmq2[j] - 1))
                                        # print(str(oc_list_all_mmq2[i])+":"+str(oc_list_all_mmq2[j]))
                                        i = j + 1
                                    else:
                                        q2mmlistoc_new.append(str(oc_list_all_mmq2[i] - 1))
                                        # print(str(oc_list_all_mmq2[i]))
                                        i = j + 1
                        lsf_new = str(q2mmlistoc_new).replace('[', '')
                        lsff_new = str(lsf_new).replace(']', '')
                        lsfff_new = str(lsff_new).replace("'", "")
                        lsffff_new = str(lsfff_new).replace(",", "")
                        q2mm_lsOrca_new = 'ActiveAtoms' + ' {' + lsffff_new + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff_new +'} '+ 'end')
                # print(f"The following {len(mm_list_sorted)} atoms are included in the active region.")
                # print(mm_list_sorted)
                print(f"The following {len(new_list_all_mmq2)} atoms are included in the active region.")
                print(new_list_all_mmq2)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} and {len(qm2mm_sortedlist)} atoms are also in the QM and QM2 system, respectively.")
            elif mm_groups.lower() == "a":
                q2mm_list = []
                for item in qm2mm_sortedlist:
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(item) == int(line[6:11]):
                                x_qm = float(line[30:38])
                                y_qm = float(line[38:46])
                                z_qm = float(line[46:54])
                                for lines in data_pdb1.split('\n')[:-1]:
                                    if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                        x_mm = float(lines[30:38])
                                        y_mm = float(lines[38:46])
                                        z_mm = float(lines[46:54])
                                        SD = ((x_mm - x_qm) ** 2 + (y_mm - y_qm) ** 2 + (z_mm - z_qm) ** 2) ** 0.5
                                        if SD <= mm_region_thickness:
                                            mm = int(lines[6:11])
                                            q2mm_list.append(mm)
                # Remove duplicates from the MM list
                q2mm_list = list(dict.fromkeys(q2mm_list))
                # Remove duplicates from the QM and MM atoms lists
                # q2mm_list = set(q2mm_list) - set(Linear_QM_Atoms_List_PDB_FORMAT)
                # Sort the MM atoms list
                q2mm_list_sorted = sorted(list(set(q2mm_list)))
                q2MM_list_ORCA_FORMAT = " ".join(str(active_atom - 1) for active_atom in q2mm_list_sorted)
                q2mm_lst = ['ActiveAtoms ', '{', q2MM_list_ORCA_FORMAT, '}', ' end']
                q2mm_lst2 = ''.join(map(str, q2mm_lst))
                ## simplify the list
                q2mm_list_sortedf = q2mm_list_sorted
                q2mmlistoc = []
                for i in range(0, len(q2mm_list_sortedf)):
                    if i <= 0:
                        for j in range(0, len(q2mm_list_sortedf)):
                            if q2mm_list_sortedf[i] in range(q2mm_list_sortedf[i] and q2mm_list_sortedf[j] + 1) and \
                                    q2mm_list_sortedf[j] + 1 not in q2mm_list_sortedf and q2mm_list_sortedf[
                                i] - 1 not in q2mm_list_sortedf:
                                if q2mm_list_sortedf[i] != min(q2mm_list_sortedf) or q2mm_list_sortedf[j] != max(q2mm_list_sortedf):
                                    if q2mm_list_sortedf[i] != q2mm_list_sortedf[j]:
                                        q2mmlistoc.append(str(q2mm_list_sortedf[i] - 1) + ":" + str(q2mm_list_sortedf[j] - 1))
                                        # print(str(mm_list_sortedf[i])+":"+str(mm_list_sortedf[j]))
                                        i = j + 1
                                    else:
                                        q2mmlistoc.append(str(q2mm_list_sortedf[i] - 1))
                                        # print(str(mm_list_sortedf[i]))
                                        i = j + 1
                        lsf = str(q2mmlistoc).replace('[', '')
                        lsff = str(lsf).replace(']', '')
                        lsfff = str(lsff).replace("'", "")
                        lsffff = str(lsfff).replace(",", "")
                        q2mm_lsOrca = 'ActiveAtoms' + ' {' + lsffff + '} ' + 'end'
                        # print('ActiveAtoms'+' {'+ lsff +'} '+ 'end')
                print(f"The following {len(q2mm_list_sorted)} atoms are included in the active region.")
                print(q2mm_list_sorted)
                print(
                    f"From the above list, {len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)} and {len(qm2mm_sortedlist)} atoms are also in the QM and QM2 system, respectively.")
            else:
                print("type 'a' or 'g'")
                break
            #########################################################################
            print(
                f"The PDB file containing the QM1, QM2, and MM active regions (QM1-QM2-MMactive.pdb) is generated.")            
            pdb_text_13 = open(pdb_file_1, "r")
            data_pdb13 = pdb_text_13.readlines()
            if mm_groups.lower() == "a":
                with open ("QM1-QM2-MMactive.pdb", "w") as output:
                    for Pline in data_pdb13:
                        if int(Pline[6:11].strip()) in q2mm_list_sorted:
                            output.write(Pline)
            if mm_groups.lower() == "g":
                with open ("QM1-QM2-MMactive.pdb", "w") as output:
                    for Pline in data_pdb13:
                        if int(Pline[6:11].strip()) in new_list_all_mmq2:
                            output.write(Pline)           
            #########################################################################            
            print("Now, we are assigning the extension shell. Be patient, this might take a while.")
            print("ORCA optimizes only the positions of active atoms in its QM/MM geometry optimization algorithm."
                  "Nevertheless, the forces exerted on these active atoms are influenced by their interactions with the"
                  " non-active atoms around them. To obtain a smooth convergence of quasi-Newton algorithms in internal"
                  " coordinates, the Hessian values between the active atoms and the directly surrounding non-active atoms"
                  " must be available. For that reason, the active atoms are extended by a shell of surrounding non-active"
                  " atoms which are also included in the geometry optimization, but whose positions are constrained. This"
                  " shell of non-active atoms can be automatically chosen by ORCA (see section 8.13.1 in the manual). This"
                  " code includes two options for choosing the extension shell. The first is the default option of ORCA"
                  " which includes those non-active atoms in the extension shell that have a distance of less than the sum"
                  " of their VDW radii plus Dist AtomsAroundOpt. The second option is to create a shell of atoms in the"
                  " non-active region at a specific distance from the active region. The user will be asked to provide the"
                  " distance.")
            # Ask the user to choose the atoms belonging to the extension shell
            Ex_sh_options =  input("Select an option to choose the atoms belonging to optRegionExt: " + "\n"
                                   "Option 1       - Do not use extended active region" + "\n"
                                   "Option 2       - Add only atoms bonded covalently to active atoms" + "\n"
                                   "Option 3       - (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1))"+ "\n"
                                   "Option 4       - Manually (using our code) define the extended active region (Default: empty list)" + "\n"
                                   "Enter your option: ")
           		    			
	    #########################################################################
            
            # ezafe shode

            #########################################################################
            # option 1
            if Ex_sh_options == "1":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmqm2mm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmqm2mm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    functional = "TPSSh"
                    basis = "def2-SVP"
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                    print("The default QM2 method is XTB")
                method = input("Do you accept the QM2 method (yes/no): ").lower()
                if method == "" or method == "yes":
                    qmmm_file.write("!QM/XTB/MM \n")
                    qmmm_file.write("%qmmm \n")
                    topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                    if topology_file == "":
                        topology_file = "prmtop.ORCAFF.prms"
                    else:
                        topology_file == topology_file
                    qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                    qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                    
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                    
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                    qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"{mm_lsOrca} \n")  
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"{mm_lsOrca_new} \n")  
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                    qmmm_file.write("Use_QM_InfoFromPDB false \n")
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"{q2mm_lsOrca} \n")  
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                    qmmm_file.write("Use_Active_InfoFromPDB false \n")
                    qmmm_file.write(f"# Option 1: no extension shell will be included in the QM/MM calculations \n")
                    qmmm_file.write("ExtendActiveRegion No \n")                
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()
            #######################    
                elif method == "no":
                    print("You have not accepted the default QM2 method")
                    Q2_method = input("Enter the name of your desired QM2 method: ")
                    if Q2_method == "QM2" or Q2_method == "qm2":
                        print("You should determine your functional and basis set")
                        functional_q2 = input("Enter the name of your desired functional: ")
                        basis_q2 = input("Enter the name of your desired basis set: ")
                        qmmm_file.write(f"!QM/QM2/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")
                        qmmm_file.write(f'QM2CUSTOMMETHOD "{functional_q2}" \n')
                        qmmm_file.write(f'QM2CUSTOMBASIS  "{basis_q2} def2/J"\n')
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 1: no extension shell will be included in the QM/MM calculations \n")
                        qmmm_file.write("ExtendActiveRegion No \n")                
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()
        
                        #####################
                    else:
                        qmmm_file.write(f"!QM/{Q2_method}/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                                            
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")  
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 1: no extension shell will be included in the QM/MM calculations \n")
                        qmmm_file.write("ExtendActiveRegion No \n")                
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()
            #########################################################################
            # option 2
            elif Ex_sh_options == "2":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmqm2mm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmqm2mm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    functional = "TPSSh"
                    basis = "def2-SVP"
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                    print("The default QM2 method is XTB")
                method = input("Do you accept the QM2 method (yes/no): ").lower()
                if method == "" or method == "yes":
                    qmmm_file.write("!QM/XTB/MM \n")
                    qmmm_file.write("%qmmm \n")
                    topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                    if topology_file == "":
                        topology_file = "prmtop.ORCAFF.prms"
                    else:
                        topology_file == topology_file
                    qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                    qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")

                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                    
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                    qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"{mm_lsOrca} \n")  
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"{mm_lsOrca_new} \n")  
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                    qmmm_file.write("Use_QM_InfoFromPDB false \n")
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"{q2mm_lsOrca} \n")  
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                    qmmm_file.write("Use_Active_InfoFromPDB false \n")
                    qmmm_file.write(f"# Option 2: add only atoms bonded covalently to active atoms \n")
                    qmmm_file.write("ExtendActiveRegion cov_bonds \n")            
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()
            #######################    
                elif method == "no":
                    print("You have not accepted the default QM2 method")
                    Q2_method = input("Enter the name of your desired QM2 method: ")
                    if Q2_method == "QM2" or Q2_method == "qm2":
                        print("You should determine your functional and basis set")
                        functional_q2 = input("Enter the name of your desired functional: ")
                        basis_q2 = input("Enter the name of your desired basis set: ")
                        qmmm_file.write(f"!QM/QM2/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")
                        qmmm_file.write(f'QM2CUSTOMMETHOD "{functional_q2}" \n')
                        qmmm_file.write(f'QM2CUSTOMBASIS  "{basis_q2} def2/J"\n')
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 2: add only atoms bonded covalently to active atoms \n")
                        qmmm_file.write("ExtendActiveRegion cov_bonds \n")                
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()
        
                        #####################
                    else:
                        qmmm_file.write(f"!QM/{Q2_method}/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")  
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 2: add only atoms bonded covalently to active atoms \n")
                        qmmm_file.write("ExtendActiveRegion cov_bonds \n")                 
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()

            ##################################################################################################
            # option 3
            elif Ex_sh_options == "3":
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmqm2mm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmqm2mm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    functional = "TPSSh"
                    basis = "def2-SVP"
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                    print("The default QM2 method is XTB")
                method = input("Do you accept the QM2 method (yes/no): ").lower()
                if method == "" or method == "yes":
                    qmmm_file.write("!QM/XTB/MM \n")
                    qmmm_file.write("%qmmm \n")
                    topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                    if topology_file == "":
                        topology_file = "prmtop.ORCAFF.prms"
                    else:
                        topology_file == topology_file
                    qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                    qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")

                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                    
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                    qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"{mm_lsOrca} \n")  
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"{mm_lsOrca_new} \n")  
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                    qmmm_file.write("Use_QM_InfoFromPDB false \n")
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"{q2mm_lsOrca} \n")  
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                    qmmm_file.write("Use_Active_InfoFromPDB false \n")
                    qmmm_file.write(f"# Option 3: (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) \n")
                    qmmm_file.write("ExtendActiveRegion distance \n")
                    distance_criterion = input("Enter the Dist_AtomsAroundOpt value (Default; 1): ")
                    if distance_criterion == "":
                        qmmm_file.write("Dist_AtomsAroundOpt 1.\n")
                    else:
                        qmmm_file.write(f"Dist_AtomsAroundOpt {distance_criterion}\n")
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()
            #######################    
                elif method == "no":
                    print("You have not accepted the default QM2 method")
                    Q2_method = input("Enter the name of your desired QM2 method: ")
                    if Q2_method == "QM2" or Q2_method == "qm2":
                        print("You should determine your functional and basis set")
                        functional_q2 = input("Enter the name of your desired functional: ")
                        basis_q2 = input("Enter the name of your desired basis set: ")
                        qmmm_file.write(f"!QM/QM2/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")
                        qmmm_file.write(f'QM2CUSTOMMETHOD "{functional_q2}" \n')
                        qmmm_file.write(f'QM2CUSTOMBASIS  "{basis_q2} def2/J"\n')
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 3: (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) \n")
                        qmmm_file.write("ExtendActiveRegion distance \n")
                        distance_criterion = input("Enter the Dist_AtomsAroundOpt value (Default; 1): ")
                        if distance_criterion == "":
                            qmmm_file.write("Dist_AtomsAroundOpt 1.\n")
                        else:
                            qmmm_file.write(f"Dist_AtomsAroundOpt {distance_criterion}\n")               
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()

                    #####################
                    else:
                        qmmm_file.write(f"!QM/{Q2_method}/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")  
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"# Option 3: (ORCA Default) Use a distance criterion (VDW distance plus Dist_AtomsAroundOpt (ORCA Default 1)) \n")
                        qmmm_file.write("ExtendActiveRegion distance \n")
                        distance_criterion = input("Enter the Dist_AtomsAroundOpt value (Default; 1): ")
                        if distance_criterion == "":
                            qmmm_file.write("Dist_AtomsAroundOpt 1.\n")
                        else:
                            qmmm_file.write(f"Dist_AtomsAroundOpt {distance_criterion}\n")                
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()

            ################################################################################################## 
            elif Ex_sh_options == "4":
                if mm_groups.lower() == "a":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in q2mm_list_sorted:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    #list_resi_num_es = []
                    #list_all_es = []
                    #for line in data_pdb1.split('\n')[:-1]:
                        #if line.startswith("ATOM") or line.startswith("HETATM"):
                            #if int(line[6:11].strip()) in es_list:
                                #atm_num = line[22:26].strip()
                                #list_resi_num_es.append(atm_num)
                                #for line in data_pdb1.split('\n')[:-1]:
                                    #if line[22:26].strip() == atm_num:
                                        #list_all_es.append(line[6:11].strip())
                                    #else:
                                        #continue
                            #else:
                                #continue
                        #else:
                            #continue
    
                    #new_list_all_es = []
                    #for one_mm in list_all_es:
                        #if int(one_mm) not in new_list_all_es:
                            #new_list_all_es.append(int(one_mm))
                        #else:
                            #continue
                    es_list_new = set(es_list) - set(q2mm_list_sorted)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                elif atoms_or_groups.lower() == "g":
                    # Ask the user for the thickness of the extension shell
                    extension_shell_thickness = input("Enter the desired thickness (in Angstroms) for the MM-fixed region surrounding the MM-active region (default is 20.0): ")
                    if extension_shell_thickness == "":
                        extension_shell_thickness = float(20.0)
                    else:
                        extension_shell_thickness = float(extension_shell_thickness)
                    # mm_list_sorted = [1,2]
                    # Define the variables
                    x_ar = 0.00
                    y_ar = 0.00
                    z_ar = 0.00
                    x_es = 0.00
                    y_es = 0.00
                    z_es = 0.00
                    SD_es = 0.00
                    # Define the new list
                    es_list_ini = []
                    for item in new_list_all_mmq2:
                        for line in data_pdb1.split('\n')[:-1]:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                if int(item) == int(line[6:11]):
                                    x_ar = float(line[30:38])
                                    y_ar = float(line[38:46])
                                    z_ar = float(line[46:54])
                                    for lines in data_pdb1.split('\n')[:-1]:
                                        if lines.startswith("ATOM") or lines.startswith("HETATM"):
                                            x_es = float(lines[30:38])
                                            y_es = float(lines[38:46])
                                            z_es = float(lines[46:54])
                                            SD_es = ((x_es - x_ar) ** 2 + (y_es - y_ar) ** 2 + (z_es - z_ar) ** 2) ** 0.5
                                            if SD_es <= extension_shell_thickness:
                                                mm = int(lines[6:11])
                                                es_list_ini.append(mm)
                    # Remove duplicates from the MM list
                    es_list = list(dict.fromkeys(es_list_ini))
                    ###############################
                    list_resi_num_es = []
                    list_all_es = []
                    for line in data_pdb1.split('\n')[:-1]:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if int(line[6:11].strip()) in es_list:
                                atm_num = line[20:26].strip()
                                list_resi_num_es.append(atm_num)
                                for line in data_pdb1.split('\n')[:-1]:
                                    if line[20:26].strip() == atm_num:
                                        list_all_es.append(line[6:11].strip())
                                    else:
                                        continue
                            else:
                                continue
                        else:
                            continue
    
                    new_list_all_es = []
                    for one_mm in list_all_es:
                        if int(one_mm) not in new_list_all_es:
                            new_list_all_es.append(int(one_mm))
                        else:
                            continue
                    #es_list_new = set(new_list_all_es)                
                    es_list_new = set(new_list_all_es) - set(new_list_all_mmq2)
                    # Sort the MM atoms list
                    es_list_sorted_PDB_format = sorted(list(set(es_list_new)))
                    ES_list_ORCA_FORMAT = " ".join(str(es_atom - 1) for es_atom in es_list_sorted_PDB_format)
                    # print(sorted(list(set(mm_list))))
                    es_lst = ['OptRegion_FixedAtoms ', '{', ES_list_ORCA_FORMAT, '}', ' end']
                    es_lst2 = ''.join(map(str, es_lst))
                    print(f"The following {len(es_list_sorted_PDB_format)} atoms are included in the extension shell.")
                    print(es_list_sorted_PDB_format)
                    # simplify the new-list_es
                    eslistoc_new = []
                    for i in range(0,len(es_list_sorted_PDB_format)):
                        if i <= 0:
                            for j in range (0,len(es_list_sorted_PDB_format)):
                                if es_list_sorted_PDB_format [i] in range(es_list_sorted_PDB_format[i] and es_list_sorted_PDB_format[j]+1) and es_list_sorted_PDB_format[j]+1 not in es_list_sorted_PDB_format and es_list_sorted_PDB_format[i]-1 not in es_list_sorted_PDB_format:
                                    if es_list_sorted_PDB_format[i] != min(es_list_sorted_PDB_format) or es_list_sorted_PDB_format[j] != max(es_list_sorted_PDB_format):
                                        if es_list_sorted_PDB_format[i] != es_list_sorted_PDB_format[j]:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1)+":"+str(es_list_sorted_PDB_format[j]-1))
                                            #print(str(es_list_sorted_PDB_format[i])+":"+str(es_list_sorted_PDB_format[j]))
                                            i = j+1
                                        else:
                                            eslistoc_new.append(str(es_list_sorted_PDB_format[i]-1))
                                            #print(str(es_list_sorted_PDB_format[i]))
                                            i = j+1
                            eslsf_new=str(eslistoc_new).replace('[','')
                            eslsff_new=str(eslsf_new).replace(']','')
                            eslsfff_new=str(eslsff_new).replace("'","")
                            eslsffff_new=str(eslsfff_new).replace(",","")
                            es_lsOrca_new = 'OptRegion_FixedAtoms'+' {'+ eslsffff_new+'} '+ 'end'
                # Assign the output PDB file for QM/MM calculations by ORCA
                #out_pdb_for_ORCA_QMMM = input(
                    #"Enter the name of your output PDB file, to assign the QM atoms and active region in its "
                    #"occupancy and B-factor columns, respectively. (out_pdb_ORCA_QMMM.pdb): ")
                #if out_pdb_for_ORCA_QMMM == "":
                    #out_pdb_for_ORCA_QMMM = "out_pdb_ORCA_QMMM.pdb"
                #else:
                    #out_pdb_for_ORCA_QMMM == out_pdb_for_ORCA_QMMM
                # We expect from the following lines to read the input PDB file (defined in the beginig of the code), line by line.
                # Read the mm_list_sorted list item by item.
                # If indexes 6–11 in a line of the input pdb file are equal to the item, replace indexes 64–67 of that line with 1.00.
                #lower = 0
                #upper = len(mm_list_sorted)
                #low_pdb = 0
                #line_in_input_pdb = data_pdb1.split("\n")[:-1]
                #up_pdb = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(line_in_input_pdb)
                #low_qm = 0
                #up_qm = len(Linear_QM_Atoms_List_PDB_FORMAT_sorted)
                #with open(out_pdb_for_ORCA_QMMM, "w") as output:
                    #for mm_atom in range(lower, upper):
                        #if mm_atom <= lower:
                            #mm_list_sorted.append("")
                            #for qm_atom in range(low_qm,up_qm):
                                #if qm_atom <= low_qm:
                                    #Linear_QM_Atoms_List_PDB_FORMAT_sorted.append("")
                                    #for line_pdb in range (low_pdb, up_pdb):
                                        #if str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip()and str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + "  1.00"+ line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(Linear_QM_Atoms_List_PDB_FORMAT_sorted[qm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:57] + "1.00" + line_in_input_pdb[line_pdb][61:] + "\n"
                                            #mm_atom += 1
                                            #qm_atom += 1
                                        #elif str(mm_list_sorted[mm_atom]) == line_in_input_pdb[line_pdb][6:11].strip():
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb][0:63] + "1.00" + line_in_input_pdb[line_pdb][67:] + "\n"
                                            #mm_atom += 1
                                        #else:
                                            #line_in_output_pdb = line_in_input_pdb[line_pdb] + "\n"
                                        #output.write(line_in_output_pdb)
                    #output.write("END )
                #########################################################################
                print(
                    f"The PDB file containing the QM1, QM2, MM active,and -fixed  regions (QM1-QM2-MMactive-MMfixed.pdb) is generated.")            
                pdb_text_14 = open(pdb_file_1, "r")
                data_pdb14 = pdb_text_14.readlines()
                if mm_groups.lower() == "a":
                    totpdb = sorted(q2mm_list_sorted + es_list_sorted_PDB_format)
                    with open ("QM1-QM2-MMactive-MMfixed.pdb", "w") as output:
                        for Pline in data_pdb14:
                            if int(Pline[6:11].strip()) in totpdb:
                                output.write(Pline)
                if mm_groups.lower() == "g":
                    totpdb = sorted(new_list_all_mmq2 + es_list_sorted_PDB_format)
                    with open ("QM1-QM2-MMactive-MMfixed.pdb", "w") as output:
                        for Pline in data_pdb14:
                            if int(Pline[6:11].strip()) in totpdb:
                                output.write(Pline)           
                #########################################################################
                qmmm_file = input("Enter the name of ORCA input file for QM/MM calculations (qmqm2mm.inp): ")
                if qmmm_file == "":
                    qmmm_file = "qmqm2mm.inp"
                else:
                    qmmm_file == qmmm_file
                qmmm_file = open(qmmm_file, "w+")
                qmmm_file.write("# ORCA's input file for performing QM/MM calculations. \n")
                qmmm_file.write("# This file is generated by a home-made Python program (pdbtoorca) by Mehdi and Maryam. \n")
                qmmm_file.write("# The program is for free. However, cite it as follows, if you use it in your project. \n")
                qmmm_file.write("# DOI: ................ \n")
                # The working directory and the SetUp time
                import os
                WorDir = os.getcwd()
                qmmm_file.write("# The working directory is: " +str(WorDir)+ "\n")
                from datetime import date,time, datetime
                today = date.today()
                now = datetime. now()
                current_time = time(now.hour, now.minute, now.second)
                qmmm_file.write("# The setup time and date are " +str(current_time)+ " and " +str(today)+ ", respectively. \n")
                print("We are setting up the level of of theory.")
                print("The default level of theory is:")
                print("RIJCOSX TPSSh def2-SVP def2/J D3BJ TIGHTSCF OPT")
                theory = input("Do you accept the default level of theory (yes/no): ").lower()
                if theory == "" or theory == "yes":
                    functional = "TPSSh"
                    basis = "def2-SVP"
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                elif theory == "no":
                    print("You have not accepted the default level of theory")
                    functional = input("Enter the name of your desired functional: ")
                    basis = input("Enter the name of your desired basis set: ")
                    qmmm_file.write(f"! RIJCOSX {functional} {basis} def2/J D3BJ TIGHTSCF Opt \n")
                    print("The default QM2 method is XTB")
                method = input("Do you accept the QM2 method (yes/no): ").lower()
                if method == "" or method == "yes":
                    qmmm_file.write("!QM/XTB/MM \n")
                    qmmm_file.write("%qmmm \n")
                    topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                    if topology_file == "":
                        topology_file = "prmtop.ORCAFF.prms"
                    else:
                        topology_file == topology_file
                    qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                    qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")

                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                    
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                    qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                    qmmm_file.write(f"# A total of {len(es_list_sorted_PDB_format)} atoms are included in the extension_shell region. \n")
                    qmmm_file.write(f"# The thickness of the extension_shell region is {extension_shell_thickness} angstroms. \n")                    
                    qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                    if atoms_or_groups.lower() == "a":
                        qmmm_file.write(f"{mm_lsOrca} \n")  
                    elif atoms_or_groups.lower() == "g":
                        qmmm_file.write(f"{mm_lsOrca_new} \n")
                    # Ask the user for the charge of the medium (QM1+QM2) system
                    charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                    if charge_of_the_medium_region == "":
                        charge_of_the_medium_region = int(0)
                    else:
                        charge_of_the_medium_region = int(charge_of_the_medium_region)
                    # Ask the user for the multiplicity of the medium (QM1+QM2) system
                    multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                    if multiplicity_of_the_medium_region == "":
                        multiplicity_of_the_medium_region = int(1)
                    else:
                        multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                    qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                    qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                    qmmm_file.write("Use_QM_InfoFromPDB false \n")
                    if mm_groups.lower() == "a":
                        qmmm_file.write(f"{q2mm_lsOrca} \n")  
                    elif mm_groups.lower() == "g":
                        qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                    qmmm_file.write("Use_Active_InfoFromPDB false \n")
                    qmmm_file.write(f"{es_lsOrca_new} \n")
                    qmmm_file.write("end \n")
                    geom = input("Do you want to do geom scan (yes/no): ").lower()
                    if geom == "y" or geom == "yes":
                        qmmm_file.write(f"% geom scan \n")  
                        atom1 = input("Enter the number of the first atom: ").lower()
                        orca_atom1 = str(int(atom1)-1)
                        atom2 = input("Enter the number of the sencond atom: ").lower()
                        orca_atom2 = str(int(atom2)-1)
                        dist1= input("Enter the initial distance: ").lower()
                        dist2 = input("Enter the final distance: ").lower()
                        num = input("Enter the number of scan proceses: ").lower()
                        qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                        qmmm_file.write("end \n")
                        qmmm_file.write("end \n")
                    # Ask the user for the number of CPU cores for parallel calculations
                    cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                    if cpu_cores == "":
                        cpu_cores = int(28)
                    else:
                        cpu_cores = int(cpu_cores)
                    qmmm_file.write("% pal \n")
                    qmmm_file.write("nprocs \n")
                    qmmm_file.write(f"{cpu_cores} \n")
                    qmmm_file.write("end \n")
                    # Ask the user for the charge of the active region
                    charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                    if charge_of_the_QM_region == "":
                        charge_of_the_QM_region = int(0)
                    else:
                        charge_of_the_QM_region = int(charge_of_the_QM_region)
                    # Ask the user for the multiplicity of the active region
                    multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                    if multiplicity_of_the_QM_region == "":
                        multiplicity_of_the_QM_region = int(1)
                    else:
                        multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                    #pdb_file_orca = 'orcapdb.pdb'
                    qmmm_file.write(
                        f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                    qmmm_file.write(" \n")
                    qmmm_file.close()
            #######################    
                elif method == "no":
                    print("You have not accepted the default QM2 method")
                    Q2_method = input("Enter the name of your desired QM2 method: ")
                    if Q2_method == "QM2" or Q2_method == "qm2":
                        print("You should determine your functional and basis set")
                        functional_q2 = input("Enter the name of your desired functional: ")
                        basis_q2 = input("Enter the name of your desired basis set: ")
                        qmmm_file.write(f"!QM/QM2/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"# A total of {len(es_list_sorted_PDB_format)} atoms are included in the extension_shell region. \n")
                        qmmm_file.write(f"# The thickness of the extension_shell region is {extension_shell_thickness} angstroms. \n") 
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")
                        qmmm_file.write(f'QM2CUSTOMMETHOD "{functional_q2}" \n')
                        qmmm_file.write(f'QM2CUSTOMBASIS  "{basis_q2} def2/J"\n')
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"{es_lsOrca_new} \n")             
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()

                    #####################
                    else:
                        qmmm_file.write(f"!QM/{Q2_method}/MM \n")
                        qmmm_file.write("%qmmm \n")
                        topology_file = input("Enter the name of topology file (prmtop.ORCAFF.prms): ")
                        if topology_file == "":
                            topology_file = "prmtop.ORCAFF.prms"
                        else:
                            topology_file == topology_file
                        qmmm_file.write(f'ORCAFFFilename "{topology_file}" \n')                
                        qmmm_file.write(f"# Excluding junctions, {len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM region. \n")
        
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(mm_list_sorted) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mm) - len(Linear_QM_Atoms_List_PDB_FORMAT)} atoms are included in the QM2 region. \n")         
                        
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"# A total of {len(q2mm_list_sorted)} atoms are included in the active region. \n")
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"# A total of {len(new_list_all_mmq2)} atoms are included in the active region. \n")
                        qmmm_file.write(f"# The thickness of the active region is {mm_region_thickness} angstroms (from the QM2 region). \n")
                        qmmm_file.write(f"# A total of {len(es_list_sorted_PDB_format)} atoms are included in the extension_shell region. \n")
                        qmmm_file.write(f"# The thickness of the extension_shell region is {extension_shell_thickness} angstroms. \n") 
                        qmmm_file.write(f"{QM_atoms_list_ORCA_format} \n")
                        if atoms_or_groups.lower() == "a":
                            qmmm_file.write(f"{mm_lsOrca} \n")  
                        elif atoms_or_groups.lower() == "g":
                            qmmm_file.write(f"{mm_lsOrca_new} \n")  
                        # Ask the user for the charge of the medium (QM1+QM2) system
                        charge_of_the_medium_region = input("Enter the charge of the medium (QM1+QM2) system (0): ")
                        if charge_of_the_medium_region == "":
                            charge_of_the_medium_region = int(0)
                        else:
                            charge_of_the_medium_region = int(charge_of_the_medium_region)
                        # Ask the user for the multiplicity of the medium (QM1+QM2) system
                        multiplicity_of_the_medium_region = input("Enter the multiplicity of the medium (QM1+QM2) system (1): ")
                        if multiplicity_of_the_medium_region == "":
                            multiplicity_of_the_medium_region = int(1)
                        else:
                            multiplicity_of_the_medium_region = int(multiplicity_of_the_medium_region)          
                        qmmm_file.write(f"Charge_Medium {charge_of_the_medium_region}   # Charge of medium (QM1+QM2) system.  \n")
                        qmmm_file.write(f"Mult_Medium   {multiplicity_of_the_medium_region}   # multiplicity of the medium system. Default 1.  \n")
                        qmmm_file.write("Use_QM_InfoFromPDB false \n")
                        if mm_groups.lower() == "a":
                            qmmm_file.write(f"{q2mm_lsOrca} \n")  
                        elif mm_groups.lower() == "g":
                            qmmm_file.write(f"{q2mm_lsOrca_new} \n")                
                        qmmm_file.write("Use_Active_InfoFromPDB false \n")
                        qmmm_file.write(f"{es_lsOrca_new} \n")          
                        qmmm_file.write("end \n")
                        geom = input("Do you want to do geom scan (yes/no): ").lower()
                        if geom == "y" or geom == "yes":
                            qmmm_file.write(f"% geom scan \n")  
                            atom1 = input("Enter the number of the first atom: ").lower()
                            orca_atom1 = str(int(atom1)-1)
                            atom2 = input("Enter the number of the sencond atom: ").lower()
                            orca_atom2 = str(int(atom2)-1)
                            dist1= input("Enter the initial distance: ").lower()
                            dist2 = input("Enter the final distance: ").lower()
                            num = input("Enter the number of scan proceses: ").lower()
                            qmmm_file.write(f"B {orca_atom1} {orca_atom2} = {dist1}, {dist2}, {num} \n")
                            qmmm_file.write("end \n")
                            qmmm_file.write("end \n")
                        # Ask the user for the number of CPU cores for parallel calculations
                        cpu_cores = input("Enter the number of CPU cores for parallel calculations (28): ")
                        if cpu_cores == "":
                            cpu_cores = int(28)
                        else:
                            cpu_cores = int(cpu_cores)
                        qmmm_file.write("% pal \n")
                        qmmm_file.write("nprocs \n")
                        qmmm_file.write(f"{cpu_cores} \n")
                        qmmm_file.write("end \n")
                        # Ask the user for the charge of the active region
                        charge_of_the_QM_region = input("Enter the charge of the QM region (0): ")
                        if charge_of_the_QM_region == "":
                            charge_of_the_QM_region = int(0)
                        else:
                            charge_of_the_QM_region = int(charge_of_the_QM_region)
                        # Ask the user for the multiplicity of the active region
                        multiplicity_of_the_QM_region = input("Enter the multiplicity of the QM region (1): ")
                        if multiplicity_of_the_QM_region == "":
                            multiplicity_of_the_QM_region = int(1)
                        else:
                            multiplicity_of_the_QM_region = int(multiplicity_of_the_QM_region)
                        #pdb_file_orca = 'orcapdb.pdb'
                        qmmm_file.write(
                            f"*pdbfile {charge_of_the_QM_region} {multiplicity_of_the_QM_region} {pdb_file_1}  # charge and mult. of the high level region\n")
                        qmmm_file.write(" \n")
                        qmmm_file.close()
            print("All necessary files for multi-scale calculations have been generated in your working directory. You can now execute another command or exit the program.")
             
            ######################################
    else:
        print("I do not understand you")
# End Time
# from datetime import datetime
end_time = datetime.now()
duration = end_time - start_time
print("""
""")
print(f"The job started at: {start_time}")
print(f"The job ended   at: {end_time}")
print(f"Time duration   is: {duration.seconds} Seconds")
