#!/usr/bin/env python3

import argparse
import tkinter as tk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo
import math
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors

class GuidemolGUI:
    """A class to represent the GUI of the app
    Consists of a canvas + scrollbar
    The canvas includes a main frame (frameM) with:
        - frame_gen3d: left frame concerning the generation of 3D sructures
        - frame_desc: first right frame for the calculation of RDKit
        descriptors
        - frame_grid_potential: second right frame for the calculation of a
        grid of electrostatic potential (EP)
        - frame_grid_vox: third right frame for the calculation of a grid of
        voxels
    """

    def __init__(self,window):
        """Configure the main window,
           Create:
            - the scrollbar, frames, buttons, labels
            Place all widgets in the GUI
        """
        self.curr_mol_text = tk.StringVar()
        self.addHs_flag = tk.IntVar()
        self.window = window
        self.window.geometry("620x700")
        self.window.title(
            "GUIDEMOL: GUI for molecular processing and 3D descriptors with RDKit"
        )
        scrollbar = tk.Scrollbar(orient="vertical")
        scrollbar.grid(column=1, row=0, sticky="ns")
        canvas = tk.Canvas(yscrollcommand=scrollbar.set)
        scrollbar.config(command=canvas.yview)
        frameM = tk.Frame(canvas)
        frame_gen3d=self.create_frame_gen3d(frameM)
        frame_gen3d.grid(row=0, column=0, sticky="ns", rowspan=8)
        lbl_title_desc = tk.Label(
                                frameM,
                                text="Molecular descriptors",
                                background='white'
                                )
        lbl_title_desc.grid(row=0, column=1, pady=15, padx=5, sticky="ew")
        frame_desc=self.create_frame_desc(frameM)
        frame_desc.grid(row=1, column=1, sticky="ew", padx=5)
        frame_grid_potential=self.create_frame_grid_potential(frameM)
        frame_grid_potential.grid(row=2, column=1, sticky="ew", padx=5)
        frame_grid_vox=self.create_frame_grid_vox(frameM)
        frame_grid_vox.grid(row=3, column=1, sticky="ew", padx=5)
        lbl_curr_mol = tk.Label(frameM, textvariable=self.curr_mol_text)
        chk_addHs = tk.Checkbutton(
                        frameM,
                        text="Add hydrogens",
                        variable=self.addHs_flag
                        )
        self.addHs_flag.set(0)
        chk_addHs.grid(row=4, column=1, padx=10, pady=5, sticky="ew")
        btn_calcdesc = tk.Button(
                                frameM,
                                text=" Calculate descriptors ",
                                command=self.calcdesc
                                )
        btn_calcdesc.grid(row=6, column=1, padx=10, pady=5, sticky="ew")
        lbl_curr_mol.grid(row=7, column=1, padx=5, pady=0, sticky="ew")
        canvas.create_window(0,0, window=frameM, anchor="nw")
        canvas.grid(column=0, row=0, sticky="nswe")
        canvas.configure(scrollregion = (0,0,620,750))
        window.columnconfigure(0, weight=1)
        window.rowconfigure(0, weight=1)

    def create_frame_desc(self, window):
        """Create the frame concerning the calculation of RDKit descriptors:
            - create the checkbuttons and the associated variables
            - place the widgets on the frame (frameDesc)
            - return the frame (frameDesc)
           Parameter: the frame to place the widgets
        """
        frameDesc = tk.Frame(
                        window,
                        bg='LightCyan2',
                        borderwidth=1,
                        relief='raised'
                        )
        self.RDF_flag = tk.IntVar()
        chk_RDF = tk.Checkbutton(
                        frameDesc,text="RDF",
                        variable=self.RDF_flag,
                        bg='LightCyan2'
                        )
        self.RDF_flag.set(0)
        self.MORSE_flag = tk.IntVar()
        chk_MORSE = tk.Checkbutton(
                        frameDesc,
                        text="MORSE",
                        variable=self.MORSE_flag,
                        bg='LightCyan2'
                        )
        self.MORSE_flag.set(0)
        self.WHIM_flag = tk.IntVar()
        chk_WHIM = tk.Checkbutton(
                        frameDesc,
                        text="WHIM",
                        variable=self.WHIM_flag,
                        bg='LightCyan2'
                        )
        self.WHIM_flag.set(0)
        self.AUTOCORR3D_flag = tk.IntVar()
        chk_AUTOCORR3D = tk.Checkbutton(
                                frameDesc,
                                text="AUTOCORR3D",
                                variable=self.AUTOCORR3D_flag,
                                bg='LightCyan2'
                                )
        self.AUTOCORR3D_flag.set(0)
        self.GETAWAY_flag = tk.IntVar()
        chk_GETAWAY = tk.Checkbutton(
                        frameDesc,
                        text="GETAWAY",
                        variable=self.GETAWAY_flag,
                        bg='LightCyan2'
                        )
        self.GETAWAY_flag.set(0)
        self.PEOE_flag = tk.IntVar()
        chk_PEOE = tk.Checkbutton(
                        frameDesc,
                        text="PEOE",
                        variable=self.PEOE_flag,
                        bg='LightCyan2'
                        )
        self.PEOE_flag.set(0)
        self.SMR_flag = tk.IntVar()
        chk_SMR = tk.Checkbutton(
                        frameDesc,
                        text="SMR",
                        variable=self.SMR_flag,
                        bg='LightCyan2'
                        )
        self.SMR_flag.set(0)
        self.MW_flag = tk.IntVar()
        chk_MW = tk.Checkbutton(
                        frameDesc,
                        text="MW",
                        variable=self.MW_flag,
                        bg='LightCyan2'
                        )
        self.MW_flag.set(0)
        chk_RDF.grid(row=2, column=0, padx=10, pady=5, sticky="w")
        chk_MORSE.grid(row=2, column=1, padx=0, pady=5, sticky="w")
        chk_WHIM.grid(row=3, column=0, padx=10, pady=5, sticky="w")
        chk_AUTOCORR3D.grid(row=3, column=1, padx=0, pady=5, sticky="w")
        chk_GETAWAY.grid(row=4, column=0, padx=10, pady=5, sticky="w")
        chk_PEOE.grid(row=4, column=1, padx=0, pady=5, sticky="w")
        chk_SMR.grid(row=5, column=0, padx=10, pady=5, sticky="w")
        chk_MW.grid(row=5, column=1, padx=0, pady=5, sticky="w")
        return frameDesc

    def create_frame_gen3d(self,window):
        """Create the frame concerning the generation of 3D structures:
            - create the checkbuttons, the associated variables and a button
            - place the widgets on the frame (frameGen3D)
            - return the frame (frameGen3D)
           Parameter: the frame to place the widgets
        """
        self.Metals_flag = tk.IntVar()
        self.Fragment_flag = tk.IntVar()
        self.Uncharge_flag = tk.IntVar()
        self.Optimize_flag = tk.IntVar()
        frameGen3D = tk.Frame(
                        window, bg='LightBlue', borderwidth=1, relief='raised'
                        )
        lbl_title_gen3D = tk.Label(
                                frameGen3D,
                                text="Get 3D structures",
                                background='white'
                                )
        chk_Metals = tk.Checkbutton(
                                frameGen3D,
                                text="Remove metals",
                                variable=self.Metals_flag,
                                bg='LightBlue'
                                )
        self.Metals_flag.set(1)
        chk_Fragment = tk.Checkbutton(
                                frameGen3D,
                                text="Remove smaller fragments",
                                variable=self.Fragment_flag,
                                bg='LightBlue'
                                )
        self.Fragment_flag.set(1)
        chk_Uncharge = tk.Checkbutton(
                                frameGen3D,
                                text="Uncharge",
                                variable=self.Uncharge_flag,
                                bg='LightBlue'
                                )
        self.Uncharge_flag.set(1)
        chk_Optimize = tk.Checkbutton(
                                frameGen3D,
                                text="Optimize 3D with MMFF94",
                                variable=self.Optimize_flag,
                                bg='LightBlue'
                                )
        self.Optimize_flag.set(1)
        btn_gen3D = tk.Button(
                        frameGen3D,
                        text=" 3D from SMILES ",
                        command=self.gen3DfromSMILES
                        )
        lbl_title_gen3D.grid(row=0, column=0, pady=15, sticky="ew")
        chk_Metals.grid(row=2, column=0, padx=10, pady=5, sticky="w")
        chk_Fragment.grid(row=3, column=0, padx=10, pady=5, sticky="w")
        chk_Uncharge.grid(row=4, column=0, padx=10, pady=5, sticky="w")
        chk_Optimize.grid(row=5, column=0, padx=10, pady=5, sticky="w")
        btn_gen3D.grid(row=6, column=0, padx=10, pady=15)
        return frameGen3D

    def create_frame_grid_potential(self, window):
        """Create the frame concerning the grid of EP:
            - create the checkbuttons, radiobuttons, labels, entries
              and the associated variables (to specify the grid)
            - place the widgets on the frame (frameGrid)
            - return the frame (frameGrid)
           Parameter: the frame to place the widgets
        """
        frameGrid = tk.Frame(
                     window, bg='LightCyan2', borderwidth=1, relief='raised'
                     )
        self.VGrid_flag = tk.IntVar()
        chk_PG = tk.Checkbutton(
                        frameGrid,
                        text="Grid of potential",
                        variable=self.VGrid_flag,
                        bg='LightCyan2'
                        )
        self.VGrid_flag.set(0)
        self.ChargeType = tk.StringVar()
        gasteigerCharge = tk.Radiobutton(
                                frameGrid,
                                text="Gasteiger charges",
                                variable=self.ChargeType,
                                value='Gasteiger',
                                bg='LightCyan2'
                                )
        MMFFCharge = tk.Radiobutton(
                                frameGrid,
                                text="MMFF charges",
                                variable=self.ChargeType,
                                value='MMFF',
                                bg='LightCyan2'
                                )
        self.ChargeType.set('Gasteiger')
        lbl_center = tk.Label(
                        frameGrid,
                        text="Center of grid (x,y,z):",
                        background='LightCyan2'
                        )
        self.centerx_str = tk.StringVar()
        entryx= tk.Entry(
                        frameGrid,
                        textvariable=self.centerx_str,
                        width=6,
                        bg='white'
                        )
        self.centerx_str.set(0)
        self.centery_str = tk.StringVar()
        entryy= tk.Entry(
                        frameGrid,
                        textvariable=self.centery_str,
                        width=6,
                        bg='white'
                        )
        self.centery_str.set(0)
        self.centerz_str = tk.StringVar()
        entryz= tk.Entry(
                        frameGrid,
                        textvariable=self.centerz_str,
                        width=6,
                        bg='white'
                        )
        self.centerz_str.set(0)
        lbl_size = tk.Label(
                        frameGrid,
                        text="Size of grid (x,y,z):",
                        background='LightCyan2'
                        )
        self.csizex_str = tk.StringVar()
        entry_sizex= tk.Entry(
                        frameGrid,
                        textvariable=self.csizex_str,
                        width=6,
                        bg='white'
                        )
        self.csizex_str.set(10)
        self.csizey_str = tk.StringVar()
        entry_sizey= tk.Entry(
                        frameGrid,
                        textvariable=self.csizey_str,
                        width=6,
                        bg='white'
                        )
        self.csizey_str.set(10)
        self.csizez_str = tk.StringVar()
        entry_sizez= tk.Entry(
                        frameGrid,
                        textvariable=self.csizez_str,
                        width=6,
                        bg='white'
                        )
        self.csizez_str.set(10)
        lbl_resol = tk.Label(
                        frameGrid,
                        text="Resolution of grid:",
                        background='LightCyan2'
                        )
        self.cubesize_str = tk.StringVar()
        entry_cubesize = tk.Entry(
                                frameGrid,
                                textvariable=self.cubesize_str,
                                width=6,
                                bg='white'
                                )
        self.cubesize_str.set(1)
        lbl_atmreach = tk.Label(
                        frameGrid,
                        text="Cutoff distance:",
                        background='LightCyan2'
                        )
        self.atmreach_str = tk.StringVar()
        entry_atmreach = tk.Entry(
                                frameGrid,
                                textvariable=self.atmreach_str,
                                width=6,
                                bg='white'
                                )
        self.atmreach_str.set(3.5)
        chk_PG.grid(row=1, column=0, pady=5, sticky="w", columnspan=2)
        gasteigerCharge.grid(row=2, column=1, padx=25, pady=2, sticky="w")
        MMFFCharge.grid(row=3, column=1, padx=25, pady=2, sticky="w")
        lbl_center.grid(row=4, column=1, padx=25, pady=5, sticky="w")
        entryx.grid(row=4, column=2, pady=5, padx=2, sticky="w")
        entryy.grid(row=4, column=3, pady=5, padx=2, sticky="w")
        entryz.grid(row=4, column=4, pady=5, padx=2, sticky="w")
        lbl_size.grid(row=5, column=1, padx=25, pady=5, sticky="w")
        entry_sizex.grid(row=5, column=2, pady=5, padx=2, sticky="w")
        entry_sizey.grid(row=5, column=3, pady=5, padx=2, sticky="w")
        entry_sizez.grid(row=5, column=4, pady=5, padx=2, sticky="w")
        lbl_resol.grid(row=6, column=1, padx=25, pady=5, sticky="w")
        entry_cubesize.grid(row=6, column=2, padx=2, pady=5, sticky="w")
        lbl_atmreach.grid(row=7, column=1, padx=25, pady=5, sticky="w")
        entry_atmreach.grid(row=7, column=2, padx=2, pady=5, sticky="w")
        return frameGrid

    def create_frame_grid_vox(self, window):
        """Create the frame concerning the grid of voxels:
            - create the checkbuttons, radiobuttons, labels, entries
              and the associated variables (to specify the grid)
            - place the widgets on the frame (frameGrid)
            - return the frame (frameVox)
           Parameter: the frame to place the widgets
        """
        frameVox = tk.Frame(
                      window, bg='LightCyan2', borderwidth=1, relief='raised'
                      )
        self.VoxGrid_flag = tk.IntVar()
        chk_VoxGrid = tk.Checkbutton(
                                     frameVox,
                                     text="Grid of voxels",
                                     variable=self.VoxGrid_flag,
                                     bg='LightCyan2'
                                     )
        self.VoxGrid_flag.set(0)
        self.VoxProp = tk.StringVar()
        VoxValenceElectrons = tk.Radiobutton(
                                             frameVox,
                                             text="Atomic numbers",
                                             variable=self.VoxProp,
                                             value='AtomicNumber',
                                             bg='LightCyan2'
                                             )
        VoxMR = tk.Radiobutton(
                               frameVox,
                               text="MR contrib",
                               variable=self.VoxProp,
                               value='MR',
                               bg='LightCyan2'
                               )
        VoxLogp = tk.Radiobutton(
                                 frameVox,
                                 text="LogP contrib",
                                 variable=self.VoxProp,
                                 value='LogP',
                                 bg='LightCyan2'
                                 )
        VoxgasteigerCharge = tk.Radiobutton(
                                            frameVox,
                                            text="Gasteiger charges",
                                            variable=self.VoxProp,
                                            value='Gasteiger',
                                            bg='LightCyan2'
                                            )
        VoxMMFFCharge = tk.Radiobutton(
                                       frameVox,
                                       text="MMFF charges",
                                       variable=self.VoxProp,
                                       value='MMFF',
                                       bg='LightCyan2'
                                       )
        self.VoxProp.set('AtomicNumber')
        lbl_voxcenter = tk.Label(
                                frameVox,
                                text="Center of grid (x,y,z):",
                                background='LightCyan2'
                                )
        self.voxcenterx_str = tk.StringVar()
        voxentryx= tk.Entry(
                        frameVox,
                        textvariable=self.voxcenterx_str,
                        width=6,
                        bg='white'
                        )
        self.voxcenterx_str.set(0)
        self.voxcentery_str = tk.StringVar()
        voxentryy= tk.Entry(
                        frameVox,
                        textvariable=self.voxcentery_str,
                        width=6,
                        bg='white'
                        )
        self.voxcentery_str.set(0)
        self.voxcenterz_str = tk.StringVar()
        voxentryz= tk.Entry(
                        frameVox,
                        textvariable=self.voxcenterz_str,
                        width=6,
                        bg='white'
                        )
        self.voxcenterz_str.set(0)
        lbl_voxsize = tk.Label(
                        frameVox,text="Size of grid (x,y,z):",
                        background='LightCyan2'
                        )
        self.voxcsizex_str = tk.StringVar()
        voxentry_sizex= tk.Entry(
                                frameVox,
                                textvariable=self.voxcsizex_str,
                                width=6,
                                bg='white'
                                )
        self.voxcsizex_str.set(10)
        self.voxcsizey_str = tk.StringVar()
        voxentry_sizey= tk.Entry(
                                frameVox,
                                textvariable=self.voxcsizey_str,
                                width=6,
                                bg='white'
                                )
        self.voxcsizey_str.set(10)
        self.voxcsizez_str = tk.StringVar()
        voxentry_sizez= tk.Entry(
                                frameVox,
                                textvariable=self.voxcsizez_str,
                                width=6,
                                bg='white'
                                )
        self.voxcsizez_str.set(10)
        lbl_voxresol = tk.Label(
                                frameVox,
                                text="Resolution of grid:",
                                background='LightCyan2'
                                )
        self.voxcubesize_str = tk.StringVar()
        voxentry_cubesize = tk.Entry(
                                frameVox,
                                textvariable=self.voxcubesize_str,
                                width=6,
                                bg='white'
                                )
        self.voxcubesize_str.set(1)
        lbl_voxsigma = tk.Label(
                                frameVox,
                                text="Voxel sigma:",
                                background='LightCyan2'
                                )
        self.voxsigma_str = tk.StringVar()
        voxentry_sigma = tk.Entry(
                                frameVox,
                                textvariable=self.voxsigma_str,
                                width=6,
                                bg='white'
                                )
        self.voxsigma_str.set(2)
        chk_VoxGrid.grid(row=1, column=0, pady=5, sticky="w", columnspan=2)
        VoxValenceElectrons.grid(row=2, column=1, padx=25, pady=2, sticky="w")
        VoxMR.grid(row=3, column=1, padx=25, pady=2, sticky="w")
        VoxLogp.grid(row=4, column=1, padx=25, pady=2, sticky="w")
        VoxgasteigerCharge.grid(
                                row=2, column=2, padx=0, pady=2,
                                sticky="w",
                                columnspan=3
                                )
        VoxMMFFCharge.grid(
                        row=3, column=2, padx=0, pady=2,
                        sticky="w",
                        columnspan=3
                        )
        lbl_voxcenter.grid(row=5, column=1, padx=25, pady=5, sticky="w")
        voxentryx.grid(row=5, column=2, padx=2, pady=5, sticky="w")
        voxentryy.grid(row=5, column=3, padx=2, pady=5, sticky="w")
        voxentryz.grid(row=5, column=4, padx=2, pady=5, sticky="w")
        lbl_voxsize.grid(row=6, column=1, padx=25, pady=5, sticky="w")
        voxentry_sizex.grid(row=6, column=2, padx=2, pady=5, sticky="w")
        voxentry_sizey.grid(row=6, column=3, padx=2, pady=5, sticky="w")
        voxentry_sizez.grid(row=6, column=4, padx=2, pady=5, sticky="w")
        lbl_voxresol.grid(row=7, column=1, padx=25, pady=5, sticky="w")
        voxentry_cubesize.grid(row=7, column=2, padx=2, pady=5, sticky="w")
        lbl_voxsigma.grid(row=8, column=1, padx=25, pady=5, sticky="w")
        voxentry_sigma.grid(row=8, column=2, padx=2, pady=5, sticky="w")
        return frameVox

    def gen3DfromSMILES(self):
        """Read molecules in the SMILES format from a file, apply
           transformations on the molecules according to the status of flag
           variables (specified in the GUI), generate the 3D models of
           molecules, save the structures to an MDL SDFile, register messages
           in a log file.
        """
        md = rdMolStandardize.MetalDisconnector()
        lfc = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        u = rdMolStandardize.Uncharger()
        inputfile = fd.askopenfilename(
            title='Choose the input file with the SMILES strings'
        )
        showinfo(title='Input file',message=inputfile)
        suppl = Chem.rdmolfiles.SmilesMolSupplier(
                  inputfile, titleLine=False, delimiter = " \t"
                  )
        foutput = fd.asksaveasfile(title='Choose the output file ')
        showinfo(title='Output file',message=foutput.name)
        w = Chem.SDWriter(foutput.name)
        logoutput = open(foutput.name+"_log.txt", "w")
        logoutput.close()
        mol_id=0
        for mol in suppl:
            mol_id += 1
            self.curr_mol_text.set(mol_id)
            self.window.update_idletasks()
            if mol is None:
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(str(mol_id)+": Could not be read.\n")
                logoutput.close()
                continue
            if self.Metals_flag.get()==1:
                mol = md.Disconnect(mol)
                mol = lfc.choose(mol)
            if self.Fragment_flag.get()==1:
                mol = lfc.choose(mol)
            if self.Uncharge_flag.get()==1:
                mol = u.uncharge(mol)
            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol,useRandomCoords=True,randomSeed=1)
                if self.Optimize_flag.get()==1:
                    AllChem.MMFFOptimizeMolecule(
                                                 mol,
                                                 mmffVariant='MMFF94',
                                                 maxIters=5000
                                                 )
                w.write(mol)
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(str(mol_id)+": ok\n")
                logoutput.close()
            except:
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(str(mol_id)+": No conformer was generated.\n")
                logoutput.close()
        w.close()

    def calcdesc(self):
        """Read molecules in the MDL SDFile format from a file, calculate
           descriptors according to the status of flag variables (specified in
           the GUI), save the descriptors to a csv file, register messages in
           a log file.
        """
        inputfile = fd.askopenfilename(
            title='Choose the input SDF file with the 3D structure'
        )
        showinfo(title='Input file',message=inputfile)
        suppl = Chem.SDMolSupplier(inputfile, removeHs = False)
        foutput = fd.asksaveasfile(title='Choose the output file ')
        showinfo(title='Output file',message=foutput.name)
        logoutput = open(foutput.name+"_log.txt", "w")
        logoutput.close()
        mol_id=0
        for mol in suppl:
            mol_id += 1
            self.curr_mol_text.set(mol_id)
            self.window.update_idletasks()
            if mol is None:
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(str(mol_id)+": Could not be read.\n")
                logoutput.close()
                if mol_id == 1: break
                continue
            #Calculate RDKit molecular descriptors
            try:
                if self.addHs_flag.get()==1:
                    mol = Chem.AddHs(mol,addCoords=True)
                if self.RDF_flag.get()==1:
                    descriptorsRDF=rdMolDescriptors.CalcRDF(mol)
                if self.MORSE_flag.get()==1:
                    descriptorsMORSE=rdMolDescriptors.CalcMORSE(mol)
                if self.MW_flag.get()==1:
                    descriptorsMW=rdMolDescriptors.CalcExactMolWt(mol)
                if self.WHIM_flag.get()==1:
                    descriptorsWHIM=rdMolDescriptors.CalcWHIM(mol)
                if self.AUTOCORR3D_flag.get()==1:
                    descriptorsAUTOCORR3D=rdMolDescriptors.CalcAUTOCORR3D(mol)
                if self.GETAWAY_flag.get()==1:
                    descriptorsGETAWAY=rdMolDescriptors.CalcGETAWAY(
                        mol,
                        precision=0.001
                        )
                if self.PEOE_flag.get()==1:
                    descriptorsPEOE=rdMolDescriptors.PEOE_VSA_(mol)
                if self.SMR_flag.get()==1:
                    descriptorsSMR=rdMolDescriptors.SMR_VSA_(mol)

                #Print descriptors labels in the first line
                if mol_id==1:
                    if self.RDF_flag.get()==1:
                        for i in range(len(descriptorsRDF)):
                            foutput.write(f"RDF{i+1},")
                    if self.MORSE_flag.get()==1:
                        for i in range(len(descriptorsMORSE)):
                            foutput.write(f"MORSE{i+1},")
                    if self.WHIM_flag.get()==1:
                        for i in range(len(descriptorsWHIM)):
                            foutput.write(f"WHIM{i+1},")
                    if self.AUTOCORR3D_flag.get()==1:
                        for i in range(len(descriptorsAUTOCORR3D)):
                            foutput.write(f"AUTOCORR3D{i+1},")
                    if self.GETAWAY_flag.get()==1:
                        for i in range(len(descriptorsGETAWAY)):
                            foutput.write(f"GETAWAY{i+1},")
                    if self.PEOE_flag.get()==1:
                        for i in range(len(descriptorsPEOE)):
                            foutput.write(f"PEOE{i+1},")
                    if self.SMR_flag.get()==1:
                        for i in range(len(descriptorsSMR)):
                            foutput.write(f"SMR{i+1},")
                    if self.MW_flag.get()==1:
                        foutput.write("MW,")
                    if self.VGrid_flag.get()==1:
                        centerx=float(self.centerx_str.get())
                        centery=float(self.centery_str.get())
                        centerz=float(self.centerz_str.get())
                        cube_size=float(self.cubesize_str.get())
                        atmreach=float(self.atmreach_str.get())
                        csizex=float(self.csizex_str.get())/2
                        csizey=float(self.csizey_str.get())/2
                        csizez=float(self.csizez_str.get())/2
                        grid_sizex, grid_sizey, grid_sizez =(
                            grid_coords_from_xyz(
                                csizex,csizey,csizez,0,0,0,cube_size
                            )
                        )
                        for colz in range(-grid_sizez,grid_sizez+1):
                            for coly in range(-grid_sizey,grid_sizey+1):
                                for colx in range(-grid_sizex,grid_sizex+1):
                                    foutput.write(
                                     "P{:g}_{:g}_{:g},".format(
                                                colx*cube_size+centerx,
                                                coly*cube_size+centery,
                                                colz*cube_size+centerz
                                                )
                                    )
                    if self.VoxGrid_flag.get()==1:
                        voxcenterx=float(self.voxcenterx_str.get())
                        voxcentery=float(self.voxcentery_str.get())
                        voxcenterz=float(self.voxcenterz_str.get())
                        voxcube_size=float(self.voxcubesize_str.get())
                        voxsigma=float(self.voxsigma_str.get())
                        voxcsizex=float(self.voxcsizex_str.get())/2
                        voxcsizey=float(self.voxcsizey_str.get())/2
                        voxcsizez=float(self.voxcsizez_str.get())/2
                        voxgrid_sizex, voxgrid_sizey, voxgrid_sizez =(
                            grid_coords_from_xyz(
                                                 voxcsizex,
                                                 voxcsizey,
                                                 voxcsizez,
                                                 0,
                                                 0,
                                                 0,
                                                 voxcube_size
                                                 )
                        )
                        for colz in range(-voxgrid_sizez,voxgrid_sizez+1):
                            for coly in range(-voxgrid_sizey,voxgrid_sizey+1):
                                for colx in range(
                                                  -voxgrid_sizex,
                                                  voxgrid_sizex+1
                                                  ):
                                    foutput.write(
                                        "V{:g}_{:g}_{:g},".format(
                                                colx*voxcube_size+voxcenterx,
                                                coly*voxcube_size+voxcentery,
                                                colz*voxcube_size+voxcenterz
                                                )
                                    )
                    foutput.write("NAME,ID\n")

                #Print RDKit descriptors values
                if self.RDF_flag.get()==1:
                    for desc in descriptorsRDF:
                        foutput.write("%g," %float(desc))
                if self.MORSE_flag.get()==1:
                    for desc in descriptorsMORSE:
                        foutput.write("%g," %float(desc))
                if self.WHIM_flag.get()==1:
                    for desc in descriptorsWHIM:
                        foutput.write("%g," %float(desc))
                if self.AUTOCORR3D_flag.get()==1:
                    for desc in descriptorsAUTOCORR3D:
                        foutput.write("%g," %float(desc))
                if self.GETAWAY_flag.get()==1:
                    for desc in descriptorsGETAWAY:
                        foutput.write("%g," %float(desc))
                if self.PEOE_flag.get()==1:
                    for desc in descriptorsPEOE:
                        foutput.write("%g," %float(desc))
                if self.SMR_flag.get()==1:
                    for desc in descriptorsSMR:
                        foutput.write("%g," %float(desc))
                if self.MW_flag.get()==1:
                    foutput.write(f"{descriptorsMW},")

                #Calculate grid descriptors and print them
                if self.VGrid_flag.get()==1:
                    gop = Grid_of_potential(
                            mol,
                            self.ChargeType.get(),
                            centerx,
                            centery,
                            centerz,
                            cube_size,
                            atmreach,
                            grid_sizex,
                            grid_sizey,
                            grid_sizez
                            )
                    for colz in range(-grid_sizez,grid_sizez+1):
                        for coly in range(-grid_sizey,grid_sizey+1):
                            for colx in range(-grid_sizex,grid_sizex+1):
                                foutput.write(
                                    "{:g},"
                                    .format(gop.gridV[colx][coly][colz])
                                )
                if self.VoxGrid_flag.get()==1:
                    gov = Grid_of_voxels(
                            mol,
                            self.VoxProp.get(),
                            voxcenterx,
                            voxcentery,
                            voxcenterz,
                            voxcube_size,
                            voxsigma,
                            voxgrid_sizex,
                            voxgrid_sizey,
                            voxgrid_sizez
                            )
                    for colz in range(-voxgrid_sizez,voxgrid_sizez+1):
                        for coly in range(-voxgrid_sizey,voxgrid_sizey+1):
                            for colx in range(-voxgrid_sizex,voxgrid_sizex+1):
                                foutput.write(
                                    "{:g},"
                                    .format(gov.gridVox[colx][coly][colz])
                                )
                mol = Chem.RemoveHs(mol)
                if mol.GetProp("_Name") != "":
                    foutput.write(mol.GetProp("_Name"))
                else:
                    foutput.write(Chem.MolToSmiles(mol))
                foutput.write(f",{mol_id}\n")
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(f"{mol_id}: ok\n")
                logoutput.close()
            except:
                logoutput = open(foutput.name+"_log.txt", "a")
                logoutput.write(f"{mol_id}: Descriptors not calculated.\n")
                logoutput.close()
                if mol_id == 1: break
        foutput.close()



class Grid_of_potential:
    """Represents a 3D grid of electrostatic potential (EP) generated by a
    molecule.

    Parameters:
        mol - the molecule (RDKit object)
        prop - atomic property
        centerx - x coordinate of the center of the grid
        centery - y coordinate of the center of the grid
        centerz - z coordinate of the center of the grid
        (the center of the grid is the center of the cube at the center)
        cube_size - size (in Angstrom) of the cubes making the grid
        atmreach - maximum distance of an atom to the cubes to be
                   influenced by that atom
        grid_sizex - size of the grid (in Angstrom) in the x direction
        grid_sizey - size of the grid (in Angstrom) in the y direction
        grid_sizez - size of the grid (in Angstrom) in the y direction

    Variable:
        gridV: list to store the EP at each cube

    Methods:
        __init__
        calc_grid
    """

    def __init__(
                 self,
                 mol,
                 prop,
                 centerx,
                 centery,
                 centerz,
                 cube_size,
                 atmreach,
                 grid_sizex,
                 grid_sizey,
                 grid_sizez
                 ):
        """Initialize variables,
           calculate the number of cubes (subgrid_size) in each side
            of the subgrid surrounding an atom where EP values are calculated,
           calculate the atomic charges (RDKit),
           calculate the EP at the cubes of the grid.
        """
        self.mol=mol
        self.prop=prop
        self.centerx= centerx
        self.centery= centery
        self.centerz= centerz
        self.cube_size= cube_size
        self.atmreach= atmreach
        self.grid_sizex= grid_sizex
        self.grid_sizey= grid_sizey
        self.grid_sizez= grid_sizez
        self.gridV = [[ [0 for colz in range(-grid_sizez,grid_sizez+1)]
                           for coly in range(-grid_sizey,grid_sizey+1)]
                           for colx in range(-grid_sizex,grid_sizex+1)]
        self.subgrid_size=round(atmreach/cube_size)
        if prop=="Gasteiger":
            AllChem.ComputeGasteigerCharges(mol)
        if prop=="MMFF":
            self.mmffprops=AllChem.MMFFGetMoleculeProperties(mol)
        self.calc_grid()

    def calc_grid(self):
        """For each atom in the molecule:
           - get the atomic partial charge,
           - get the atomic coordinates (atm_coords),
           - get the atomic coordinates in the grid scale
             (atm_gx, atm_gy, atm_gz),
           - calculate the EP generated by the atom at each cube of the
             subgrid around it,
           - add the EP to the value in the cube.
           (in the grid scale the coordinates of the center cube are [0][0][0])
        """
        for atm in self.mol.GetAtoms():
            if self.prop=="Gasteiger":
                atm_charge=float(atm.GetProp('_GasteigerCharge'))
            if self.prop=="MMFF":
                atm_charge=self.mmffprops.GetMMFFPartialCharge(atm.GetIdx())
            atm_coords=self.mol.GetConformer().GetAtomPosition(atm.GetIdx())
            atm_gx, atm_gy, atm_gz = grid_coords_from_xyz(
                                                          atm_coords.x,
                                                          atm_coords.y,
                                                          atm_coords.z,
                                                          self.centerx,
                                                          self.centery,
                                                          self.centerz,
                                                          self.cube_size
                                                          )
            for runx in range(-self.subgrid_size,self.subgrid_size+1):
                if (atm_gx+runx<-self.grid_sizex
                    or
                    atm_gx+runx>self.grid_sizex
                   ): continue
                for runy in range(
                                -self.subgrid_size,
                                self.subgrid_size+1
                                ):
                    if (
                        atm_gy+runy<-self.grid_sizey
                        or
                        atm_gy+runy>self.grid_sizey
                        ): continue
                    for runz in range(
                                    -self.subgrid_size,
                                    self.subgrid_size+1
                                    ):
                        if (
                            atm_gz+runz<-self.grid_sizez
                            or
                            atm_gz+runz>self.grid_sizez
                            ): continue
                        pointgrid=(
                        (atm_gx+runx)*self.cube_size+self.centerx,
                        (atm_gy+runy)*self.cube_size+self.centery,
                        (atm_gz+runz)*self.cube_size+self.centerz
                        )
                        if math.dist(
                            (atm_coords.x,atm_coords.y,atm_coords.z),
                            pointgrid
                            )>self.atmreach: continue
                        if math.dist(
                            (atm_coords.x,atm_coords.y,atm_coords.z),
                            pointgrid
                            )<1:
                            self.gridV[atm_gx+runx][
                                    atm_gy+runy][
                                    atm_gz+runz]+= atm_charge
                        else:
                            self.gridV[atm_gx+runx][
                                    atm_gy+runy][
                                    atm_gz+runz]+= (
                                            atm_charge
                                            /math.dist(
                                                (
                                                atm_coords.x,
                                                atm_coords.y,
                                                atm_coords.z
                                                ),
                                                pointgrid
                                            )
                                    )


class Grid_of_voxels:
    """Represents a 3D grid of voxels obtained from the 3D structure of a
    molecule.

    Parameters:
        mol - the molecule (RDKit object)
        prop - atomic property
        voxcenterx - x coordinate of the center of the grid
        voxcentery - y coordinate of the center of the grid
        voxcenterz - z coordinate of the center of the grid
        (the center of the grid is the center of the cube in its center)
        voxcube_size - size (in Angstrom) of the cubes making the grid
        voxsigma - the sigma parameter in the voxelisation Gaussian kernel
        voxgrid_sizex - size of the grid (in Angstrom) in the x direction
        voxgrid_sizey - size of the grid (in Angstrom) in the y direction
        voxgrid_sizez - size of the grid (in Angstrom) in the y direction

    Variable:
        gridVox: list to store the values at each voxel

    Methods:
        __init__
        calc_grid
    """

    def __init__(
                 self,
                 mol,
                 prop,
                 voxcenterx,
                 voxcentery,
                 voxcenterz,
                 voxcube_size,
                 voxsigma,
                 voxgrid_sizex,
                 voxgrid_sizey,
                 voxgrid_sizez
                 ):
        """Initialize variables,
           calculate the number of cubes in each direction (voxsubgrid_size)
             around an atom which are influenced by it,
           calculate the atomic properties (RDKit),
           calculate the atomic contributions to the voxels of the grid.
        """
        self.mol= mol
        self.prop= prop
        self.voxcenterx= voxcenterx
        self.voxcentery= voxcentery
        self.voxcenterz= voxcenterz
        self.voxcube_size= voxcube_size
        self.voxsigma= voxsigma
        self.voxgrid_sizex= voxgrid_sizex
        self.voxgrid_sizey= voxgrid_sizey
        self.voxgrid_sizez= voxgrid_sizez
        self.gridVox= [[ [0 for colz in range(-voxgrid_sizez,voxgrid_sizez+1)]
                            for coly in range(-voxgrid_sizey,voxgrid_sizey+1)]
                            for colx in range(-voxgrid_sizex,voxgrid_sizex+1)]
        self.voxsubgrid_size=round(3/voxcube_size)
        if prop=="Gasteiger":
            AllChem.ComputeGasteigerCharges(mol)
        if prop=="MMFF":
            self.mmffprops=AllChem.MMFFGetMoleculeProperties(mol)
        if prop=="MR" or prop=="LogP":
            self.atom_contributes = rdMolDescriptors._CalcCrippenContribs(mol)
        self.calc_grid()

    def calc_grid(self):
        """For each atom in the molecule:
           - get the atomic property,
           - get the atomic coordinates (voxatm_coords),
           - get the atomic coordinates in the grid scale
             (voxatm_gx, voxatm_gy, voxatm_gz),
           - calculate the contribution of the atom to each voxel (cube)
             around it (up to 3 Angstrom),
           - add the contribution to the value in the cube.
           (in the grid scale the coordinates of the center cube are
           [0][0][0])
        """
        for atm in self.mol.GetAtoms():
            if self.prop=="Gasteiger":
                voxatm_prop=float(atm.GetProp('_GasteigerCharge'))
            if self.prop=="MMFF":
                voxatm_prop=self.mmffprops.GetMMFFPartialCharge(atm.GetIdx())
            if self.prop=="AtomicNumber":
                voxatm_prop=atm.GetAtomicNum()
            if self.prop=="MR":
                voxatm_prop=float(self.atom_contributes[atm.GetIdx()][1])
            if self.prop=="LogP":
                voxatm_prop=float(self.atom_contributes[atm.GetIdx()][0])
            voxatm_coords=self.mol.GetConformer().GetAtomPosition(atm.GetIdx())
            voxatm_gx, voxatm_gy, voxatm_gz = grid_coords_from_xyz(
                                                voxatm_coords.x,
                                                voxatm_coords.y,
                                                voxatm_coords.z,
                                                self.voxcenterx,
                                                self.voxcentery,
                                                self.voxcenterz,
                                                self.voxcube_size
                                                )
            for runx in range(-self.voxsubgrid_size,self.voxsubgrid_size+1):
                if (
                   voxatm_gx+runx<-self.voxgrid_sizex
                   or
                   voxatm_gx+runx>self.voxgrid_sizex
                   ): continue
                for runy in range(
                                 -self.voxsubgrid_size,
                                 self.voxsubgrid_size+1
                                 ):
                    if (
                       voxatm_gy+runy<-self.voxgrid_sizey
                       or
                       voxatm_gy+runy>self.voxgrid_sizey
                       ): continue
                    for runz in range(
                                     -self.voxsubgrid_size,
                                     self.voxsubgrid_size+1
                                     ):
                        if (
                           voxatm_gz+runz<-self.voxgrid_sizez
                           or
                           voxatm_gz+runz>self.voxgrid_sizez
                           ): continue
                        voxpointgrid=(
                          (voxatm_gx+runx)
                            *self.voxcube_size+self.voxcenterx,
                          (voxatm_gy+runy)
                            *self.voxcube_size+self.voxcentery,
                          (voxatm_gz+runz)
                            *self.voxcube_size+self.voxcenterz
                        )
                        if math.dist(
                            (voxatm_coords.x,voxatm_coords.y,voxatm_coords.z),
                            voxpointgrid
                            )>3: continue
                        self.gridVox[voxatm_gx+runx][
                                     voxatm_gy+runy][
                                     voxatm_gz+runz] += (
                                        voxatm_prop*
                                        math.exp(
                                            -math.pow(
                                                    math.dist(
                                                        (voxatm_coords.x,
                                                            voxatm_coords.y,
                                                            voxatm_coords.z)
                                                        ,voxpointgrid
                                                        )
                                                    ,2
                                                    )
                                            /2/math.pow(self.voxsigma,2)
                                            )
                                        )

def grid_coords_from_xyz(
                  cx,
                  cy,
                  cz,
                  centerx,
                  centery,
                  centerz,
                  cube_size
                  ):
    """Convert the Cartesian coordinates to the grid scale,
       return the coordinates in the grid scale

       Parameters:
            cx - x coordinate of a point in the 3D space
            cy - y coordinate of a point in the 3D space
            cz - z coordinate of a point in the 3D space
            centerx - x coordinate of the grid center in the 3D space
            centery - y coordinate of the grid center in the 3D space
            centerz - z coordinate of the grid center in the 3D space
            cube_size - size (in Angstrom) of the cubes making the grid
       Variables:
            gx, gy, gz - coordinates of the point in the grid
       (in the grid scale the coordinates of the center cube are [0][0][0])
    """
    gx = round((cx-centerx)/cube_size)
    gy = round((cy-centery)/cube_size)
    gz = round((cz-centerz)/cube_size)
    return gx, gy, gz

def operate_cli():
    """Command-line interface for EP and voxel grids:
       Read the arguments, print the labels,
       calculate the grid according to the arguments,
       print the grid values
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='Input file')
    parser.add_argument('-o', help='Output file (optional)')
    parser.add_argument(
        '-g',
        choices = ['pot','vox'],
        required=True,
        help='Type of grid'
        )
    parser.add_argument(
        '-p',
        choices = ['Gasteiger','MMFF','MR','LogP','AtomicNumber'],
        required=True,
        help='Atomic property'
        )
    parser.add_argument(
        '-cx',
        type=float,
        required=True,
        help='x coordinate of grid center'
        )
    parser.add_argument(
        '-cy',
        type=float,
        required=True,
        help='y coordinate of grid center'
        )
    parser.add_argument(
        '-cz',
        type=float,
        required=True,
        help='z coordinate of grid center'
        )
    parser.add_argument(
        '-r',
        type=float,
        required=True,
        help='Cube size (Angstrom)'
        )
    parser.add_argument(
        '-w',
        type=float,
        required=True,
        help='Atom reach (cutoff distance or sigma)'
        )
    parser.add_argument(
        '-sx',
        type=float,
        required=True,
        help='Size of grid (xx axys)'
        )
    parser.add_argument(
        '-sy',
        type=float,
        required=True,
        help='Size of grid (yy axys)'
        )
    parser.add_argument(
        '-sz',
        type=float,
        required=True,
        help='Size of grid (zz axys)'
        )
    args = parser.parse_args()
    if args.g=='pot' and args.p not in ['Gasteiger','MMFF']:
        print ("Invalid atomic property for the electrostatic potential.")
        return
    suppl = Chem.SDMolSupplier(args.i, removeHs = False)
    if args.o != None:
        logoutput = open(args.o+"_log.txt", "w")
        logoutput.close()
    mol_id=0
    for mol in suppl:
        mol_id += 1
        if mol is None:
            if args.o is None:
                print (str(mol_id)+": Could not be read.\n")
            else:
                logoutput = open(args.o+"_log.txt", "a")
                logoutput.write(str(mol_id)+": Could not be read.\n")
                logoutput.close()
            if mol_id == 1: break
            continue
        if mol_id==1:  # Print descriptor labels in first row
            if args.o != None:
                output_file = open(args.o, "w")
            if args.g=='pot':
                grid_sizex, grid_sizey, grid_sizez = grid_coords_from_xyz(
                                                                args.sx/2,
                                                                args.sy/2,
                                                                args.sz/2,
                                                                0,
                                                                0,
                                                                0,
                                                                args.r
                                                                )
                for colz in range(-grid_sizez,grid_sizez+1):
                    for coly in range(-grid_sizey,grid_sizey+1):
                         for colx in range(-grid_sizex,grid_sizex+1):
                            if args.o is None:
                                sys.stdout.write(
                                 "P{:g}_{:g}_{:g},".format(colx*args.r+args.cx,
                                                           coly*args.r+args.cy,
                                                           colz*args.r+args.cz)
                                 )
                            else:
                                output_file.write(
                                 "P{:g}_{:g}_{:g},".format(colx*args.r+args.cx,
                                                          coly*args.r+args.cy,
                                                          colz*args.r+args.cz)
                                 )
            if args.g=='vox':
                voxgrid_sizex,voxgrid_sizey,voxgrid_sizez=grid_coords_from_xyz(
                                                                args.sx/2,
                                                                args.sy/2,
                                                                args.sz/2,
                                                                0,
                                                                0,
                                                                0,
                                                                args.r
                                                                )
                for colz in range(-voxgrid_sizez,voxgrid_sizez+1):
                    for coly in range(-voxgrid_sizey,voxgrid_sizey+1):
                        for colx in range(-voxgrid_sizex,voxgrid_sizex+1):
                            if args.o is None:
                                sys.stdout.write(
                                 "V{:g}_{:g}_{:g},".format(colx*args.r+args.cx,
                                                           coly*args.r+args.cy,
                                                           colz*args.r+args.cz)
                                 )
                            else:
                                output_file.write(
                                 "V{:g}_{:g}_{:g},".format(colx*args.r+args.cx,
                                                           coly*args.r+args.cy,
                                                           colz*args.r+args.cz)
                                 )
            if args.o != None:
                output_file.write("NAME,ID\n")
            else:
                sys.stdout.write("NAME,ID\n")
        if args.g=='pot': #Calculate and print descriptors
            gop = Grid_of_potential(
                    mol,
                    args.p,
                    args.cx,
                    args.cy,
                    args.cz,
                    args.r,
                    args.w,
                    grid_sizex,
                    grid_sizey,
                    grid_sizez)
            for colz in range(-grid_sizez,grid_sizez+1):
                for coly in range(-grid_sizey,grid_sizey+1):
                    for colx in range(-grid_sizex,grid_sizex+1):
                        if args.o is None:
                            sys.stdout.write(
                                "{:g},".format(gop.gridV[colx][coly][colz])
                                )
                        else:
                            output_file.write(
                                "{:g},".format(gop.gridV[colx][coly][colz])
                                )
        if args.g=='vox': #Calculate and print descriptors
            gov = Grid_of_voxels(
                    mol,
                    args.p,
                    args.cx,
                    args.cy,
                    args.cz,
                    args.r,
                    args.w,
                    voxgrid_sizex,
                    voxgrid_sizey,
                    voxgrid_sizez)
            for colz in range(-voxgrid_sizez,voxgrid_sizez+1):
                for coly in range(-voxgrid_sizey,voxgrid_sizey+1):
                    for colx in range(-voxgrid_sizex,voxgrid_sizex+1):
                        if args.o is None:
                            sys.stdout.write(
                             "{:g},".format(gov.gridVox[colx][coly][colz])
                             )
                        else:
                            output_file.write(
                             "{:g},".format(gov.gridVox[colx][coly][colz])
                             )
        mol = Chem.RemoveHs(mol)
        if mol.GetProp("_Name") != "":
            if args.o is None:
                sys.stdout.write(mol.GetProp("_Name"))
            else:
                output_file.write(mol.GetProp("_Name"))
        else:
            if args.o is None:
                sys.stdout.write(Chem.MolToSmiles(mol))
                sys.stdout.write(f",{mol_id}\n")
            else:
                output_file.write(Chem.MolToSmiles(mol))
                output_file.write(f",{mol_id}\n")
        if args.o != None:
            logoutput = open(args.o+"_log.txt", "a")
            logoutput.write(f"{mol_id}: ok\n")
            logoutput.close()
    if args.o != None:
        output_file.close()

if __name__ == "__main__":
    if len(sys.argv)==1:
        towindow = tk.Tk()
        TheGUI = GuidemolGUI(towindow)
        towindow.mainloop()
    else: operate_cli()
