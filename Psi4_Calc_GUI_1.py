
# -*- coding: utf-8 -*-
from concurrent.futures import thread
from fileinput import filename
import tkinter as tk
from tkinter import BOTH, ttk
from tkinter import StringVar, messagebox
from turtle import onclick
from tkinter import filedialog
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from cProfile import label
import seaborn as sns
from tkinter.ttk import Entry
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from  rdkit import rdBase, Chem
from rdkit.Chem import AllChem
import py3Dmol as p3d

def data_import(file):
    # file select by filedialog 
    type=[('mol file','*.mol'), ('xyz file', '*.xyz')]
    file=filedialog.askopenfilename(filetypes=type)
    filename_sv.set(file)
    fname=filename_sv.get()

def Geom_Opt():
    # Geometry optimization
    # getting meth(calculation method), func(function), base(basis set) 
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Display 'Start'
    Label10=ttk.Label(text=u'Geometry optimization was started', font=("Times","12"))
    Label10.place(x=20, y=570)
    
    # geometry optimize calculation (method in HF and MP2, function in DFT)
    if meth=='HF':
        
        psi4.optimize(meth+'/'+base, molecule=molgeom)
        
    elif meth=='DFT':
        psi4.optimize(func+'/'+base, molecule=molgeom)
        
    elif meth=='MP2':
        psi4.optimize(meth+'/'+base, molecule=molgeom)
    
    optimized_geom = molgeom.save_string_xyz_file()
    
    # Display 'Calculation finish'
    Label11=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label11.place(x=20, y=590)
    # Save optimized geometry
    svname=savename_sv.get()
    with open (svname+'-optgeom.xyz','w') as f:
        f.write(optimized_geom)
                
def Vib_Calc():
    # Frequency calculation
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Display 'Start'
    Label10=ttk.Label(text=u'Calculations(Frequency) were started', font=("Times","12"))
    Label10.place(x=20, y=570)
    
    # Frequency calculation (method in HF and MP2, function in DFT)
    if meth=='HF':        
        energy, wfn = psi4.frequency(meth+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='DFT':
        energy, wfn = psi4.frequency(func+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='MP2':
        energy, wfn = psi4.frequency(meth+'/'+base, molecule=molgeom, return_wfn=True)
    
    # making frequency list
    freqdata=[]

    for i in range(30):
        freqdatum=wfn.frequencies().get(0,i)
        freqdata.append(freqdatum)

    freqlist = np.array(freqdata)
    freqposi=freqlist[(freqlist>1)&(freqlist<5000)]
    freqnega=freqlist[(freqlist<0)]
    
    FreqP=pd.DataFrame(np.array(freqposi), columns=['frequency/nm'])
    FreqN=pd.DataFrame(np.array(freqnega), columns=['frequency/nm'])
    pd.options.display.precision =2
    
    # Display 'Calculation finish'
    Label10.place_forget()
    Label11=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label11.place(x=20, y=570)

    # Sub window 1 (Display calculation results)

    sub_window1=tk.Toplevel()
    sub_window1.title('Calculations results')
    sub_window1.geometry('720x540')

    if FreqN.empty:
        LabelS_1=ttk.Label(sub_window1, text=u'No imaginary frequency', font=("Arial","12","bold"))
        LabelS_1.place(x=10, y=10, width=200)
        print  ("No imaginary frequency")   
    else :
        LabelS_2=ttk.Label(sub_window1, text=FreqN, font=("Arial","12"))
        LabelS_2.place(x=10, y=10, width=200)
        print (FreqN)

    LabelS_3=ttk.Label(sub_window1, text=FreqP, font=("Arial","12"))
    LabelS_3.place(x=30, y=50, width=200)
    print (FreqP) 
    

def MO_anal():
    # Molecular orbital analysis
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    outputname_fchk=savename_sv.get()+'.fchk'
    
    # Display 'Start'
    Label12=ttk.Label(text=u'Calculations(MO analysis) were started', font=("Times","12"))
    Label12.place(x=20, y=570)
    
    # Molecular orbital analysis (method in HF and MP2, function in DFT)
    if meth=='HF':        
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='DFT':
        energy, wfn = psi4.energy(func+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='MP2':
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
    
    # Save 'fchk' file
    psi4.fchk(wfn, outputname_fchk)
    
    # Display 'Calculation finish'
    Label12.place_forget()
    Label13=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label13.place(x=20, y=570)
    
    # Display HOMO LUMO level
    # orbital_level : number of display level
    orbital_level=5
    orblev=np.array(list(range(((orbital_level)+1))))

    HOMOdata=[]
    LUMOdata=[]
    Hindex=[]
    Lindex=[]
    LUMO_idx=wfn.nalpha()
    HOMO_idx=LUMO_idx-1

    for i in orblev:

        nexthomo=HOMO_idx-(i)
        nextlumo=LUMO_idx+(i)
        nhomolev=wfn.epsilon_a_subset("AO","ALL").np[nexthomo]
        nlumolev=wfn.epsilon_a_subset("AO","ALL").np[nextlumo]
    
        print('homo-'+str(i),round(nhomolev,5), ' a.u.','level:', nexthomo)
        print('lumo+'+str(i),round(nlumolev,5), ' a.u.', 'level:', nextlumo)
    
        HOMOdata.append(nhomolev)
        LUMOdata.append(nlumolev)
        Hindex.append('homo-'+str(i))
        Lindex.append('lumo+'+str(i))

    HOMOdf=pd.DataFrame(HOMOdata, index=Hindex, columns=['energy/a.u.'])
    LUMOdf=pd.DataFrame(LUMOdata, index=Lindex, columns=['energy/a.u.'])

    sub_window1=tk.Toplevel()
    sub_window1.title('Calculations results')
    sub_window1.geometry('720x540')
    
    LabelS_4=ttk.Label(sub_window1, text='HOMO, LUMO Energy', font=("Arial","14", 'bold'))
    LabelS_4.place(x=10, y=10, width=300)
    
    LabelS_5=ttk.Label(sub_window1, text=HOMOdf, font=("Arial","12"))
    LabelS_5.place(x=10, y=30, width=150)
    print (HOMOdf)
    
    LabelS_6=ttk.Label(sub_window1, text=LUMOdf, font=("Arial","12"))
    LabelS_6.place(x=200, y=30, width=150)
    print (LUMOdf)

def gaussBand(x, band, strength, stdev): 
    "Produces a Gaussian curve"
    bandshape =  constant* (strength/(1.0/stdev)) * np.exp(-((((1.0/x)-(1.0/band))**2)/(1.0/stdev)**2))
    return bandshape


def UV_VIS_Spec():
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Display 'start'
    Label14=ttk.Label(text=u'Calculations(UV_VIS_Spectrum) were started', font=("Times","12"))
    Label14.place(x=20, y=570) 
    
    #　TDSCF calculation (method in HF and MP2, function in DFT)

    psi4.set_options({
    'save_jk':True,
    })

    if meth=='HF':        
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
        res = tdscf_excitations (wfn, states =10)
    elif meth=='DFT':
        energy, wfn = psi4.energy(func+'/'+base, molecule=molgeom, return_wfn=True)
        res = tdscf_excitations (wfn, states =10)
    elif meth=='MP2':
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
        res = tdscf_excitations (wfn, states =10)
    
    # Display 'finish'
    Label14.place_forget()
    Label15=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label15.place(x=20, y=570)
    
    # Calculation of spectrum data from energy and osstrength (ref. Hill Research Group: http://www.grant-hill.group.shef.ac.uk/plot-uv.html)
    
    exenergy = [r["EXCITATION ENERGY"] for r in res]
    osstrength = [r["OSCILLATOR STRENGTH (LEN)"] for r in res]

    arrayex = np.array(exenergy)
    arrayos = np.array(osstrength)
    arraynm0 =[(1240/(arrayex*27.21162))]
    arraynm = np.array(arraynm0)

    arrayexr = np.reshape(arrayex,(-1,1))
    arrayosr = np.reshape(arrayos,(-1,1))
    arraynmr = np.reshape(arraynm,(-1,1))
    
    listexr=np.round(arrayexr, decimals=5).tolist()
    listosr=np.round(arrayosr, decimals=5).tolist()
    listnmr=np.round(arraynmr, decimals=2).tolist()
    
    ExEnergydf=pd.DataFrame({'Excitation Energy /au     ': listexr, 'Excitation Energy /nm     ':listnmr, 'Oscilleator Strength     ':listosr})  

    # Display energy data and spectrum
    sub_window1=tk.Toplevel()
    sub_window1.title('Calculations results (Energy Data)')
    sub_window1.geometry('720x540')
    
    LabelS_7=ttk.Label(sub_window1, text='Excitation Energy and Oscilleator Strength', font=("Arial","14", 'bold'))
    LabelS_7.place(x=10, y=10, width=300)
    
    LabelS_8=ttk.Label(sub_window1, text=ExEnergydf, font=("Arial","12"))
    LabelS_8.place(x=10, y=50)

    start=50
    finish=400
    points=300

    # En/hc =2000cm-1 = 5000 nm
    stdev = 5000
    # constant a'/hc
    global constant
    constant = 0.0003
    # Excitation energies in nm
    bands = arraynmr
    # Oscillator strengths (dimensionless)
    f = arrayosr
    # Basic check that we have the same number of bands and oscillator strengths
    if len(bands) != len(f):
        print('Number of bands does not match the number of oscillator strengths.')

    # Information on producing spectral curves (Gaussian) is adapted from Hitachi.resercher Kinka seminer
    # Gaussian curves are for fit for UV/Vis.

    x = np.linspace(start,finish,points)
    composite = 0
    for count,peak in enumerate(bands):
        thispeak = gaussBand(x, peak, f[count], stdev)

        composite += thispeak
    
    sub_window2=tk.Toplevel()
    sub_window2.title('Calculations results (Spectrum)')
    sub_window2.geometry('720x540')
    
    fig = Figure()
    ax = fig.add_subplot(1,1,1)
    
    fig_canvas = FigureCanvasTkAgg (fig, master=sub_window2)
    fig_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(fig_canvas, sub_window2)
    toolbar.update()
    fig_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    ax.plot(x,composite)
    ax.set_xlabel('$\lambda$ / nm')
    ax.set_ylabel('Intensity / unit')


def Mul_Charge():
    # Calculation of Mulliken Charges
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Display 'start'
    Label15=ttk.Label(text=u'Calculations(Mulliken analysis) were started', font=("Times","12"))
    Label15.place(x=20, y=570)
    
    # Molecular orbital analysis (method in HF and MP2, function in DFT)
    if meth=='HF':        
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='DFT':
        energy, wfn = psi4.energy(func+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='MP2':
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
    
    # Display 'finish'
    Label15.place_forget()
    Label16=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label16.place(x=20, y=570)   
    
    # Calculation of Mulliken charges
    psi4.oeprop(wfn,'MULLIKEN_CHARGES')
    mulliken = np.array(wfn.atomic_point_charges())

    atom_symbol=[]
    mulliken_data=[]

    for n in range(molgeom.natom()):
        print (molgeom.symbol(n), round(mulliken[n],4))
        atom_symbol.append(molgeom.symbol(n))
        mulliken_data.append(round(mulliken[n],4))

    Mullikendf=pd.DataFrame(({'Atom': atom_symbol, 'Mulliken Charge':mulliken_data}))

    # Display results of Mullken data
    sub_window3=tk.Toplevel()
    sub_window3.title('Calculations results (Mulliken Data)')
    sub_window3.geometry('720x540')
    
    LabelS_9=ttk.Label(sub_window3, text='Mulliken Charges', font=("Arial","14", 'bold'))
    LabelS_9.place(x=10, y=10, width=300)
    
    LabelS_10=ttk.Label(sub_window3, text=Mullikendf, font=("Arial","12"))
    LabelS_10.place(x=10, y=50)

def Dip_Moment():
    #Calculation of Dopole moment
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Display 'start'
    Label17=ttk.Label(text=u'Calculations(Dipole analysis) were started', font=("Times","12"))
    Label17.place(x=20, y=570)
    
    # Molecular orbital analysis (method in HF and MP2, function in DFT)
    if meth=='HF':        
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='DFT':
        energy, wfn = psi4.energy(func+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='MP2':
        energy, wfn = psi4.energy(meth+'/'+base, molecule=molgeom, return_wfn=True)
    
    # Display 'finish'
    Label17.place_forget()
    Label18=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label18.place(x=20, y=540)   
    
    # Calculation of dipole moment
    psi4.oeprop(wfn,'DIPOLE')
    dipole_vec = psi4.variable('scf dipole')

    # Display results of dipole monment from x, y, x vector data
    Dipmom=round(np.sqrt(np.sum(dipole_vec **2)),3)
    Label19=ttk.Label(text=f'Dipole Moment= {Dipmom} D', font=("Times","12"))
    Label19.place(x=20, y=570)   

def Polar():
    # Calculation of Polarizability 
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    outputname=savename_sv.get()+'.txt'
    
    # Display 'start'
    Label20=ttk.Label(text=u'Calculations(polarizability analysis) were started', font=("Times","12"))
    Label20.place(x=20, y=570)
    
    # Calculation of polarizability at cc2 
    
    polar = psi4.properties('cc2/'+base, properties=['polarizability'], molecule=molgeom)
    
    # Display 'finish'
    Label20.place_forget()
    Label21=ttk.Label(text=u'Calculation was finished', font=("Times","12"))
    Label21.place(x=20, y=540)   

    with open (outputname) as f:
        alltext= f.readlines()
        endgyou = len(alltext)
        endtxt = alltext[endgyou-1].strip()
    print (endtxt)

    Label21=ttk.Label(text=f'Polalizability  {endtxt} ', font=("Times","12"))
    Label21.place(x=20, y=570)   

def calc_run(self):
    
    chr=charge.get()
    mul=multi.get()
    global fname
    fname=filename_sv.get()
    
    
    #  Conformer search by rdkit then make xyz file
    with open (fname) as f:
        if '.mol' in fname:
            mol_H=f.read()
            
            confs = AllChem.EmbedMultipleConfs(mol_H, 10, pruneRmsThresh=1)
            prop = AllChem.MMFFGetMoleculeProperties(mol_H)

            energy=[]
            for conf in confs:
                mmff = AllChem.MMFFGetMoleculeForceField(mol_H, prop,confId=conf)
                mmff.Minimize()
                energy.append((mmff.CalcEnergy(), conf))
    
                conflist = np.array(energy)
                sortconf=conflist[np.argsort(conflist[:,0]) ]

                stconfid=sortconf[0,1]
                stconfid2=int(stconfid)
                stconfgeom=mol_H.GetConformer(stconfid2)

            xyz = chr, mul
            for atom, (x,y,z) in zip(mol_H.GetAtoms(), stconfgeom.GetPositions()):
                xyz += '\n'
                xyz += '{}\t{}\t{}\t{}'.format(atom.GetSymbol(), x, y, z)
                
            #  Setting input file
            global molgeom
            molgeom = psi4.geometry(xyz)

        elif '.xyz' in fname:
            xyz=f.read()
            
            print (xyz)
            #  Setting input file
            molgeom = psi4.geometry(xyz)

        
    # Set calculation method
    thre=threads.get()
    mem=memory.get()
    psi4.set_num_threads(nthread=thre)
    psi4.set_memory(mem*1000000)
    
    # Set output files
    outputname=savename_sv.get()+'.txt'
    psi4.set_output_file(outputname)
    
    # Select calculation task
    task=task_sv.get()
    if task=='Geometry optimization':
        Geom_Opt()
        return molgeom
        
    elif task=='Vibration analysis':
        Vib_Calc()
        return molgeom
    
    elif task=='Molecular orbitals analysis':
        MO_anal()
        return molgeom
        
    elif task=='UV-Vis spectrum':
        UV_VIS_Spec()
        
    elif task=='Muliken charges':
        Mul_Charge()
        
    elif task=='Dipole moment':
        Dip_Moment()
        
    elif task=='Polarizability':
        Polar()
        
        
# finish program
def scry_finish():
    exit()

#　Tkinter main 

root = tk.Tk()
root.title("Psi4 Calculation Setup")
root.geometry('800x650')

# Select molecule file
Label1=ttk.Label(root, text=u'Molecule',font=("Times","14","bold"))
Label1.place(x=20, y=60)

filename_sv = tk.StringVar()
filenameEntry = ttk.Entry(width=60, text="", textvariable=filename_sv)
filenameEntry.place(x=20, y= 90)

Button1 = ttk.Button(text=u'Select',width=10)
Button1.bind("<Button-1>", data_import) 
Button1.place(x=600, y=90)

# Select calculation task
Label2=ttk.Label(text=u'Task', font=("Times","14","bold"))
Label2.place(x=20, y=140)
task_sv = tk.StringVar()
task_contents=('Geometry optimization','Vibration analysis', 'Molecular orbitals analysis', 'UV-Vis spectrum','Muliken charges','Dipole moment', 'Polarizability')
comboBox2=ttk.Combobox(root, height=5, width=20, state='readonly', values=task_contents, textvariable=task_sv)
comboBox2.place(x=20, y=170)

# Select calculation methods
Label3=ttk.Label(text=u'Calculation Method', font=("Times","14","bold"))
Label3.place(x=20, y=220)

Label3_1=ttk.Label(root, text=u'Method',font=("Times","12"))
Label3_1.place(x=20, y=240)
method_sv = tk.StringVar()
methods=('HF','DFT', 'MP2')
comboBox3_1=ttk.Combobox(root, height=5, width=10, state='readonly', values=methods, textvariable=method_sv)
comboBox3_1.place(x=20, y=260)

Label3_2=ttk.Label(root, text=u'Function',font=("Times","12"))
Label3_2.place(x=200, y=240)
function_sv = tk.StringVar()
functions=('','b3lyp','cam-b3lyp', 'edf2','m06', 'pbe','wb97x-d')
comboBox3_2=ttk.Combobox(root, height=5, width=10, state='readonly', values=functions, textvariable=function_sv)
comboBox3_2.place(x=200, y=260)

Label3_3=ttk.Label(root, text=u'Basis set',font=("Times","12"))
Label3_3.place(x=400, y=240)
baseset_sv = tk.StringVar()
base_sets=('3-21g','6-31g', '6-31g(d)','6-311g', 'aug-cc-pvtz')
comboBox3_3=ttk.Combobox(root, height=5, width=10, state='readonly', values=base_sets, textvariable=baseset_sv)
comboBox3_3.place(x=400, y=260)

# Select options
Label4=ttk.Label(text=u'Options', font=("Times","14","bold"))
Label4.place(x=20, y=300)

Label4_1=ttk.Label(root, text=u'Tread',font=("Times","12"))
Label4_1.place(x=30, y=320)
threads=tk.IntVar(value=2)
textBox1_1=ttk.Entry(root, width=5, textvariable=threads)
textBox1_1.place(x=30, y=340)

Label4_2=ttk.Label(root, text=u'Memory/MB',font=("Times","12"))
Label4_2.place(x=130, y=320)
memory=tk.IntVar(value=500)
textBox1_2=ttk.Entry(root, width=5, textvariable=memory)
textBox1_2.place(x=130, y=340)

Label4_3=ttk.Label(root, text=u'Charge',font=("Times","12"))
Label4_3.place(x=30, y=370)
charge=tk.IntVar(value=0)
textBox2_1=ttk.Entry(root, width=5, textvariable=charge)
textBox2_1.place(x=30, y=390)

Label4_4=ttk.Label(root, text=u'Multiplicity',font=("Times","12"))
Label4_4.place(x=130, y=370)
multi=tk.IntVar(value=1)
textBox2_2=ttk.Entry(root, width=5, textvariable=multi)
textBox2_2.place(x=130, y=390)

# Input the name of calculated output files
Label5=ttk.Label(text=u'Name of output files', font=("Times","14","bold"))
Label5.place(x=20, y=450)

savename_sv = tk.StringVar()
textBox3=ttk.Entry(root, width=30, textvariable=savename_sv)
textBox3.place(x=30, y=490)

# Calculation Run
Button2=ttk.Button(text=u'Run',width=20)
Button2.bind("<Button-1>", calc_run) 
Button2.place(x=300, y=530)

Label6=ttk.Label(text=u'Finish the program')
Label6.place(x=630, y=550)
Button3 = ttk.Button(text=u'Quit',width=10, command=scry_finish)
Button3.place(x=630, y=570)


root.mainloop()