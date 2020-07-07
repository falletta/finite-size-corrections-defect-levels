# ===========================================================================
# FINITE-SIZE CORRECTION OF DEFECT ENERGY LEVELS INVOLVING IONIC POLARIZATION
# ===========================================================================

# Libraries
import sys
import os
import numpy as np
import scipy as sp
import scipy.integrate as integrate
import matplotlib
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({'font.size': 24,'legend.fontsize': 20,'legend.handlelength': 0.5})

#Constants
hartree2eV = 27.211396132
bohr2A = 0.529177249
A2bohr = 1.8897259886

# Gaussian function
def gaussian(x, sig, xmax):
    return 1/((2.0*np.pi*sig**2)**0.5) * np.exp( -0.5 * (min(x, xmax-x)/sig)**2 )
    
# Gaussian in K-space for unitary charge
def rho(k,σ):
	return np.exp(-0.5 * σ**2 * k**2)	

# Convolution of two arrays with same length
def convolution(f,g):
    dim=len(f)
    res=np.zeros(dim)																																								
    for i in range(dim):
        s=0
        for j in range(dim):
            s+=f[j]*g[i-j]
        res[i]=s				
    return res

# Trasform fractional coordinates to cartesian coordinates
def from_fractional_to_cartesian(v, p):
    [a,b,c] = p[0:3]
    [alpha,beta,gamma] = [x*np.pi/180 for x in p[3:] ]
    omega = a*b*c*np.sqrt(1-np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    Mfc=[[a, b*np.cos(gamma), c*np.cos(beta)],
         [0.0, b*np.sin(gamma), c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/(np.sin(gamma))], 
         [0.0, 0.0, omega/(a*b*np.sin(gamma))]]
    return np.dot(Mfc, v)

# Parse input file
def parser(infile):
    with open(infile, 'r') as f:
        data = f.readlines()
    for line in data:
        s = line.split()
        if len(s) > 1:
            if s[0] == "system":
                system = s[1]
            if s[0] == "latt_param[A,deg]":
                p = [float(x) for x in s[1:]]
            if s[0] == "r_defect[A]":
                rdef = [float(x) for x in s[1:]]
            if s[0] == "E_cutoff[Ry]":
                Ecut = float(s[1])
            if s[0] == "direction":
                dir = s[1]
            if s[0] == "sigma[Bohr]":
                σ = float(s[1])
            if s[0] == "eps_0":
                ε0 = float(s[1])
            if s[0] == "eps_inf":
                εinf = float(s[1])  
            if s[0] == "alignment":
                if s[1] in ["Yes", "yes", "1"]:
                    align = True
                else:
                    align = False
    states = []; Vtots  = {}; r_idir = 0
    for i,line in enumerate(data):
        s = line.split()
        if len(s) > 1:
            if s[0] == "N_states":
                N = int(s[1])
                for line2 in data[i+2:i+2+N]:
                    s2 = line2.split()
                    qC = float(s2[0])
                    qR = float(s2[1])
                    states.append((qC,qR))
                    if align:
                        Vtots[(qC,qR)] = Potential(s2[2],dir)
                        if qC == 0 and qR == 0:
                            r_idir = Vtots[(0,0)].r_idir
    if align:
        if (0,0) not in states or (0,0) not in Vtots.keys():
            print("Error: the potential of the state (0,0) is missing")
            sys.exit()
    parameters = [system, p, rdef, r_idir, Ecut, dir, σ, ε0, εinf, align]
    return [parameters, states, Vtots]
    

# Electrostatic potential
class Potential:
    
    def __init__ (self, cubefile, dir):
        
        # Parse cubefile
        data = []
        with open(cubefile, 'r') as f:
            lines = f.readlines()
        N = int(lines[2].split()[0])
        nx, ny, nz = [int(line.split()[0]) for line in lines[3:6]]
        a = np.sqrt(sum([float(x)**2 for x in lines[3].split()[1:]]))*bohr2A
        b = np.sqrt(sum([float(x)**2 for x in lines[4].split()[1:]]))*bohr2A
        c = np.sqrt(sum([float(x)**2 for x in lines[5].split()[1:]]))*bohr2A
        for line in lines[(6+N):]:
            data.extend([float(x) for x in line.split()])
            
        # xyz grid
        x = np.arange(nx)*a
        y = np.arange(ny)*b
        z = np.arange(nz)*c
        
        # Potential
        V_3D = sp.zeros((nx, ny, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    V_3D[i,j,k] = data[i*ny*nz + j*nz + k] * hartree2eV
        
        # Projected potential
        if dir == "x":
            i1 = 1; i2 = 1; normalization = ny*nz; self.r_idir = x
        if dir == "y":
            i1 = 0; i2 = 1; normalization = nx*nz; self.r_idir = y
        if dir == "z":
            i1 = 0; i2 = 0;	normalization = nx*ny; self.r_idir = z
        V = sp.sum(sp.sum(V_3D, axis=i1), axis=i2) / normalization
        
        # Invert sign
        self.V = -V                
        
# FNV corrections
def FNV_corr(parameters):
    
    # inputs and unit conversions
    [system, p, rdef, r_idir, Ecut, dir, σ, ε0, εinf, align] = parameters
    if dir == "x":
        idir = 0
    if dir == "y":
        idir = 1
    if dir == "z":
        idir = 2
    rdef = [x*A2bohr for x in rdef]
    p = [float(x)*A2bohr for x in p[0:3]]+[float(x) for x in p[3:6]]
    a = p[0]; b = p[1]; c = p[2]
    Gmax = np.sqrt(2*Ecut)
    
    # kgrid (lengths in Bohr)
    ax = from_fractional_to_cartesian([1,0,0], p)
    ay = from_fractional_to_cartesian([0,1,0], p)
    az = from_fractional_to_cartesian([0,0,1], p)
    vol = np.dot(ax, np.cross(ay, az))
    bx = 2*np.pi*np.cross(ay, az)/vol
    by = 2*np.pi*np.cross(az, ax)/vol
    bz = 2*np.pi*np.cross(ax, ay)/vol
    volg = np.dot(bx, np.cross(by, bz))
    nx = 2*int(round(Gmax/(2*np.pi/a)))						
    ny = 2*int(round(Gmax/(2*np.pi/b)))						
    nz = 2*int(round(Gmax/(2*np.pi/c)))						
    gx = np.fft.fftfreq(nx, 1/(nx*np.sqrt(np.dot(bx, bx))))
    gy = np.fft.fftfreq(ny, 1/(ny*np.sqrt(np.dot(by, by))))
    gz = np.fft.fftfreq(nz, 1/(nz*np.sqrt(np.dot(bz, bz))))
    n = [nx, ny, nz]
    g = [gx, gy, gz]
    
    # Isolated defect energy (in eV), for q=1 and ε=1
    Eiso = 1.0/(np.pi)*integrate.quad(lambda k: (rho(k,σ))**2, 0, Gmax)[0]*hartree2eV	
    
    # Periodic defect energy (in eV), for q=1 and ε=1
    Eper = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not (i == 0 and j == 0 and k == 0):
                    G = np.sqrt(gx[i]**2 + gy[j]**2 + gz[k]**2)
                    if G <= Gmax:
                        Eper += (2*np.pi)/(vol) * (rho(G,σ))**2 / (G**2)
    Eper *= hartree2eV
    
    # Model potential for a Gaussian charge distribution (in eV)
    Vgauss_r = []
    if align:
        Vgauss_g = np.zeros(n[idir],dtype=complex)
        for i in range(1,n[idir],1):
            Vgauss_g[i] = (4*np.pi)*rho(g[idir][i],σ)/(g[idir][i]**2)*np.exp(-1j*g[idir][i]*rdef[idir])
        Vgauss_g *= hartree2eV
        Vgauss_r_tmp = np.fft.ifft(Vgauss_g).real*n[idir]/vol
        Vgauss_r_int = interp1d(np.linspace(0, p[idir]*bohr2A, n[idir]), Vgauss_r_tmp)        
        Vgauss_r = np.array(Vgauss_r_int(r_idir))
    
    return [Eiso, Eper, Vgauss_r]
    
# Falletta-Wiktor-Pasquarello corrections
class FWP_corr:
    
    #Initialization
    def __init__ (self, state, parameters, Eiso, Eper, Vgauss, Vtots):
        
        #Parse inputs
        [self.system, self.p, self.rdef, self.r_idir, self.Ecut, 
         self.dir, self.σ, self.ε0, self.εinf, self.align] = parameters
        self.qC     = state[0]
        self.qR     = state[1]
        self.qpol   = -self.qR*(1-self.εinf/self.ε0)
        self.label  = '$({:+.2f},\mathbf{{R}}_{{{:+.2f}}})$'.format(self.qC,self.qR)
        self.Vtots  = Vtots
        self.Eiso   = Eiso
        self.Eper   = Eper
        self.Vgauss = Vgauss
        if self.dir == "x":
            self.idir = 0
        if self.dir == "y":
            self.idir = 1
        if self.dir == "z":
            self.idir = 2
        if self.align:
            self.i_ΔV = np.argmin(self.Vgauss)
                
    # Corrections for total energies
    def corrections_total_energy(self):
        E_madelung = (self.Eiso-self.Eper)*(self.qR**2/self.ε0 - (self.qR+self.qpol)**2/self.εinf + (self.qC+self.qpol)**2/self.εinf ) 
        E_alignment = 0
        if self.align:
            E_alignment = self.calculate_alignment()
        self.Ecor = E_madelung - E_alignment
        
    # Corrections for single-particle levels
    def corrections_KS_level(self):
        E_madelung_KS  = (self.Eiso-self.Eper)*(self.qC+self.qpol)**2/self.εinf
        E_alignment_KS = 0
        if self.align:
            V_KS      = self.Vtots[(self.qC,self.qR)].V - self.Vtots[(0,0)].V
            Vgauss_KS = self.Vgauss * (self.qC + self.qpol)/self.εinf
            E_alignment_KS = (self.qC+self.qpol) * self.calculate_alignment_FNV(V_KS, Vgauss_KS, "KS", plot_flag=False)
        self.εKS_cor =  - 2 * (E_madelung_KS - E_alignment_KS) / (self.qC + self.qpol)
            
    # Alignment for state (qC,qR)
    def calculate_alignment(self):
    
        # DFT defect potentials wrt the (0,0) state
        V1 = self.Vtots[(self.qR,self.qR)].V - self.Vtots[(0,0)].V
        V2 = self.Vtots[(self.qR,self.qR)].V - self.Vtots[(0,0)].V
        V3 = self.Vtots[(self.qC,self.qR)].V - self.Vtots[(0,0)].V
        
        # Model defect potentials
        Vgauss1 = self.Vgauss *  self.qR/self.ε0 
        Vgauss2 = self.Vgauss * (self.qR+self.qpol)/self.εinf
        Vgauss3 = self.Vgauss * (self.qC+self.qpol)/self.εinf
        
        # Labels
        label1 = ": $ΔV(q',ε_0)$"
        label2 = ": $ΔV(q'+ q'_\mathrm{pol},ε_{\infty})$"
        label3 = ": $ΔV(q+q'_\mathrm{pol},ε_{\infty})$"
        
        #Compute single alignments
        qΔV_1 =  self.qR            * self.calculate_alignment_FNV(V1, Vgauss1, label1)
        qΔV_2 = (self.qR+self.qpol) * self.calculate_alignment_FNV(V2, Vgauss2, label2)
        qΔV_3 = (self.qC+self.qpol) * self.calculate_alignment_FNV(V3, Vgauss3, label3)
        qΔV   = qΔV_1 - qΔV_2 + qΔV_3
        
        return qΔV
            
    # Alignment for state with qC = qR
    def calculate_alignment_FNV(self, V, Vgauss, extralabel, plot_flag=True):
        
        #Convolution of V with a gaussian of the same width as the model potential
        N = len(V)
        V_convolution = [0]*N
        for i in range(N):
            V_convolution[i] = gaussian(self.r_idir[i], self.σ, self.r_idir[-1]) / self.p[self.idir]
        V_convolution = convolution(V, V_convolution)
        area_below_gaussian = integrate.quad(lambda x: gaussian(x,self.σ,self.r_idir[-1]), 0, self.r_idir[-1])[0]
        V_convolution = [x/area_below_gaussian for x in V_convolution]
        
        #Useful quantities
        r_idir_meas = self.r_idir[self.i_ΔV]
        ΔV = V_convolution-Vgauss
        ΔV_far = ΔV[self.i_ΔV]
        ymin = min(V_convolution[self.i_ΔV],Vgauss[self.i_ΔV])
        ymax = max(V_convolution[self.i_ΔV],Vgauss[self.i_ΔV])
        length = ( max(max(max(V),max(Vgauss)), max(V_convolution)) - min(min(min(V),min(Vgauss)), min(V_convolution)) ) / 10
        
        #Plot the alignment
        if plot_flag:
            f  = plt.figure(figsize=(6,6), dpi=60)
            ax = plt.gca()
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(2.5)
            ax.xaxis.set_tick_params(which='major', width=3.0, length=12, direction="in")
            ax.yaxis.set_tick_params(which='major', width=3.0, length=12, direction="in")
            plt.gcf().subplots_adjust(left=0.15, bottom=0.19, top=0.79, right=0.95)
            plt.title(self.label+extralabel, fontsize=24)
            plt.xlabel(self.dir+" (Å)")
            plt.xlim(self.r_idir[0], self.r_idir[-1])
            plt.plot(self.r_idir, Vgauss, label='$V_\mathrm{m}$', linewidth=2.5)
            plt.plot(self.r_idir, V, label="$V_\mathrm{DFT}$", linewidth=2.5) 
            plt.plot(self.r_idir, V_convolution, ls=':', label="$\\langle V_\mathrm{DFT} \\rangle$", linewidth=2.5)
            plt.annotate('', xy=(r_idir_meas,ymax+length), xycoords='data',xytext=(r_idir_meas, ymax), textcoords='data',arrowprops=dict(zorder=100, arrowstyle='<-', color='red', linewidth=2))
            plt.annotate('', xy=(r_idir_meas,ymin-length), xycoords='data',xytext=(r_idir_meas, ymin), textcoords='data',arrowprops=dict(zorder=100, arrowstyle='<-', color='red', linewidth=2))
            plt.annotate('ΔV = {:1.2f} eV'.format(ΔV_far), xy=(max(self.r_idir)/2, min(min(V), min(Vgauss))), xycoords='data',xytext=(-10, +10), color='red', textcoords='offset points')
            plt.legend(frameon=False)
            pdf.savefig()
            plt.close()
    
        return ΔV_far

# Executing the script
print("Calculating FWP finite-size corrections")
[parameters, states, Vtots] = parser(sys.argv[1])
[Eiso, Eper, Vgauss] = FNV_corr(parameters)
with PdfPages(parameters[0]+".pdf") as pdf:
    for state in states:
        if state != (0,0):
            fs_defect = FWP_corr(state, parameters, Eiso, Eper, Vgauss, Vtots)
            fs_defect.corrections_total_energy()
            fs_defect.corrections_KS_level()
            print('\nState: ({:+.2f}, R{:+.2f})'.format(state[0],state[1]))
            print('Ecor    = {:+.6f} eV'.format(fs_defect.Ecor))
            print('eps_cor = {:+.6f} eV'.format(fs_defect.εKS_cor))
if not parameters[-1]:
    os.remove(parameters[0]+".pdf")