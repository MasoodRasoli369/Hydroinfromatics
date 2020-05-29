# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:40:15 2020

@author: mra013
"""
#File='..\\inputs\Extra_Test.txt'
#File='..\\inputs\Static_Test_1.txt'
#File='..\\inputs\Static_Test_2.txt'
#File='..\\inputs\Static_Test_3.txt'
#File='..\\inputs\Steady_Test.txt'
#File='..\\inputs\Transient_Test_A_1.txt'
#File='..\\inputs\Transient_Test_A_2.txt'
#File='..\\inputs\Transient_Test_A_3.txt'
#File='..\\inputs\Transient_Test_B_1.txt'
File='..\\inputs\Transient_Test_B_2.txt'
#File='..\\inputs\Transient_Test_B_3.txt'
import numpy as np
#import os
# =============================================================================
# # Importing Data Part    
# =============================================================================
global L,T,Del_x,Del_t
# Read the Parameters from the file
try:
    L=np.genfromtxt('{}'.format(File),skip_header=10, usecols=1, max_rows=1)
    T=np.genfromtxt('{}'.format(File),skip_header=11, usecols=1, max_rows=1)
    Del_x=np.genfromtxt('{}'.format(File),skip_header=12, usecols=1, max_rows=1)
    Del_t=np.genfromtxt('{}'.format(File),skip_header=13, usecols=1, max_rows=1)

except:
    print('A value is missing in Length and Time of Discretization Table')

global Q_i, h_i
# Read initial conditions from the file.
try:
    Q_i=np.genfromtxt('{}'.format(File),skip_header=17, usecols=1, max_rows=1)
    h_i=np.genfromtxt('{}'.format(File),skip_header=18, usecols=1, max_rows=1)
except:
    print('A value is missing in Initial Conditions Table')  


global Chezy,b,S0
# Read the specifications of channel from file.
try:
    Chezy=np.genfromtxt('{}'.format(File),skip_header=22, usecols=1, max_rows=1)     # It is chezy constant of channel.
    b=np.genfromtxt('{}'.format(File),skip_header=23, usecols=1, max_rows=1)         # It is width of channel.
    S0=np.genfromtxt('{}'.format(File),skip_header=24, usecols=1, max_rows=1)     # It is chezy constant of channel.
except:
    print('A value is missing in Channel Specifications Table')          


global tita,kasi
try:
    # Read weighting factors of the scheme
    tita=np.genfromtxt('{}'.format(File),skip_header=28, usecols=1, max_rows=1)     # It is chezy constant of channel.
    kasi=np.genfromtxt('{}'.format(File),skip_header=29, usecols=1, max_rows=1)         # It is width of channel.
except:
    print('A value is missing in Scheme Weighting Factor Table')          
    
global bc_Type, lbc_Type, rbc_Type
# Read the type of boundary conditions from the file
try:
    bc_Type=str(np.genfromtxt('{}'.format(File),skip_header=34, usecols=1, max_rows=1, dtype='str')) 
    lbc_Type=str(np.genfromtxt('{}'.format(File),skip_header=39, usecols=1, max_rows=1, dtype='str'))
    rbc_Type=str(np.genfromtxt('{}'.format(File),skip_header=41, usecols=1, max_rows=1, dtype='str'))
except:
    print('A value is missing in Boundary Conditions Type Tables')          
            
    
global Q_lbc_C,h_lbc_C,Q_rbc_C,h_rbc_C,Q_lbc_V, h_lbc_V,Q_rbc_V,h_rbc_V,t_lbc_V,t_rbc_V
# Read Boundry Conditions from the file
try:
    if bc_Type=='C' and lbc_Type=='Q':
        Q_lbc_C=np.genfromtxt('{}'.format(File),skip_header=47, usecols=1, max_rows=1, delimiter=':')
    if bc_Type=='C' and lbc_Type=='H':
        h_lbc_C=np.genfromtxt('{}'.format(File),skip_header=47, usecols=1, max_rows=1, delimiter=':')
    if bc_Type=='C' and rbc_Type=='Q':
        Q_rbc_C=np.genfromtxt('{}'.format(File),skip_header=51, usecols=1, max_rows=1, delimiter=':')
    if bc_Type=='C' and rbc_Type=='H':
        h_rbc_C=np.genfromtxt('{}'.format(File),skip_header=51, usecols=1, max_rows=1, delimiter=':')
    if bc_Type=='V' and lbc_Type=='Q':
        Q_lbc_V=np.genfromtxt('{}'.format(File),skip_header=58, usecols=1, max_rows=10,)
    if bc_Type=='V' and lbc_Type=='H':
        h_lbc_V=np.genfromtxt('{}'.format(File),skip_header=58, usecols=1, max_rows=10)
    if bc_Type=='V' and rbc_Type=='Q':
        Q_rbc_V=np.genfromtxt('{}'.format(File),skip_header=72, usecols=1, max_rows=10,)
    if bc_Type=='V' and rbc_Type=='H':
        h_rbc_V=np.genfromtxt('{}'.format(File),skip_header=72, usecols=1, max_rows=10)
    if bc_Type=='V':
        t_lbc_V=np.genfromtxt('{}'.format(File),skip_header=58, usecols=0, max_rows=10)
    if bc_Type=='V':
        t_rbc_V=np.genfromtxt('{}'.format(File),skip_header=72, usecols=0, max_rows=10)

except:
    print('A value is missing in Boundary Conditions Tables')          


global Iter, Err
# Read number of iterations and Error from the file.
Iter=np.genfromtxt('{}'.format(File),skip_header=85, usecols=1, max_rows=1)
Err=np.genfromtxt('{}'.format(File),skip_header=89, usecols=1, max_rows=1)

global g, beta
g=9.81
beta=1

# =============================================================================
# # Calculation Part
# =============================================================================
import numpy as np
# Define the constant coefficients (Continuity Equation)
A1j=(-tita/Del_x)
B1j=(b*(1-kasi)/Del_t)
C1j=tita/Del_x
D1j=b*kasi/Del_t

# Calculate number of space steps
M=int(L/Del_x)+1
# Calculate number of time steps
N=int(T/Del_t)+1

# --------------------- Initiating Matrices for initial values ----------------
#------------------------------------------------------------------------------    
# Make a matrix of Q
Mtrx_Q=np.zeros((N,M))
# Insert initial values of Q to matrix
for i in range(M):
    Mtrx_Q[0,i]=Q_i

# Make a matrix of h
Mtrx_h=np.zeros((N,M))
# Insert initial values of h to matrix
for i in range(M):
    Mtrx_h[0,i]=h_i
    
# Make a matrix for the Q and h
Mtrx_Qh=np.zeros((2*M,N))    

# Insert initial Q and h to the matrix
for n in range(N):
    for j in range(M):
        Mtrx_Qh[j*2,n]=Mtrx_Q[n,j]
for n in range(N):
    for j in range(M):
        Mtrx_Qh[j*2+1,n]=Mtrx_h[n,j]
 
# ------------------------- Creating Necessary Matrices ----------------------
#-----------------------------------------------------------------------------            
# Create matrices out of loop
# Create a matrix for E1j
Mtrx_E1=np.zeros((N,M))   
# Create a matrix for the coefficients of the calculations
Mtrx_Coff=np.zeros((N,2*M,2*M)) 
# Create a matrix of free terms
Mtrx_FTerm=np.zeros((2*M,N))
# Create a matrix for hj,n+1/2
Mtrx_hjn12=np.zeros((N,M))
# Create a matrix for area
Mtrx_Area=np.zeros((N,M))
# Create a matrix of Ajn12
Mtrx_Areajn12=np.zeros((N,M))
# Create a matrix for Areaj12n12
Mtrx_Areaj12n12=np.zeros((N,M))
# Create a matrix for A2
Mtrx_A2=np.zeros((N,M))
# Create a matrix for C2
Mtrx_C2=np.zeros((N,M))
# Create a matrix for B2
Mtrx_B2=np.zeros((N,M))
# Create a matrix for D2
Mtrx_D2=np.zeros((N,M))
# Create a matrix for E2
Mtrx_E2=np.zeros((N,M))
   
# -------------------- Inserting Coefficient to Coefficients Matrix ----------
#-----------------------------------------------------------------------------
# Insert the constant coefficients to coefficient matrix
for n in range(N):
    for i in np.arange(1,Mtrx_Coff.shape[1]-2,2):
        Mtrx_Coff[n,i,i-1]=A1j    
        Mtrx_Coff[n,i,i]=B1j    
        Mtrx_Coff[n,i,i+1]=C1j    
        Mtrx_Coff[n,i,i+2]=D1j    

# Insert boundry conditions factors to coefficient matrix
# Insert left boundry conditions factors
for n in range(N):
    if lbc_Type=='Q': 
        Mtrx_Coff[n,0,0]=1

    elif lbc_Type=='H':
        Mtrx_Coff[n,0,1]=1
# Insert right boundry conditions factors to coefficient matrix        
for n in range(N):
    if rbc_Type=='H':
        Mtrx_Coff[n,2*M-1,2*M-1]=1
    elif rbc_Type=='Q': 
        Mtrx_Coff[n,2*M-1,2*M-2]=1

# ----------------------- Insert Boundry Conditions to Free Term Matix--------
#-----------------------------------------------------------------------------
        
# Insert Constant Boundry conditions to Free Term Matrix
if bc_Type=='C':
    # Insert Left Boundry conditions to Free Term Matrix
    for n in range(N):
        if bc_Type=='C':
            if lbc_Type=='Q':
                Mtrx_FTerm[0,n]=Q_lbc_C    
            elif lbc_Type=='H':
                Mtrx_FTerm[0,n]=h_lbc_C
        pass
    # Insert right boundry condition to Free Term Matrix
    for n in range(N):
        if bc_Type=='C':            
            if rbc_Type=='Q':
                Mtrx_FTerm[-1,n]=Q_rbc_C
            elif rbc_Type=='H':
                Mtrx_FTerm[-1,n]=h_rbc_C
        pass

# Make a time series for the boundry conditions
t=np.arange(0,T+Del_t,Del_t)
# Interpolate variable boundry conditions for the time series
if bc_Type=='V':
    if lbc_Type=='Q':
        Q_lbc_V_t=np.interp(t,t_lbc_V,Q_lbc_V)
        pass
    elif lbc_Type=='H':
        h_lbc_V_t=np.interp(t,t_lbc_V,h_lbc_V)
        pass
    if rbc_Type=='Q':
        Q_rbc_V_t=np.interp(t,t_rbc_V,Q_rbc_V)
        pass
    elif rbc_Type=='H':
        h_rbc_V_t=np.interp(t,t_rbc_V,h_rbc_V)
        pass
    
# Insert variable boundry conditions to free term matrix
for n in range(N):
    if bc_Type=='V':
        if lbc_Type=='Q':
            Mtrx_FTerm[0,n]=Q_lbc_V_t[n]
            pass
        elif lbc_Type=='H':
            Mtrx_FTerm[0,n]=h_lbc_V_t[n]
            pass
        if rbc_Type=='Q':
            Mtrx_FTerm[-1,n]=Q_rbc_V_t[n]
            pass
        elif rbc_Type=='H':
            Mtrx_FTerm[-1,n]=h_rbc_V_t[n]
            pass  
    
                
    
#---------------- Next Time Steps Calculations -------------------------------
# ----------------------------------------------------------------------------    
for n in range(N-1):
# Calculate Matrix E1
    for j in range(M-1):
        Mtrx_E1[n,j]=-(1-tita)*(Mtrx_Q[n,j+1]-Mtrx_Q[n,j])/Del_x+b*((1-kasi)*Mtrx_h[n,j]/Del_t+kasi*Mtrx_h[n,j+1]/Del_t)
             
# Insert the values of E1 to the Free Term matrix.
#for n in range(N):
    for j in range(M-1):
        Mtrx_FTerm[j*2+1,n]=Mtrx_E1[n,j] 
  
# Assume next values of h
#for n in range(N-1):
    for j in range(M):
        Mtrx_h[n+1,j]=Mtrx_h[n,j]

# Do the iterations
    for iterations in range(int(Iter)):
    # Calculate the values of matrix of hjn12
    #for n in range(N-1):
        for j in range(M):
            Mtrx_hjn12[n,j]=(Mtrx_h[n+1,j]+Mtrx_h[n,j])/2
    
    # Calculate the values of the Area matrix
    #for n in range(N-1):
        for j in range(M):
            Mtrx_Area[n,j]=Mtrx_h[n,j]*b
            Mtrx_Area[n+1,j]=Mtrx_h[n+1,j]*b
    
    # Calculate the next values of  Mtrx_Areajn12
    #for n in range(N-1):
        for j in range(M):
            Mtrx_Areajn12[n,j]=b*Mtrx_hjn12[n,j]
    
    
    # Calculate the next values of  Mtrx_Areaj12n12
    #for n in range(N-1):
        for j in range(M-1):
            Mtrx_Areaj12n12[n,j]=(Mtrx_Area[n,j]+Mtrx_Area[n,j+1]+Mtrx_Area[n+1,j]+Mtrx_Area[n+1,j+1])/4
         
      
    # Calculate matrix of A2
    #for n in range(N-1):                       
        for j in range(M):
            Mtrx_A2[n,j]=((1-kasi)/Del_t-beta*Mtrx_Q[n,j]/(Del_x*Mtrx_Areajn12[n,j])+g*Mtrx_Areajn12[n,j]*(1-kasi)*abs(Mtrx_Q[n,j])/(Chezy*b*(Mtrx_hjn12[n,j])**(3/2))**2)
    
    # Put A2 in the coefficient matrix
    #for n in range(N-1):
        for j in range(Mtrx_A2.shape[1]-1):
                Mtrx_Coff[n,j*2+2,j*2]=Mtrx_A2[n,j] 
    
    # Calculate matrix of C2
    #for n in range(N-1):                       
        for j in range(M-1):
            Mtrx_C2[n,j]=(kasi/Del_t+beta*Mtrx_Q[n,j+1]/(Del_x*Mtrx_Areajn12[n,j+1])+g*Mtrx_Areajn12[n,j+1]*kasi*abs(Mtrx_Q[n,j+1])/(Chezy*b*(Mtrx_hjn12[n,j+1])**(3/2))**2)
    
    
    # Put C2 in the coefficient matrix
    #for n in range(N-1):
        for j in range(Mtrx_C2.shape[1]-1):
                Mtrx_Coff[n,j*2+2,j*2+2]=Mtrx_C2[n,j] 
     
    # Calculate the next values of the matrix B2
    #for n in range(N-1):
        for j in range(M-1):
            Mtrx_B2[n,j]=-(g*tita/Del_x)*Mtrx_Areaj12n12[n,j]
    
    # Put B2 in the coefficient matrix
    #for n in range(N-1):
        for j in range(Mtrx_B2.shape[1]-1):
                Mtrx_Coff[n,j*2+2,j*2+1]=Mtrx_B2[n,j] 
    
    
    # Calculate the next values of the matrix D2
    #for n in range(N-1):
        for j in range(M-1):
            Mtrx_D2[n,j]=g*tita/Del_x*Mtrx_Areaj12n12[n,j]
    
    
    # Put D2 in the coefficient matrix
    #for n in range(N-1):
        for j in range(Mtrx_D2.shape[1]-1):
                Mtrx_Coff[n,j*2+2,j*2+3]=Mtrx_D2[n,j] 
    
    
    # Calculate matrix of E2
    #for n in range(N-1):
        for j in range(M-1):
            Mtrx_E2[n,j]=(1-kasi)/Del_t*Mtrx_Q[n,j]+kasi/Del_t*Mtrx_Q[n,j+1]-(1-tita)*g*Mtrx_Areaj12n12[n,j]/Del_x*(Mtrx_h[n,j+1]-Mtrx_h[n,j])+g*S0*Mtrx_Areaj12n12[n,j]#+beta*(1-tita)/(Del_x*Mtrx_Area[n,j])*(Mtrx_Q[n,j])**2-beta*(1-tita)/(Del_x*Mtrx_Area[n,j+1])*(Mtrx_Q[n,j+1])**2
            
    # Put values of E2 in the Free Terms Matrix
    #for n in range(N-1):                            # Make it sure!
        for j in np.arange(M-1):
            Mtrx_FTerm[j*2+2,n]=Mtrx_E2[n,j] 
             
        # Inverse the coefficient matrix           
        Mtrx_Inverse=np.linalg.inv(Mtrx_Coff[n,:,:])
        #for j in range(2*M):
        Mtrx_Qh_N=np.dot(Mtrx_Inverse,Mtrx_FTerm[:,n])
        
        # Set the error limit for the calculations
        Error=1
        for i in range(2*M):
            # Make a signal for small error
            if abs(Mtrx_Qh[i,n+1]-Mtrx_Qh_N[i])<=int(Err):          # Made it better                 
                Error=0
            elif abs(Mtrx_Qh[i,n+1]-Mtrx_Qh_N[i])>int(Err):          # Made it better                 
                Error=1
                # Insert the new calculated Q and h to Mtrx_Qh
                #for n in range(1,N):
                for i in range(2*M):
           
                    Mtrx_Qh[i,n+1]=Mtrx_Qh_N[i]
                        
           
                # Insert values from Mtrx_Qh to Mtrx_h and Mtrx_Q
                #for n in range(1,N):
                for j in range(M):
                    Mtrx_Q[n+1,j]=Mtrx_Qh[j*2,n+1]
                    
                #for n in range(1,N):
                for j in range(M):
                    Mtrx_h[n+1,j]=Mtrx_Qh[j*2+1,n+1]
                   
                # Rewrite the constant boundry conditions on the Q and h matrices
                if bc_Type=='C':
                    if  bc_Type=='C' and lbc_Type=='Q':
                        Mtrx_Q[n,0]=Q_lbc_C
                        pass
                    if bc_Type=='C' and lbc_Type=='H':
                        Mtrx_h[n,0]=h_lbc_C
                        pass
                    if bc_Type=='C' and rbc_Type=='Q':
                        Mtrx_Q[n,-1]=Q_rbc_C  
                        pass                      
                    if bc_Type=='C' and rbc_Type=='H':
                       Mtrx_h[n,-1]=h_rbc_C
                       pass
                               
                # Rewrite the variable boundry conditions on the Q and h matrices
                if bc_Type=='V':
                    if  lbc_Type=='Q':
                        Mtrx_Q[n,0]=Q_lbc_V_t[n]
                        pass
                    if lbc_Type=='H':
                        Mtrx_h[n,0]=h_lbc_V_t[n]
                        pass
                    if rbc_Type=='Q':
                        Mtrx_Q[n,-1]=Q_rbc_V_t[n]  
                        pass                      
                    if rbc_Type=='H':
                       Mtrx_h[n,-1]=h_rbc_V_t[n]
                       pass
                   
        # Break the loop if error is less
        if Error==0:
            break 
                  
        # Insert the new calculated Q and h to Mtrx_Qh
        #for n in range(1,N):
        for i in range(2*M):
   
            Mtrx_Qh[i,n+1]=Mtrx_Qh_N[i]
                
   
        # Insert values from Mtrx_Qh to Mtrx_h and Mtrx_Q
        #for n in range(1,N):
        for j in range(M):
            Mtrx_Q[n+1,j]=Mtrx_Qh[j*2,n+1]
            
        #for n in range(1,N):
        for j in range(M):
            Mtrx_h[n+1,j]=Mtrx_Qh[j*2+1,n+1]
           
        # Rewrite the constant boundry conditions on the Q and h matrices
        if bc_Type=='C':
            if  bc_Type=='C' and lbc_Type=='Q':
                Mtrx_Q[n,0]=Q_lbc_C
                pass
            if bc_Type=='C' and lbc_Type=='H':
                Mtrx_h[n,0]=h_lbc_C
                pass
            if bc_Type=='C' and rbc_Type=='Q':
                Mtrx_Q[n,-1]=Q_rbc_C  
                pass                      
            if bc_Type=='C' and rbc_Type=='H':
               Mtrx_h[n,-1]=h_rbc_C
               pass
           
        # Rewrite the variable boundry conditions on the Q and h matrices
        if bc_Type=='V':
            if  lbc_Type=='Q':
                Mtrx_Q[n,0]=Q_lbc_V_t[n]
                pass
            if lbc_Type=='H':
                Mtrx_h[n,0]=h_lbc_V_t[n]
                pass
            if rbc_Type=='Q':
                Mtrx_Q[n,-1]=Q_rbc_V_t[n]  
                pass                      
            if rbc_Type=='H':
               Mtrx_h[n,-1]=h_rbc_V_t[n]
               pass
           
           
# =============================================================================
# # Plotting Part   
# =============================================================================
# Plot the results
import matplotlib.pyplot as plt

 # Plot the results
# Make the space steps for plotting
x1=np.arange(0,L,Del_x)
x2=np.arange(0,L+Del_x,Del_x)

t1=np.arange(0,T,Del_t)
t2=np.arange(0,T+Del_t,Del_t)

fig1=plt.figure(figsize=(12,8))
plt.title('Longitudnal Water Depth Profile in Time', size=20, color='red')
plt.grid()
plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Water Depth (m)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)

for i in np.arange(1,N,int(N/6)):
    try:
        h_plot=Mtrx_h[i,:]
        plt.plot(x1,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
       h_plot=Mtrx_h[i,:]
       plt.plot(x2,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
       pass
plt.legend(loc='upper right')
Name1=input('Give the name to the H graph file:')
fig1.savefig('..\\report\\{}.jpg'.format(Name1))


fig2=plt.figure(figsize=(12,8))
plt.title('Longitudnal Discharge Profile in Time', size=20, color='red')
plt.grid()
plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Discharge (m3/s)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)
for i in np.arange(1,N,int(N/6)):
    try:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x1,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x2,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass
plt.legend(loc='upper right')

# Save the figure
Name2=input('Give the name to the Q graph file:')
fig2.savefig('..\\report\\{}.jpg'.format(Name2))

# Make a figure for both Q and H
fig3=plt.figure(figsize=(18,14))
plt.subplot(2,1,1)
plt.title('1D Sain Venant Equation Solution by Preismmann Scheme', size=30, color='red')

for i in np.arange(1,N,int(N/6)):
    try:
        h_plot=Mtrx_h[i,:]
        plt.plot(x1,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
       h_plot=Mtrx_h[i,:]
       plt.plot(x2,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
       pass
plt.grid()
#plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Water Depth (m)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)
plt.legend(loc='upper right')

plt.subplot(2,1,2)
#plt.title('Longitudnal Discharge Profile in Time', size=20, color='red')
for i in np.arange(1,N,int(N/6)):
    try:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x1,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x2,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass
plt.grid()
plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Discharge (m3/s)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)
plt.legend(loc='upper right')

# Save the figure
Name3=input('Give the name to the Q-h graph file:')
fig3.savefig('..\\report\\{}.jpg'.format(Name3))

# =============================================================================
# # Saving Result Matrices to TXT file Part
# =============================================================================
np.savetxt('..\\report\Matrix_h.txt', Mtrx_h, fmt='%1.1e')
np.savetxt('..\\report\Matrix_Q.txt', Mtrx_Q, fmt='%1.1e')
   


# Do the animations
from celluloid import Camera 

# Animate both at one
fig3=plt.figure(figsize=(18,14))
Camera3=Camera(fig3)
Camera4=Camera(fig3)
plt.subplot(2,1,1)
plt.grid()
#plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Water Depth (m)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)
plt.title('1D Sain Venant Equation Solution by Preismmann Scheme', size=30, color='red')

for i in np.arange(1,N,1):
    try:
        h_plot=Mtrx_h[i,:]
        plt.plot(x1,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
       h_plot=Mtrx_h[i,:]
       plt.plot(x2,h_plot, label='T={0:.0f} sec'.format(i*Del_t))
       pass
    #plt.legend(fig3,'Time={} seconds'.format(i*Del_t), loc='upper center')
    Camera3.snap()
    



plt.subplot(2,1,2)
plt.grid()
plt.xlabel('Chainage (m)', size=15)
plt.ylabel('Discharge (m3/s)', size=15)
plt.xticks(size=12)
plt.yticks(size=12)

for i in np.arange(1,N,1):
    try:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x1,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass

    except ValueError:
        Q_plot=Mtrx_Q[i,:]
        plt.plot(x2,Q_plot, label='T={0:.0f} sec'.format(i*Del_t))
        pass
    Camera4.snap()


global animation3, animation4
animation3=Camera3.animate()
animation4=Camera4.animate()
   








Mtrx_Result=np.ndarray(shape=(2,N,M))
Mtrx_Result[0,:,:]=Mtrx_Q
Mtrx_Result[1,:,:]=Mtrx_h
#return Mtrx_Result
#%%

#a=Calc_FSF('..\\inputs\Transient_Test_A_1.txt')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    