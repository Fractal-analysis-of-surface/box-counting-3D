# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:39:12 2019

@author: Tomáš Ič
"""
#kniznice
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import tkinter
import tkinter as tk
import os.path 
import math
from tkinter import messagebox
from tkinter import filedialog
from tkinter.ttk import *
from tkinter import ttk
from tkinter import *

cond=1
cond1=1
while cond<10:     ## Volbou hodnoty napr. 10 nastavíme podmienku najmenšej prijatelnej hustoty gridu 10 m
    if cond1!=1:    ## Hlasenie podmienky na nevhodnú hustotu gridu súboru
        print("Error (wrong size of grid) ...")
        print("Actual size of grid:" ,cond, "(must be more than 9)",)
        
    cond1+=1    
#Pouzivatelske rozhranie GUI
    a = tkinter.Tk() # Definovanie pola GUI
    a.title("BCdim")# Nazov tabulky
    a.sourceFolder = ''
    a.sourceFile = ''
    a.bb = ''

    def callbackFunc(event):# Funkcia na dopytovanie vybraneho poctu tteracii
        a.bb=comboExample.get()
        print (comboExample.get())
        global rrr
        rrr = a.bb
        
    def chooseFile():# Funkcia vyberu vstupneho suboru
        a.sourceFile = filedialog.askopenfilename(filetypes=('* files','* .txt'))
        label=Label(text=a.sourceFile)
        label.grid(column=3, row=0)  
        label.configure(text=a.sourceFile)

    def close_window (): # Funkcia zatvorenia okna GUI
        a.destroy()  
    
    b_chooseFile = tkinter.Button(a, text = "Load", width = 15, height = 1, \
                              command = chooseFile)# Parametre tlacidla Load
    b_chooseFile.grid(column=1, row=0, sticky = W) # Velkost tlacidla Load

    lbl = tkinter.Label(a, text="Path to file:").grid(column=0, row=0) #Popis
    lbl2 = tkinter.Label(a, text="Number of iteration:").grid(column=0, row=1)#Popis
 
    comboExample = ttk.Combobox(a, values=[1, 2, 3, 4, 5, 6, 7 ], \
                            width=15)# Prednastavené hodnotyodnoty iteracii
    comboExample.grid(column=1, row=1)# Velkost pola s hodnotami iteracii
    comboExample.current(1)
    comboExample.bind("<<ComboboxSelected>>", \
                  callbackFunc)# Zaznamenanie volby iteracii
    comboExample.set("Choose iteration")# Predvolena hodnota v poli iteracii

    frame = Frame(a)
    frame.grid(column=1, row=3)
    button = Button (frame, text = "Get start", command = close_window, \
                 width=15)# Tlacidlo spustenia vypoctu / vypnutie GUI
    button.grid(column=1, row=2)# Parametre tlacidla Get start

    a.mainloop() # Ukoncenie skriptu definujuceho GUI
    print(a.sourceFile )

###############################################################################
# Nacitanie vstupnych dat
    filename=a.sourceFile # Nazov suboru vstupných dát
    adress=os.path.dirname(filename)# Cesta k vstupným datam
    p=int(rrr) # Pocet iteracii  

    print("Loading data ...")
    
    Data=np.loadtxt(filename,dtype=float) ## Do Matice data sa načítajú hodnoty zo vstupného súboru
    
#Spracovanie dát            
    length = len(data) # Pocet prvkov gridu
    i=int(length**(1/2)) # Pocet voxelov pozdlz osi X alebo Y
    mincoord=[0]*3 # Minimalne hodnoty súradníc X,Y,Z 
    mincoord = data.min(axis=0)
    maxcoord=[0]*3 # Maximalne hodnoty súradníc X,Y,Z
    maxcoord = data.max(axis=0)
    vox=(maxcoord[1]-mincoord[1])/((length**(1/2))-1) # Velkosť voxla
    cond=round(vox,2)

S_Matrix=np.zeros(shape = (i,i), dtype = 'float') ## Rastrová matica S_Matrix obsahuje výšky bodov povrchu

if math.ceil((maxcoord[2]-mincoord[2])/vox)<=i:
    iz=i
elif math.ceil((maxcoord[2]-mincoord[2])/vox)>i:
    iz=math.ceil((maxcoord[2]-mincoord[2])/vox)

i = int(iz) 

# Nacitanie udajov o výsje do poli matice S_Matrix
mmm=0
for m in range(0,i,1): # os y
    for n in range(0,i,1): # os x
        S_Matrix[m][n]=data[(mmm*i)+n][2]
        continue
    mmm += 1 
    
print("Creating Matrix ...")
Matrix=np.zeros(shape = (i,i), dtype = 'int') ## Rastrová matica Matrix obsahuje výšky bodov povrchu prepočítané na poschodie voxela
for m in range(0, i, 1):
    for n in range(0, i, 1):
        if (S_Matrix[m][n]-mincoord[2])%vox==0 and S_Matrix[m][n]==mincoord[2]:
            Matrix[m][n]= ((S_Matrix[m][n]-mincoord[2])//vox)+1
            continue 
        if (S_Matrix[m][n]-mincoord[2])%vox ==0 and S_Matrix[m][n]!=mincoord[2]:
            Matrix[m][n]= ((S_Matrix[m][n]-mincoord[2])//vox)
            continue 
        if (S_Matrix[m][n]-mincoord[2])%vox !=0 :
            Matrix[m][n]= ((S_Matrix[m][n]-mincoord[2])//vox)+1
            continue            
            
N_s = [0]*p 
saa = [0]*p
dif_coord=maxcoord[1]-mincoord[1] # Rozdiel maimalnych a minimalnych suradnic

##Zaciatok iteracneho vypoctu
for q in range(0, p, 1):    # Hodnotu iteracie predstavuje q
    pp = (2)**(q+1) # Pocet boxov v smere X alebo Y
    dbox = (dif_coord+vox)/pp # Velkost boxu v iteracii

    if ((maxcoord[2]-mincoord[2])%dbox)==0:
        zpp=(2)**(q+1)
    elif ((maxcoord[2]-mincoord[2])%dbox)!=0:
        zpp=((maxcoord[2]-mincoord[2])//dbox)+1
    zpp = int(zpp)
    

    index1 = (i/pp)
    saa [q] = (1 / pp) ## Pomer zmensenia velkosti boxu voci velkosti boxu v 0. iteracii
    
    print ("Iteration number: "+ str(q))
# Generovanie matic MinMatrix & MaxMatrix    
    MinMatrix=np.zeros(shape = (pp, pp), dtype = 'float') # Polia matice obsahuju hodnoty najvyssich pochodí plnych boxov v ramci stlpca
    MaxMatrix=np.zeros(shape = (pp, pp), dtype = 'float') # Polia matice obsahuju hodnoty najnizsich pochodí plnych boxov v ramci stlpca
    
    if index1==int(index1):
        for k in range(0,pp,1):
            for l in range(0,pp,1):
                indexa=Matrix[index1*l][index1*k]
                indexb=Matrix[index1*l][index1*k]
                for m in range(0,index1,1):
                    for n in range(0,index1,1):
                        if indexa<Matrix[(index1*l)+n][(index1*k)+m]:
                            indexa=Matrix[(index1*l)+n][(index1*k)+m]
                        if indexb>Matrix[(index1*l)+n][(index1*k)+m]:
                            indexb=Matrix[(index1*l)+n][(index1*k)+m]
                MinMatrix[l][k]= indexb
                MaxMatrix[l][k]= indexa
                
    if index1!=int(index1):
        for k in range(0,pp,1):
            for l in range(0,pp,1):
                indexa=Matrix[math.ceil(k*index1)][math.ceil(l*index1)]
                indexb=Matrix[math.ceil(k*index1)][math.ceil(l*index1)]
                m=math.ceil(k*index1)
                while m<(index1*(k+1)):
                    n=math.ceil(l*index1)
                    while n<(index1*(l+1)):                       
                        if indexa<Matrix[n][m]:
                            indexa=Matrix[n][m]  
                        if indexb>Matrix[n][m]:
                            indexb=Matrix[n][m]
                        n=n+1
                    m=m+1
                MaxMatrix[l][k]= indexa
                MinMatrix[l][k]= indexb
          
###############################################################################
# Napocitavanie plnych boxov     
    for j in range(0, int(zpp), 1):   # Posun v smere osi Z  
        for k in range(0, int(pp), 1): # Posun v smere osi Y
            for l in range(0, int(pp), 1): # Posun v smere osi Y
                
                ### Podmienky na hladanie oplnzch boxov
                if MinMatrix[l][k]*vox<=(j+1)*dbox and  \
                    MinMatrix[l][k]*vox>(j)*dbox:
                    N_s [q] = N_s [q] + 1 
                    continue
                if MaxMatrix[l][k]*vox<=(j+1)*dbox and \
                    MaxMatrix[l][k]*vox>(j)*dbox:
                    N_s [q] = N_s [q] + 1
                    continue
                
                ### Podmienky na hladanie medzilahlych plnych boxov
                if MaxMatrix[l][k]*vox>(j+1)*dbox and \
                    MinMatrix[l][k]*vox<(j)*dbox:
                    N_s [q] = N_s [q] + 1
                    continue
                
                #stred
                if l!=0 and k!=0 and l!=pp-1 and k!=pp-1 and \
                    MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k+1])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                    
                #lavo hore
                if l==0 and k==0 and MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l+1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k+1])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #lavo dole
                if l==0 and k==pp-1 and MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #pravo hore
                if l==pp-1 and k==0 and MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox):
                    N_s [q] = N_s [q] + 1 
                    continue
                #pravo dole
                if l==pp-1 and k==pp-1 and MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #lava hrana
                if l==0 and k!=0 and k!=pp-1 and \
                    MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k+1])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #prava hrana
                if l==pp-1 and k!=0 and k!=pp-1 and \
                    MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #horna hrana
                if l!=0 and l!=pp-1 and k==0 and \
                    MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k+1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k+1])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
                #dolna hrana
                if l!=0 and l!=pp-1 and k==pp-1 and \
                    MinMatrix[l][k]*vox>(j)*dbox and \
                    ((j)*dbox>=(MaxMatrix[l-1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k-1])*vox or \
                    (j)*dbox>=(MaxMatrix[l-1][k])*vox or \
                    (j)*dbox>=(MaxMatrix[l+1][k])*vox):
                    N_s [q] = N_s [q] + 1
                    continue
    
                    
                    
# Vypocet box-counting dimenzie
N = [1] # Pole s ulozenym poctom plnych boxov v ramci iteracie 
ss = [1] # Pole s ulozenymi velkostami boxov v ramci iteracie
N.extend(N_s)
ss.extend(saa)
 
# Vypocet suradnic bodov v log-log grafe
x = abs (np.log (ss))
y = np.log (N)

# Vypocet fraktalnej dimenzie pre jednotlive iteracie
fd = [0]*p
for w in range(0,p,1):
    fd [w]= (y[w+1]-y[w])/(x[w+1]-x[w])   
D = np.asarray(fd)   
         
#####################################################
            
# MNŠ - linearna regresia
A = np.vstack([x, np.ones(len(x))]).T
m1, cc = np.linalg.lstsq(A, y)[0] # hodnota m1 je celkova fraktalna dimenzia

# Graficke znazornenie log-log grafu
fig = plt.figure()
plt.plot(x,y, 'ro', x, m1*x+cc,'b--')
plt.ylabel('log N(s)')
plt.xlabel('log(1/s)')
plt.suptitle('Log-log graf určenia fraktálnej dimenzie', fontsize=16)
plt.text(1.0, 1.0, r'$D_B$' ' = ''{:0.2f}'.format(m1) , {'color': 'black', \
         'fontsize': 12})
plt.show()
fig.savefig(adress +"/"+'Loglog.png') # Ulozenie grafu na adrese vstupnych dat

# Ulozenie suboru s vysledkami jednotlivych iteracii vo formate *.txt
index2 = [0]*p
for t in range(0,p,1):
    index2 [t] = t+1
    
zip(index2,D)

import csv
with open (adress +"/" +'boxcouting.txt', 'w+') as func:
    writer = csv.writer(func,delimiter='\t')
    writer.writerows(zip(index2,D))
