import thermoplotting as tp
import matplotlib.pyplot as plt
import os
import numpy
import glob

refT=1200
maxT=2200
minT=10
#maxmu=-0.36486
maxmu=-0.388
minmu=-1.7274
Tstep=5
mustep=0.005

#Translate .txt to .thin for reading. Careful with what you sampled!! casm outputs different things depending
#on what you're sampling

dirs=["heating", "cooling", "coolref"]

#for directory in dirs:
#    print "Preparing input for "+directory+"..."
#    os.chdir("./"+directory)
#    names=glob.glob("*.txt")
#    tp.gc.runs.prepare_input(names,2)
#    os.chdir("./..")

#Read in all the data from .thin files

datablocks={}

for directory in dirs:
    print "Loading data for "+directory+"..."
    os.chdir("./"+directory)
    block=tp.thermoio.clobber()
    os.chdir("./..")
    datablocks[directory]=block

print ""


#Trim off edges of data. Cooling and heating span the same space, while coolref
#is one temperature increment above, and has the chemical potential going much
#further. For a T-mu diagram:
#
#   HHHHHHHHHHHH            CCCCCCCCCCCC            RRRRRRRRRRRRRRRRRRRRRRRRRRRRR
#   HHHHHHHHHHHH            CCCCCCCCCCCC
#   HHHHHHHHHHHH            CCCCCCCCCCCC
#   HHHHHHHHHHHH            CCCCCCCCCCCC
#

print "Trim and concatenate high temperature reference..."
datablocks["coolref"]=tp.gc.runs.slice_copy_T(datablocks["coolref"], maxT, 2)
datablocks["coolref"]=tp.gc.runs.ceil_copy_mu(datablocks["coolref"],1, minmu, 2)
remainingmuref=tp.gc.runs.slice_copy_T(datablocks["cooling"], maxT, 2)
remainingmuref=tp.gc.runs.ceil_copy_mu(remainingmuref, 1, maxmu, 2)
remainingmuref=tp.gc.runs.floor_copy_mu(remainingmuref, 1, minmu, 2)
datablocks["coolref"]=numpy.vstack((datablocks["coolref"],remainingmuref))

for block in ["heating","cooling"]:
    print "Trimming data for "+block+" block"
    datablocks[block]=tp.gc.runs.ceil_copy_mu(datablocks[block], 1, maxmu, 2)
    datablocks[block]=tp.gc.runs.floor_copy_mu(datablocks[block], 1, minmu, 2)
    datablocks[block]=tp.gc.runs.ceil_copy_T(datablocks[block], maxT, 2)
    datablocks[block]=tp.gc.runs.floor_copy_T(datablocks[block], minT, 2)
print ""

#Organize the data so that it can be properly integrated
print "Sorting data for high temperature reference integration step..."
datablocks["coolref"]=tp.gc.runs.sort_mu_copy(datablocks["coolref"],1,2)
#datablocks["cooling"]=tp.gc.runs.sort_T_copy(datablocks["cooling"],2)
#datablocks["cooling"]=tp.gc.access.backwards(datablocks["cooling"])
#datablocks["heating"]=tp.gc.runs.sort_T_copy(datablocks["heating"],2)
print ""


#Fill integration values for the high temperature reference
print "Integrating for high temperature reference..."
startref=tp.gc.access.energy(datablocks["coolref"],2)[0]
tp.gc.integrate.mu(datablocks["coolref"], 1, startref, 2)

#plt.scatter(tp.gc.access.mu(datablocks["coolref"], 1, 2), tp.gc.access.free_energy(datablocks["coolref"],2))
#plt.show()


#Integrate down the temperature
endmu=max(tp.gc.access.mu(datablocks["cooling"],1,2))
startmu=min(tp.gc.access.mu(datablocks["cooling"],1,2))

intercool=[]
interheat=[]

print "Prepare for integration..."
muval=startmu
criticals=[]
while muval<maxmu:

    coolrefslice=tp.gc.runs.slice_copy_mu(datablocks["coolref"],1,muval,2)
    coolrefvalue=tp.gc.access.free_energy(coolrefslice,2)

    #print muval
    #print coolrefvalue

    coolrun=tp.gc.runs.slice_copy_mu(datablocks["cooling"], 1, muval, 2)
    coolrun=tp.gc.runs.sort_T_copy(coolrun, 2)
    coolrun=tp.gc.access.backwards(coolrun)


    heatrun=tp.gc.runs.slice_copy_mu(datablocks["heating"],1,muval,2)
    heatrun=tp.gc.runs.sort_T_copy(heatrun,2)
    heatrefvalue=tp.gc.access.low_T_expansion(heatrun,2)[0]
    

    tp.gc.integrate.beta(coolrun, coolrefvalue, 2)
    tp.gc.integrate.beta(heatrun, heatrefvalue, 2)

    coolrun=tp.gc.access.backwards(coolrun)

    #plt.figure(1)
    #plt.scatter(tp.gc.access.temperature(coolrun, 2), tp.gc.access.free_energy(coolrun,2))
    #plt.scatter(tp.gc.access.temperature(heatrun, 2), tp.gc.access.free_energy(heatrun,2), color='red')
    #plt.show()

    print muval
    #plt.figure(0)
    #coolheatcap=tp.gc.runs.heat_capacity(coolrun,2)
    #plt.scatter(tp.gc.access.temperature(coolrun,2), coolheatcap)
    #plt.figure(1)
    #plt.scatter(tp.gc.access.species(coolrun,1,2)/2048, tp.gc.access.temperature(coolrun,2))
    #plt.show()

    plt.figure(0)
    gibbsheat=tp.gc.runs.gibbs(heatrun,2)
    #plt.scatter(tp.gc.access.species(heatrun,1,2)/2048,gibbsheat, color='red')
    plt.scatter(tp.gc.access.temperature(heatrun,2)/2048,gibbsheat, color='red')
    #gibbscool=tp.gc.runs.gibbs(coolrun,2)
    #plt.scatter(tp.gc.access.species(coolrun,1,2)/2048,gibbscool)
    plt.show()


    intercool.append(tp.gc.runs.sort_T_copy(coolrun, 2))
    interheat.append(tp.gc.runs.sort_T_copy(heatrun, 2))

    criticals=criticals+(tp.gc.integrate.cross(coolrun,heatrun, 2, 0.0001))

    muval+=5*mustep
print "Done integrating."

datablocks["cooling"]=numpy.vstack(intercool)
datablocks["heating"]=numpy.vstack(interheat)

criticals=numpy.vstack(criticals)
plt.figure(0)
plt.scatter(tp.gc.access.species(criticals,1,2)/2048,tp.gc.access.temperature(criticals,2))
plt.figure(1)
plt.scatter(tp.gc.access.mu(criticals,1,2),tp.gc.access.temperature(criticals,2))
plt.show()
