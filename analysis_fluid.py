import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from tqdm import tqdm


#pdb.set_trace()
#table = ssv.loadf('0030_m.ssv')
maximum_error = 0
for index in tqdm(range(0,3000,30)):
    

    step_string = "%04d" % (index,)
    data1 = pd.read_csv('../result_collection/rocket_reflection_top_bottom/reference_solution_rocket_4000_600_push_700/'+step_string+'_m.ssv',sep='\t')
    data2 = pd.read_csv('../result_collection/rocket_reflection_top_bottom/backtracking_of_error/reference_1500_600_push_700/'+step_string+'_m.ssv',sep='\t')

    vx1 = data1["ux"]
    vy1 = data1["uy"]

    pressure1 = data1["press"]
    rho1 = data1["rho"]
    obstacle = data1["obsval"]

    xpos1 = data1["x"]
    ypos1 = data1["y"]
    xposn = xpos1.to_numpy()
    sizex = int(xposn[-1])+1
    sizey = (xpos1.shape[0])//sizex
    xposn.resize((sizey,sizex))

    ux_numpy1 = vx1.to_numpy()
    uy_numpy1 = vy1.to_numpy()
    obsval = obstacle.to_numpy()
    obsval.resize((sizey,sizex))

    ux_numpy1.resize((sizey,sizex))
    uy_numpy1.resize((sizey,sizex))

    vx2 = data2["ux"]
    vy2 = data2["uy"]

    pressure2 = data2["press"]
    rho2 = data2["rho"]
    obstacle = data1["obsval"]

    xpos2 = data2["x"]
    ypos2 = data2["y"]

    ux_numpy2 = vx2.to_numpy()
    uy_numpy2 = vy2.to_numpy()

    #print("Test the shape:", ux_numpy2.size)
    obsval = obstacle.to_numpy()
    obsval.resize((sizey,sizex))

    ux_numpy2.resize((sizey,sizex))
    uy_numpy2.resize((sizey,sizex))



    diff = (ux_numpy2-ux_numpy1)*(ux_numpy2-ux_numpy1) + (uy_numpy2-uy_numpy1)*(uy_numpy2-uy_numpy1) 

    x = np.linspace(0,599,600)
    y = np.linspace(0,599,600)

    X,Y = np.meshgrid(x,y)
    
    err = (diff.sum())**0.5
    if(err>maximum_error):
        maximum_error = err

    #plt.title("Plot of the difference after "+step_string+" time steps.\n The maximal error is: "+str(err))
    #plt.pcolormesh(X,Y,diff,cmap= cm.plasma)
    #plt.savefig("../result_collection/pml_with_stable_fluid/difference_rocket/difference_step_"+step_string+".png")
    #pdb.set_trace()

print("The maximum error is: ",maximum_error)

"""
plt.quiver(ux,uy, scale=1.5, headwidth=1.)

#plt.savefig("pictures/obstacle_"+step_string+"_m.png", format="png", dpi=200)
#plt.show()
print("Current step: "+step_string)
"""