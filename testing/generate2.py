import numpy as np

#I'm dumb, so I'm reducing the problem to 2D so I can see what's happening

##Make fake data
# an N x 5 array containing a regular mesh representing the stimulus params
stim_params=np.mgrid[10:25,20:22].reshape(2,-1).T

# an N x 3 array representing the output values for each simulation run
stimnum=15*2
output_vals=np.arange(stimnum*3).reshape(stimnum,3)

# shuffle the rows for a bit of added realism
shuf=np.random.permutation(stim_params.shape[0])

stim_params=stim_params[shuf]
output_vals=output_vals[shuf]


##Now we have to arrays, one with the independent variables (stim_params) that
##will represent T and mu values, the other (output_vals) is all the measurements
##you made.
##You can use lexical sort to get the indices of the independent variables in the
##right order, then apply the same indexes to the measurement array

# get the number of unique values for each stimulus parameter
#Due to float precision, you might have to round to the nearest decimal via round
params_shape=tuple(np.unique(col).shape[0] for col in stim_params.T)

# get the set of row indices that will sort the stimulus parameters in ascending
# order, starting with the final column
indx=np.lexsort(stim_params[:,::-1].T)

# sort and reshape the stimulus parameters:
sorted_params=stim_params[indx].T.reshape((2,)+params_shape)

# sort and reshape the output values
sorted_output=output_vals[indx].T.reshape((3,)+params_shape)


###What do the dimensions mean?
## array of stimulus parameters, with dimensions (n_params, p1, p2, p3, p4, p5)
#print(sorted_params.shape)
#
## to check that the sorting worked as expected, we can look at the values of the
## 5th parameter when all the others are held constant at 0:
#print(sorted_params[4,0,0,0,0,:])
#
## ... and the 1st parameter when we hold all the others constant:
#print(sorted_params[0,:,0,0,0,0])
#
## ... now let the 1st and 2nd parameters covary:
#print(sorted_params[:2, :, :, 0, 0, 0])
#
###The same indexing logic applies to the sorted simulation outputs:
## array of outputs, with dimensions (n_outputs, p1, p2, p3, p4, p5)
#print(sorted_output.shape)
#
## the first output variable whilst holding the first 4 simulation parameters
## constant at 0:
#print(sorted_output[0, 0, 0, 0, 0, :])
