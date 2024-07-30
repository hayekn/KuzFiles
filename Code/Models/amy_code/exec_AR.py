from final_graphs_AR import *
y_sol = xppaut_model(t, y0, noise=True) #Solve initially so that all graphs come from the same simulation

# subplots(y_sol,time = 120, noise=False, anim=False)
        # ARGUMENTS: time=120, noise=False, save_sub=False, save_ani = False, anim=False
        # DESCRIPTION: Create a subplot of entire system's activity; option for animation to cycle through simulations with various parameter values

# td_vect(t,y0, save=False)
        #DESCRIPTION: create a 3-D vector field and trajectory. x and y are binge/seek activity. z is time (up to 40)

# ind_plots(y_sol['Int'], 'Subcort')
        #ARGUMENTS: y, graph, csTOvta=csTOvta, t=t, save=False
        #DESCRIPTION: individual plots for: PFC, Insula, STR, VTA, Alc, Param, Cortex, Subcort

# DA_graphs(1.35, 1.5, 1.65, 2, save=False)
        #ARGUMENTS: fail, low, med, high DA conc, save
        #DESCRIPTION: Individual nac, vta, alc, subplots of combined, DA effect on total alc and peak nac 

# td_vect_ani(t,y0)
        #DESCRIPTION: 3d vector field animation in which the trajectory moves forward in time. 

# param_array = np.linspace(1,4,1000) #CS to VTA
# da_animation(param_array, save = False)
        #ARGUMENTS: param array for CS to VTA connection
        #DESCRIPTION: animation showing VTA, NAc and Alc activity as the csTOvta strength is altered

# vector_field(y0,y_traj,t=np.linspace(0,50,1000), save=False)
        #ARGUMENTS: y_traj = initial conditions for trajectory
        #DESCRIPTION: creates an animation of the vector fields, nullclines, and trajectory of seek binge system