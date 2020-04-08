"""Module contains functions for greating the turbine design GUI"""

# Import required modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from turbine import turbine, optimise
from prof_gen import Profile

# Setting global variables
Po1_ = 145*10**5
To1_ = 950
mdot_ = 16
Omega_ = 6782*2*np.pi/60
W_ = 17*10**6
t_ = 0.0003
g_ = 0.0003
phi_ = 0.4
psi_ = 1
Lambda_ = 0.5
AR_ = 1
dho_ = 1
n_ = 10
ptc_ = 1
ain_ = 0
init_turbine = turbine(Po1_, To1_, mdot_, Omega_, W_, t_, g_, phi_, psi_,
                       Lambda_, AR_, dho_, n_, ptc_, ain_)
new_turbine = init_turbine


def b2b_data(turbine_data):
    """Return blade-to-blade profiles for the whole turbine"""

    # Extract values from the turbine function output
    angs = [[i[0], i[1], i[2], i[4]] for i in turbine_data[10]]
    chords = [[i[4], i[5]] for i in turbine_data[5]]
    ptcs = [[i[6], i[7]] for i in turbine_data[5]]
    tTE = turbine_data[11][5]
    n = turbine_data[11][12]
    # Number of points per blade plot
    points = 500
    TE_points = 200
    # Scale the outputs for m/cm/mm etc
    scale = 1000
    # Initialise values
    L = 0
    # Array size set by the maximum number of stages
    max_stages = 30
    x = np.zeros([4*max_stages, 2*points+TE_points])
    y = np.zeros([4*max_stages, 2*points+TE_points])
    # Loop over every stage
    for i in range(n):
        # Extract stage parameters
        a1, a2, b2, b3 = angs[i]
        Cxst, Cxro = chords[i]
        ptcst, ptcro = ptcs[i]
        # Determine start point depending on b2b spacing
        if i == 0:
            x0st = 0
            x0ro = Cxst*1.25+Cxro*0.25
        else:
            x0st = L+Cxst*0.25
            x0ro = L+Cxst*1.5+Cxro*0.25
        # Calculate the length along the turbine
        if i == 0:
            L += Cxst*1.25+Cxro*1.5
        else:
            L += Cxst*1.5+Cxro*1.5
        # Pass parameters to profile function for stator
        *_, X, Y = Profile(a1, a2, tTE, Cxst, points, TE_points)
        # Store results
        x[4*i] = [j+x0st for j in X]
        y[4*i] = Y
        x[4*i+1] = [j+x0st for j in X]
        y[4*i+1] = [j+Cxst*ptcst for j in Y]
        # Pass parameters to profile function for rotor
        *_, X, Y = Profile(b2, b3, tTE, Cxro, points)
        # Store results
        x[4*i+2] = [j+x0ro for j in X]
        y[4*i+2] = Y
        x[4*i+3] = [j+x0ro for j in X]
        y[4*i+3] = [j+Cxro*ptcro for j in Y]

    return scale*x, scale*y, scale, max_stages


def b2b_plot(turbine_data):
    """Plot blade-to-blade profiles for the whole turbine"""

    # Use data from turbine function to get data for profiles
    data = b2b_data(turbine_data)
    # Find the number of stages to plot
    n = turbine_data[11][12]
    # Plot each bit of data, which are in the form of 2D arrays
    plt.figure()
    plt.plot(data[0][:4*n].T, data[1][:4*n].T, 'black', lw=2)
    plt.axis('equal')
    plt.xlabel('Distance along turbine (mm)')
    plt.ylabel('Tangential distance (mm)')
    plt.show()


def b2b_variable(turbine_data=init_turbine):
    """Plot blade-to-blade profiles with variable loadings"""

    # Set the global variables
    global new_turbine
    new_turbine = turbine_data
    # Save the start for resetting
    original_turbine = new_turbine
    # Extract new parameters
    Po1, To1, mdot, Omega, W, t, g, phi, psi, Lambda,\
        AR, dho, n, ptc, ain, gas = new_turbine[11]
    phi01 = phi[0]
    phi02 = phi[-1]
    psi01 = psi[0]
    psi02 = psi[-1]
    Lambda01 = Lambda[0]
    Lambda02 = Lambda[-1]
    dho01 = dho[0]
    dho02 = dho[-1]
    AR01 = AR[0]
    AR02 = AR[-1]
    ptc01 = ptc[0]
    ptc02 = ptc[-1]
    # Use function to get data for b2b plot
    data = b2b_data(new_turbine)
    scale = data[2]
    max_stages = data[3]
    Cxst1 = new_turbine[5][0][4]
    Cxron = new_turbine[5][-1][5]
    Cxmax = np.amax([i[4] for i in new_turbine[5]])
    # Initialise plots
    _, ax = plt.subplots()
    plt.subplots_adjust(left=0.06, right=0.86, top=0.85, bottom=0.2)
    plt.subplots_adjust(bottom=0.25)
    blades = plt.plot(data[0].T, data[1].T, 'black', lw=2)
    plt.axis('equal')
    ax.set_xbound(np.amin(data[0][:4*n])-scale*Cxst1,
                  np.amax(data[0][:4*n])+scale*Cxron)
    ax.set_ybound(np.amin(data[1][:4*n])-scale*Cxmax,
                  np.amax(data[1][:4*n])+scale*Cxmax)
    plt.xlabel('Distance along turbine (mm)')
    plt.ylabel('Tangential distance (mm)')
    # Create axes for position of sliders [left, bottom, width, height]
    axphi1 = plt.axes([0.06, 0.16, 0.32, 0.02])
    axphi2 = plt.axes([0.54, 0.16, 0.32, 0.02])
    axpsi1 = plt.axes([0.06, 0.13, 0.32, 0.02])
    axpsi2 = plt.axes([0.54, 0.13, 0.32, 0.02])
    axLambda1 = plt.axes([0.06, 0.1, 0.32, 0.02])
    axLambda2 = plt.axes([0.54, 0.1, 0.32, 0.02])
    axdho1 = plt.axes([0.06, 0.07, 0.32, 0.02])
    axdho2 = plt.axes([0.54, 0.07, 0.32, 0.02])
    axAR1 = plt.axes([0.06, 0.04, 0.32, 0.02])
    axAR2 = plt.axes([0.54, 0.04, 0.32, 0.02])
    axptc1 = plt.axes([0.06, 0.01, 0.32, 0.02])
    axptc2 = plt.axes([0.54, 0.01, 0.32, 0.02])
    # Create sliders for each variable
    sphi1 = Slider(axphi1, '$\\phi_1$', 0.1, 1.0, valinit=phi01)
    sphi2 = Slider(axphi2, '$\\phi_2$', 0.1, 1.0, valinit=phi02)
    spsi1 = Slider(axpsi1, '$\\psi_1$', 0.2, 3.0, valinit=psi01)
    spsi2 = Slider(axpsi2, '$\\psi_2$', 0.2, 3.0, valinit=psi02)
    sLambda1 = Slider(axLambda1, '$\\Lambda_1$', 0.01, 0.99, valinit=Lambda01)
    sLambda2 = Slider(axLambda2, '$\\Lambda_2$', 0.01, 0.99, valinit=Lambda02)
    sdho1 = Slider(axdho1, '$\\Delta h_{01}$', 0.5, 3, valinit=dho01)
    sdho2 = Slider(axdho2, '$\\Delta h_{02}$', 0.5, 3, valinit=dho02)
    sAR1 = Slider(axAR1, '$AR_1$', 0.2, 5.0, valinit=AR01)
    sAR2 = Slider(axAR2, '$AR_2$', 0.2, 5.0, valinit=AR02)
    sptc1 = Slider(axptc1, '$p/C_{x1}$', 0.3, 1.5, valinit=ptc01)
    sptc2 = Slider(axptc2, '$p/C_{x2}$', 0.3, 1.5, valinit=ptc02)
    # Can probably put p/Cx in a box as it wont vary
    # Text box showing turbine efficiency
    effax = plt.axes([0.885, 0.65, 0.087, 0.04])
    eff = TextBox(effax, '',
                  'Efficiency: {}%'.format(np.round(100*new_turbine[0], 2)),
                  color='1.0')
    # Text box showing mass
    massax = plt.axes([0.885, 0.55, 0.087, 0.04])
    mass = TextBox(massax, '',
                   'Mass: {} kg'.format(np.round(new_turbine[2], 2)),
                   color='1.0')
    # Text box showing total blades
    nbax = plt.axes([0.885, 0.45, 0.087, 0.04])
    nb = TextBox(nbax, '',
                 'No. Blades: {}'.format(int(new_turbine[6])),
                 color='1.0')
    # Text box showing axial rotor force
    fxax = plt.axes([0.885, 0.35, 0.087, 0.04])
    fx = TextBox(fxax, '',
                 'Rotor F$_x$: {} kN'.format(np.round(0.001*new_turbine[13],
                                                      2)),
                 color='1.0')
    # Text box showing maximum angle
    if new_turbine[12][0]:
        col = 'r'
    else:
        col = 'g'
    angle_maxax = plt.axes([0.8735, 0.75, 0.11, 0.04])
    angle_max = TextBox(angle_maxax,
                        '', 'Maximum angle: {}º'.format(
                            np.round(new_turbine[12][1], 2)),
                        color=col, hovercolor=col)
    # Text boxes that allow for other inputs to be changed
    Po1ax = plt.axes([0.14, 0.9, 0.03, 0.04])
    To1ax = plt.axes([0.22, 0.9, 0.03, 0.04])
    mdotax = plt.axes([0.3, 0.9, 0.03, 0.04])
    Omegaax = plt.axes([0.38, 0.9, 0.035, 0.04])
    Wax = plt.axes([0.46, 0.9, 0.03, 0.04])
    tax = plt.axes([0.54, 0.9, 0.03, 0.04])
    gax = plt.axes([0.62, 0.9, 0.03, 0.04])
    nax = plt.axes([0.7, 0.9, 0.03, 0.04])
    a0ax = plt.axes([0.78, 0.9, 0.03, 0.04])
    Po1box = TextBox(Po1ax, '$P_{01}$ (bar) ',
                     '{}'.format(Po1/10**5), color='1.0', label_pad=0.15)
    To1box = TextBox(To1ax, '$T_{01}$ (K) ',
                     '{}'.format(To1), color='1.0', label_pad=0.15)
    mdotbox = TextBox(mdotax, '$\\dot{m}$ (kg/s) ',
                      '{}'.format(mdot), color='1.0', label_pad=0.15)
    Omegabox = TextBox(Omegaax, '$\\Omega$ (rpm) ',
                       '{}'.format(np.round(Omega*60/(2*np.pi), 0)),
                       color='1.0', label_pad=0.15)
    Wbox = TextBox(Wax, '$\\dot{W}$ (MW) ',
                   '{}'.format(W/10**6), color='1.0', label_pad=0.15)
    tbox = TextBox(tax, '$t_{TE}$ (mm) ',
                   '{}'.format(t*10**3), color='1.0', label_pad=0.15)
    gbox = TextBox(gax, '$g$ (mm) ',
                   '{}'.format(g*10**3), color='1.0', label_pad=0.15)
    nbox = TextBox(nax, '$n_{stages}$ ',
                   '{}'.format(n), color='1.0', label_pad=0.15)
    a0box = TextBox(a0ax, '$\\alpha_{in}$ (º) ',
                    '{}'.format(ain), color='1.0', label_pad=0.15)

    # Update the figure and text when the sliders are changed
    def update(val):
        """Update the plots"""
        # Declare the turbine data to be global
        global new_turbine
        # Get data from inputs
        Po1 = float(Po1box.text)*10**5
        To1 = float(To1box.text)
        mdot = float(mdotbox.text)
        Omega = float(Omegabox.text)*2*np.pi/60
        W = float(Wbox.text)*10**6
        t = float(tbox.text)/10**3
        g = float(gbox.text)/10**3
        n = int(float(nbox.text))
        ain = float(a0box.text)
        Po1box.stop_typing()
        To1box.stop_typing()
        mdotbox.stop_typing()
        Omegabox.stop_typing()
        Wbox.stop_typing()
        tbox.stop_typing()
        gbox.stop_typing()
        nbox.stop_typing()
        a0box.stop_typing()
        phi = [sphi1.val, sphi2.val]
        psi = [spsi1.val, spsi2.val]
        Lambda = [sLambda1.val, sLambda2.val]
        dho = [sdho1.val, sdho2.val]
        AR = [sAR1.val, sAR2.val]
        ptc = [sptc1.val, sptc2.val]
        # Recalculate turbine performance
        new_turbine = turbine(Po1, To1, mdot, Omega, W, t,
                              g, phi, psi, Lambda, AR, dho, n, ptc, ain, gas)
        # Get new profile data
        data = b2b_data(new_turbine)
        # Update plots
        for i in range(max_stages):
            for j in range(4):
                blades[4*i+j].set_xdata(data[0].T[:, 4*i+j])
                blades[4*i+j].set_ydata(data[1].T[:, 4*i+j])
        # Update the figure
        Cxst1 = new_turbine[5][0][4]
        Cxron = new_turbine[5][-1][5]
        Cxmax = np.amax([i[4] for i in new_turbine[5]])
        ax.set_xbound(np.amin(data[0][:4*n])-scale*2*Cxst1,
                      np.amax(data[0][:4*n])+scale*2*Cxron)
        ax.set_ybound(np.amin(data[1][:4*n])-scale*Cxmax,
                      np.amax(data[1][:4*n])+scale*Cxmax)
        # Update the efficiency and angle boxes
        eff.set_val('Efficiency: {}%'.format(np.round(100*new_turbine[0], 2)))
        eff.stop_typing()
        angle_max.set_val(
            'Maximum angle: {}º'.format(np.round(new_turbine[12][1], 2)))
        angle_max.stop_typing()
        fx.set_val(
            'Rotor F$_x$: {} kN'.format(np.round(0.001*new_turbine[13], 2)))
        fx.stop_typing()
        mass.set_val('Mass: {} kg'.format(np.round(new_turbine[2], 2)))
        mass.stop_typing()
        nb.set_val('No. Blades: {}'.format(int(new_turbine[6])))
        nb.stop_typing()
        a0box.set_val(np.round(new_turbine[11][-1], 1))
        # Change the colour of the angle box if needed
        if new_turbine[12][0]:
            angle_max.color = 'r'
            angle_max.hovercolor = 'r'
        else:
            angle_max.color = 'g'
            angle_max.hovercolor = 'g'
    # When any of the sliders are changed, update the figure
    sphi1.on_changed(update)
    sphi2.on_changed(update)
    spsi1.on_changed(update)
    spsi2.on_changed(update)
    sLambda1.on_changed(update)
    sLambda2.on_changed(update)
    sdho1.on_changed(update)
    sdho2.on_changed(update)
    sAR1.on_changed(update)
    sAR2.on_changed(update)
    sptc1.on_changed(update)
    sptc2.on_changed(update)
    # When a new input is submitted to the text boxes, update the figure
    Po1box.on_submit(update)
    To1box.on_submit(update)
    mdotbox.on_submit(update)
    Omegabox.on_submit(update)
    Wbox.on_submit(update)
    tbox.on_submit(update)
    gbox.on_submit(update)
    nbox.on_submit(update)
    a0box.on_submit(update)
    # Reset button to return sliders to initial values
    resetax = plt.axes([0.42, 0.01, 0.08, 0.04])
    reset_button = Button(resetax, 'Reset', color='1.0', hovercolor='0.5')

    def reset(event):
        """Reset the sliders"""
        Po1, To1, mdot, Omega, W, t, g, *_, ain = original_turbine[11]
        sphi1.reset()
        sphi2.reset()
        spsi1.reset()
        spsi2.reset()
        sLambda1.reset()
        sLambda2.reset()
        sdho1.reset()
        sdho2.reset()
        sAR1.reset()
        sAR2.reset()
        sptc1.reset()
        sptc2.reset()
        Po1box.set_val(Po1/10**5)
        To1box.set_val(To1)
        mdotbox.set_val(mdot)
        Omegabox.set_val(Omega*60/(2*np.pi))
        Wbox.set_val(W/10**6)
        tbox.set_val(t*10**3)
        gbox.set_val(g*10**3)
        nbox.set_val(n)
        a0box.set_val(ain)
        update(event)
    reset_button.on_clicked(reset)
    # Button that sets the exit loading equal to inlet to
    # create repeating stages
    repeatingax = plt.axes([0.42, 0.14, 0.08, 0.04])
    repeating_button = Button(repeatingax, 'Repeating',
                              color='1.0', hovercolor='0.5')

    def repeat(event):
        """Set output loading to input"""
        sphi2.set_val(sphi1.val)
        spsi2.set_val(spsi1.val)
        sLambda2.set_val(sLambda1.val)
        sdho2.set_val(sdho1.val)
        sAR2.set_val(sAR1.val)
        sptc2.set_val(sptc1.val)
        update(event)
    repeating_button.on_clicked(repeat)
    repeating_button.on_clicked(update)
    # Button to plot the annulus of the turbine
    annulusax = plt.axes([0.42, 0.0966, 0.08, 0.04])
    annulus_button = Button(annulusax, 'Annulus', color='1.0',
                            hovercolor='0.5')

    def annulus_plot(event):
        """Plot the output when clicked"""
        annulus(new_turbine)
    annulus_button.on_clicked(annulus_plot)
    # Button to find an optimum given the current variables
    optax = plt.axes([0.42, 0.0533, 0.08, 0.04])
    opt_button = Button(optax, 'Optimise', color='1.0', hovercolor='0.5')

    def find_opt(event):
        """Find an optimum arrangement"""
        global new_turbine
        phi, psi, Lambda, AR, dho = optimise(new_turbine)
        sphi1.set_val(phi[0])
        sphi2.set_val(phi[1])
        spsi1.set_val(psi[0])
        spsi2.set_val(psi[1])
        sLambda1.set_val(Lambda[0])
        sLambda2.set_val(Lambda[1])
        sAR1.set_val(AR[0])
        sAR2.set_val(AR[1])
        sdho1.set_val(dho[0])
        sdho2.set_val(dho[1])
        update(event)
    opt_button.on_clicked(find_opt)
    # These create dummy variables to ensure buttons
    # can be referenced outside of function
    repeatingax._button = repeating_button
    resetax._button = reset_button
    annulusax._button = annulus_button
    optax._button = opt_button


def annulus(turbine_data):
    """Plot annulus size through the turbine"""

    # Extract values form turbine function output
    chords = [[i[3], i[4]] for i in turbine_data[5]]
    H1s = [i[1] for i in turbine_data[5]]
    H2s = [i[2] for i in turbine_data[5]]
    H3s = [i[3] for i in turbine_data[5]]
    rm = [i[0] for i in turbine_data[5]]
    Ro = [i[8] for i in turbine_data[5]]
    Ri = [i[9] for i in turbine_data[5]]
    n = turbine_data[11][12]
    g = turbine_data[11][6]
    # Initialise lists
    length = 0
    x_hub = [0]
    x_cas = [0]
    r_hub = [rm[0]-H1s[0]/2]
    r_cas = [rm[0]+H1s[0]/2]
    # Loop over every stage
    _, ax = plt.subplots()
    for i in range(n):
        # Extract stage parameters
        Cxst, Cxro = chords[i]
        # Hub line including shroud cavity for stator
        hub1x = length+0.125*Cxst-g/2
        hub2x = hub1x
        hub3x = length+1.375*Cxst+g/2
        hub4x = hub3x
        hub5x = length+1.5*Cxst
        hub6x = length+1.5*Cxst+1.5*Cxro
        x_hub += [hub1x, hub2x, hub3x, hub4x, hub5x, hub6x]
        hub5y = rm[i]-H2s[i]/2
        hub6y = rm[i]-H3s[i]/2
        hub1y = r_hub[6*i]+(0.125*Cxst-g/2)/(1.5*Cxst)*(hub5y-r_hub[6*i])
        hub2y = hub1y-0.001-g
        hub4y = r_hub[6*i]+(1.375*Cxst+g/2)/(1.5*Cxst)*(hub5y-r_hub[6*i])
        hub3y = hub4y-0.001-g
        r_hub += [hub1y, hub2y, hub3y, hub4y, hub5y, hub6y]
        # Casing line including shroud cavity for rotor
        cas1x = length+1.5*Cxst
        cas2x = cas1x+0.125*Cxro-g/2
        cas3x = cas2x
        cas4x = cas1x+1.375*Cxro+g/2
        cas5x = cas4x
        cas6x = cas1x+1.5*Cxro
        x_cas += [cas1x, cas2x, cas3x, cas4x, cas5x, cas6x]
        cas1y = rm[i]+H2s[i]/2
        cas6y = rm[i]+H3s[i]/2
        cas2y = cas1y+(0.125*Cxro-g/2)/(1.5*Cxro)*(cas6y-cas1y)
        cas3y = cas2y+0.001+g
        cas5y = cas1y+(1.375*Cxro+g/2)/(1.5*Cxro)*(cas6y-cas1y)
        cas4y = cas5y+0.001+g
        r_cas += [cas1y, cas2y, cas3y, cas4y, cas5y, cas6y]
        # Calculate the length along the turbine
        length += Cxst*1.5+Cxro*1.5
        # Plot stator, starting at upstream root and working anti-clockwise
        y1 = r_cas[6*i]+0.25/1.5*(cas1y-r_cas[6*i])
        y2 = r_hub[6*i]+0.25/1.5*(hub5y-r_hub[6*i])
        y3 = r_hub[6*i]+(0.125*Cxst+g/2)/(1.5*Cxst)*(hub5y-r_hub[6*i])
        y4 = y3-0.001
        y6 = r_hub[6*i]+(1.375*Cxst-g/2)/(1.5*Cxst)*(hub5y-r_hub[6*i])
        y5 = y6-0.001
        y7 = r_hub[6*i]+1.25/1.5*(hub5y-r_hub[6*i])
        y8 = r_cas[6*i]+1.25/1.5*(cas1y-r_cas[6*i])
        x1 = x_cas[6*i]+0.25*Cxst
        x2 = x_cas[6*i]+0.125*Cxst+g/2
        x3 = x_cas[6*i]+1.375*Cxst-g/2
        x4 = x_cas[6*i]+1.25*Cxst
        x = [x1, x1, x2, x2, x3, x3, x4, x4]
        y = [y1, y2, y3, y4, y5, y6, y7, y8]
        plt.plot(x, y, color='0.5', lw=0.8)
        plt.plot([x1, x4], [y1, y7], color='0.5', lw=0.8)
        plt.plot([x1, x4], [y2, y8], color='0.5', lw=0.8)
        # Plot rotor, starting at upstream root and working clockwise
        y1 = hub5y+0.25/1.5*(hub6y-hub5y)
        y2 = cas1y+0.25/1.5*(cas6y-cas1y)
        y3 = cas1y+(0.125*Cxro+g/2)/(1.5*Cxro)*(cas6y-cas1y)
        y4 = y3+0.001
        y6 = cas1y+(1.375*Cxro-g/2)/(1.5*Cxro)*(cas6y-cas1y)
        y5 = y6+0.001
        y7 = cas1y+1.25/1.5*(cas6y-cas1y)
        y8 = hub5y+1.25/1.5*(hub6y-hub5y)
        x1 = cas1x+0.25*Cxro
        x2 = cas1x+0.125*Cxro+g/2
        x3 = cas1x+1.375*Cxro-g/2
        x4 = cas1x+1.25*Cxro
        x = [x1, x1, x2, x2, x3, x3, x4, x4]
        y = [y1, y2, y3, y4, y5, y6, y7, y8]
        plt.plot(x, y, color='0.2', lw=0.8)
        plt.plot([x1, x4], [y1, y7], color='0.2', lw=0.8)
        plt.plot([x1, x4], [y2, y8], color='0.2', lw=0.8)
    # Duplicate some entries for square stages
    Ri = [val-0.001-g for val in Ri for _ in range(7)]
    Ro = [val+0.001+g for val in Ro for _ in range(7)]
    for i in range(1, n):
        r_hub.insert(6*i+i-1, r_hub[6*i+i-1])
        r_cas.insert(6*i+i-1, r_cas[6*i+i-1])
        x_hub.insert(6*i+i-1, x_hub[6*i+i-1])
        x_cas.insert(6*i+i-1, x_cas[6*i+i-1])
    # Stages at constant thickness
    for i in range(0, n):
        dr_hub = r_hub[7*i+6]-r_hub[7*i]
        Ri[7*i:7*(i+1)] = list(np.interp(x_hub[7*i:7*(i+1)],
                                         [x_hub[7*i], x_hub[7*i+6]],
                                         [Ri[7*i], Ri[7*i]+dr_hub]))
    for i, j in enumerate(Ri):
        if j < 0:
            Ri[i] = 0
    for i in range(0, n):
        dr_cas = r_cas[7*i+6]-r_cas[7*i]
        Ro[7*i:7*(i+1)] = list(np.interp(x_cas[7*i:7*(i+1)],
                                         [x_cas[7*i], x_cas[7*i+6]],
                                         [Ro[7*i], Ro[7*i]+dr_cas]))
    # Plot lines separating stages
    plt.plot([0, 0], [r_hub[0], 0], color='black', lw=1.2)
    plt.plot([0, 0], [r_cas[0], Ro[0]], color='black', lw=1.2)
    for i in range(1, n):
        min_Ri = min(Ri[7*i-1], Ri[7*i])
        max_Ro = max(Ro[7*i-1], Ro[7*i])
        plt.plot([x_hub[7*i-1], x_hub[7*i-1]],
                 [r_hub[7*i-1], min_Ri], color='black', lw=1.2)
        plt.plot([x_cas[7*i-1], x_cas[7*i-1]],
                 [r_cas[7*i-1], max_Ro], color='black', lw=1.2)
    plt.plot([x_hub[-1], x_hub[-1]], [r_hub[-1], 0], color='black', lw=1.2)
    plt.plot([x_cas[-1], x_cas[-1]], [r_cas[-1], Ro[-1]],
             color='black', lw=1.2)
    # Plot the results
    plt.plot(x_hub, r_hub, label='Hub', color=(0.0, 0.6, 0.9), linewidth=0.2)
    plt.plot(x_cas, r_cas, label='Casing', color=(0.9, 0.2, 0.0),
             linewidth=0.2)
    plt.plot(x_hub, Ri, color=(0.0, 0.6, 0.9), linewidth=0.2)
    plt.plot(x_cas, Ro, color=(0.9, 0.2, 0.0), linewidth=0.2)
    plt.fill_between(x_hub, Ri, r_hub, color=(0.0, 0.6, 0.9))
    plt.fill_between(x_cas, r_cas, Ro, color=(0.9, 0.2, 0.0))
    plt.xlabel('Length along turbine (m)')
    plt.ylabel('Radius (m)')
    plt.axis('equal')
    ax.set_ybound(0, 1.05*np.amax(Ro))
    ax.set_xbound(0, x[-1])
    # Zoom button
    zoomax = plt.axes([0.4, 0.9, 0.08, 0.04])
    zoom_button = Button(zoomax, 'Zoom', color='1.0', hovercolor='0.5')
    resetax = plt.axes([0.5, 0.9, 0.08, 0.04])
    reset_button = Button(resetax, 'Reset', color='1.0', hovercolor='0.5')

    def zoom(val):
        """Zoom in on annulus"""
        ax.set_ybound(0.95*np.amin(r_hub), 1.05*np.amax(r_cas))
        ax.set_xbound(-0.05*x[-1], 1.05*x[-1])

    def reset(val):
        """Reset view"""
        ax.axis('equal')
        ax.set_ybound(0, 1.05*np.amax(Ro))
        ax.set_xbound(0, x[-1])
    zoom_button.on_clicked(zoom)
    reset_button.on_clicked(reset)
    zoomax._button = zoom_button
    resetax._button = reset_button
    plt.show()
