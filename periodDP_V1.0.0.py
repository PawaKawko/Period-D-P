#The main idea here that we try to approximate the light curve by Fourier series with different periods
#and choose that one, for which the sum of square deviations dots from the approximation is the smallest.
#Then programm build a light curve and phase curve. All dots that are stands out from the approximation
#is cutted off. Program writes in the file the pictures of phase curves and data with cutted points

Version = "V1.0.0"

"""==========================================="""
"""IMPORTING LIBRUARIES"""
"""==========================================="""

import scipy.optimize as spo     #for the method of LS
import numpy as np               #for math stuff
import matplotlib.pyplot as plt  #for plotting
import time                      #to know time of calculations
import tkinter as tnk            #graphic interface
import os                        #to work with directories
import decimal
import matplotlib.font_manager
import warnings
warnings.filterwarnings("ignore")

"""==========================================="""
"""Path to files"""
"""==========================================="""

path_file = os.getcwd()          #constant for the path to the folder, where code is stored

"""==========================================="""
"""ERRORS"""
"""==========================================="""

def Error_1():                   #function to display an error in Manual mode that is caused by inputting not correct value of T
    window_error = tnk.Tk()
    bcg_cl = '#ffff00'
    window_error.title("Period D&P " + Version)
    w = 550
    h = 180
    window_error.geometry(str(w) + 'x' + str(h))
    window_error.config(bg=bcg_cl)
    window_error.resizable(width=False, height=False)

    lb_error = tnk.Label(window_error, font = ('Algerian', 19), text = 'Error #1', bg=bcg_cl)
    lb_describtion_1 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = 'The program has not found minimum in periodogram', bg=bcg_cl)
    lb_describtion_2 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = 'Please try another period or its error', bg=bcg_cl)

    lb_error.place(x = 200, y = 30)                #their place on the window
    lb_describtion_1.place(x = 20, y = 80)
    lb_describtion_2.place(x = 90, y = 110)
    window_error.mainloop()

def Error_2(File, Number_error):                   #function to display an error that arrive due to absense of some files
    window_error = tnk.Tk()
    bcg_cl = '#9999FF'
    window_error.title("Period D&P " + Version)
    w = 850
    h = 180
    window_error.geometry(str(w) + 'x' + str(h))
    window_error.config(bg=bcg_cl)
    window_error.resizable(width=False, height=False)

    if Number_error == 1:
        error_text = 'The program has not found ' + File
        lb_error = tnk.Label(window_error, font = ('Algerian', 24), text = 'Error #2.1', bg=bcg_cl)
    if Number_error == 2:
        error_text = 'Problem while reading ' + File
        lb_error = tnk.Label(window_error, font = ('Algerian', 24), text = 'Error #2.2', bg=bcg_cl)
    lb_describtion_1 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = error_text, bg=bcg_cl)
    lb_describtion_2 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = 'Please check and repeat', bg=bcg_cl)

    lb_error.place(x = 350, y = 30)
    lb_describtion_1.place(x = 20, y = 80)
    lb_describtion_2.place(x = 240, y = 110)
    window_error.mainloop()

"""==========================================="""
"""TRIGONOMETRIC POLYNOMIAL FUNCTIONS"""
"""==========================================="""

def sin(t, pp, n):                  #approximation of function by Fourie series (t -> x_data, pp - parameters)
    x = np.zeros(len(t))
    x += pp[0]
    for i in range(n):
        x += pp[2*i+2]*np.sin(2*np.pi*t*(i+1)/pp[1]+pp[2*i+3])      # x = SUM( A*sin(2*pi*n*t/T + phi))
    return x

def sin1(t, pp, n):                 #the same as sin(), but give you not array, but a value
    y = pp[0]
    for i in range(n):
        y += pp[2*i+2]*np.sin(2*np.pi*t/pp[1]*(i+1)+pp[2*i+3])
    return y

def Trend(t, pp):
    y = pp[0] + pp[1] * t
    return y

def Polymom(t, pp, n):
    y = pp[0]
    for i in range(1, n+1):
        y+= pp[i]*(t**i)
    return y

"""==========================================="""
"""READING DATA FROM FILE"""
"""==========================================="""

def read_data(name):                             #function to read raw data
    Name = path_file + '/data/' + name           #data is stored in the same sirectory in the folder "data"
    try:
        Data = np.loadtxt(Name)
        x = np.array(Data[:,0])
        y = np.array(Data[:,1])
        y_err = np.array(Data[:,2])
        Error_program = 0
    except FileNotFoundError:
        Error_program = 1
        x = 0
        y = 0
        y_err = 0
    except ValueError:
        Error_program = 2
        x = 0
        y = 0
        y_err = 0
    return x, y, y_err, Error_program

"""==========================================="""
"""READING PARAMETERS AND TASKS FROM FILE"""
"""==========================================="""
def read_parametrs(Parametrs_file):                     #function to read parameters for work
    try:
        parametrs = np.loadtxt(Parametrs_file)
        n_app_T = int(parametrs[0])                    #number of additions in Fourie series in function Approximation T
        n_approximation = int(parametrs[1])            #number of additions in Fourie series in function becoming perfect
        edge_appr_T = float(parametrs[2])              #to cut minimum in periodogram
        TT_min_par = float(parametrs[3])               #the minimum value of period in Periodogram
        Presize_appr_T = float(parametrs[4])           #the distance between points in the Periodogram
        ratio = float(parametrs[5])                    #size of Phase picture (x:y)
        dpi_picture = int(parametrs[6])               #quality of picture
        dots_size = int(parametrs[7])                  #size of dots ob phase curves
        Start_phase = float(parametrs[8])              #start phase of observation
        Error_program = 0
        return n_app_T, n_approximation, edge_appr_T, TT_min_par, Presize_appr_T, ratio, dpi_picture, dots_size, Start_phase, Error_program
    except FileNotFoundError:
        Error_program = 1
        return 0,0,0,0,0,0,0,0,0,Error_program
    except ValueError:
        Error_program = 2
        return 0,0,0,0,0,0,0,0,0,Error_program

def read_task(task_file):
    try:
        Task = np.genfromtxt(task_file, dtype='str')
        for value in Task:
            if not len(value.split('.')) == 2:
                raise ValueError
        Error_program_task = 0
    except FileNotFoundError:
        Error_program_task = 1
    except ValueError:
        Error_program_task = 2
    return Task, Error_program_task
"""==========================================="""
"""CALCULATING PRESIZE VALUE OF PERIOD"""
"""==========================================="""

def first_approximation(Tappr, A0, x, y, y_err, n_approximation, name, n_app_T, ans_start, dpi_picture, dots_size, ratio, I):
    p0 = np.ones(2*n_approximation + 2)             #start conditions
    p0[0] = ans_start[0]                                #first = ideal from periodogram
    p0[1] = Tappr

    if(n_approximation > n_app_T):                   #set conditions the same as the best in ApproximationT
        for i in range(2*n_app_T):
             p0[i+2] = ans_start[i+1]
    else:
        for i in range(2*n_approximation + 2):
             p0[i+2] = ans_start[i]

    fun = lambda pp: (y - sin(x, pp, n_approximation))/y_err       #core of least squares
    ans = spo.leastsq(fun, p0, full_output=1)
    sigma = np.sum((y - sin(x, ans[0], n_approximation))**2)/len(x)
    error = np.sqrt(np.diag(ans[1]*sigma))

    T_ideal = ans[0][1]
    error_T = error[1]
    ans_ideal  = ans[0]         #ideal parametrs

    order_Error = -int(np.log10(error_T))+1                  #evaluate order of Error
    save_path = path_file + '/Results/' + name + '/'         #save results in the folder "Results"

    fig = plt.figure(2 + I * 6)             #plot dots and curve
    plt.gca().invert_yaxis()                        #to invert y axis
    fig.set_size_inches(20, 7)
    plt.rc('xtick', labelsize=20)                   #size of tics
    plt.rc('ytick', labelsize=20)
    plt.plot(x, y, '.b')                            #blue dots
    plt.xlabel('BJD', fontsize = 20)                #name of axis
    plt.ylabel('$\Delta$T, mmag', fontsize = 20)
    plt.title('Light curve', fontsize = 20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.savefig(save_path + name + " light curve.png", dpi = 300)                       #without approximation
    xx = np.linspace(min(x), max(x), len(x))                                            #to plot approximation on the parts, where are not data
    plt.plot(xx, sin(xx, ans_ideal, n_approximation), '-r')
    plt.savefig(save_path + name + " light curve with approximation.png", dpi = dpi_picture)    #with approximation
    plt.close()

    return ans_ideal, np.round(T_ideal, order_Error)

def remove_trends(x, y, y_err, ans_ideal, name, n_approximation, dpi_picture, dots_size, ratio, I):
    y_new = y.copy()

    sigma = np.sqrt(np.sum((y - sin(x, ans_ideal, n_approximation))**2)/len(x))
    key = True
    for index in range(len(x)):
        Condition = np.abs(y[index] - sin1(x[index], ans_ideal, n_approximation)) > (3*sigma)
        if key and Condition:
            Index1 = index
            key = False
        if (not key) and (not Condition):
            Index2 = index
            key = True

            if (Index2 - Index1) > 2:               #removing trend
                y_trend = y[Index1:(Index2+1)]
                y_err_trend = y_err[Index1:(Index2+1)]
                x_trend = x[Index1:(Index2+1)]
                trend = y_trend - sin(x_trend, ans_ideal, n_approximation)

                p0 = [1, 1]
                fun = lambda pp: (trend - Trend(x_trend, pp))/y_err_trend
                ans = spo.leastsq(fun, p0, full_output=1)

                y_new[Index1:(Index2+1)] -= Trend(x_trend, ans[0])

    save_path = path_file + '/Results/' + name + '/'         #save results in the folder "Results"
    fig = plt.figure(3 + I*6)             #plot dots and curve
    plt.gca().invert_yaxis()
    fig.set_size_inches(20, 7)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.plot(x, y, '.g')
    plt.plot(x, y_new, '.b')
    xx = np.linspace(min(x), max(x), len(x))                                            #to plot approximation on the parts, where are not data
    plt.plot(xx, sin(xx, ans_ideal, n_approximation), '-r')
    plt.xlabel('BJD', fontsize = 20)
    plt.ylabel('$\Delta$T, mmag', fontsize = 20)
    plt.title('Light curve (trends)', fontsize = 20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.savefig(save_path + name + " light curve no trends.png", dpi = 300)
                       #without approximation
    return y

def remove_linear(x, y, y_err):
    number = 10
    key = 0
    for i in range(10):
        key += np.sign(y[1] - y[0])
    key = np.round(key/10)

    key1 = 0
    for i in range(1, len(x)):
        if not np.sign(y[i] - y[i-1]) == key:
            key1 += 1
        if key1 == 150:
            break
    if i > number:
        x_new = x[i:]
        y_new = y[i:]
        y_err_new = y_err[i:]
        return x_new, y_new, y_err_new

    else:
        return x, y, y_err

def remove_trends_2(x, y, y_err, ans_ideal, name, ftype, n_approximation, dpi_picture, dots_size, ratio, I):
    n = 3

    start = []                              #cutting in parts
    end = []
    start.append(0)
    delta = x[1] - x[0]
    for i in range(len(x)-1):
        if (x[i+1] - x[i]) > 100*delta:
            end.append(i)
            start.append(i+1)
    end.append(len(x)-1)

    save_path = path_file + '/Results/' + name + '/'
    fig, axs = plt.subplots(4, 1)
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(30, 30)
    plt.rc('ytick', labelsize=30)
    axs[0].set_title('Light curve (trends) - ' + name, fontsize = 35)
    xx = np.linspace(np.min(x), np.max(x), len(x))
    axs[0].plot(x, y, '.g')
    #axs[0].plot(xx, sin(xx, ans_ideal, n_approximation), '.r')
    plt.rc('xtick', labelsize=30)
    for i in range(4):
        axs[i].set_ylabel('$\Delta$T, mmag', fontsize = 30)
        axs[i].invert_yaxis()

    X_new = np.array([])
    Y_new = np.array([])
    Y_err_new = np.array([])

    for i in range(len(start)):
        x_part = x[start[i]:end[i]].copy()
        y_part = y[start[i]:end[i]].copy()
        y_err_part = y_err[start[i]:end[i]].copy()

        x_part, y_part, y_err_part = remove_linear(x_part, y_part, y_err_part)   # ?????????????

        if len(x_part) > n+1:
            p0 = 0.1 * np.ones(n+1)
            fun = lambda pp: (y_part - sin(x_part, ans_ideal, n_approximation) - Polymom(x_part, pp, n)) / y_err_part
            ans = spo.leastsq(fun, p0, full_output=1)

            xx = np.linspace(np.min(x_part), np.max(x_part), len(x_part))
            axs[1].plot(x_part, y_part - sin(x_part, ans_ideal, n_approximation), '.g')
            axs[1].plot(xx, Polymom(xx, ans[0], n), '.r')

            y_part -= Polymom(x_part, ans[0], n)
            axs[2].plot(x_part, y_part - sin(x_part, ans_ideal, n_approximation), '.g')
        else:
            axs[1].plot(x_part, y_part - sin(x_part, ans_ideal, n_approximation), '.g')
            axs[2].plot(x_part, y_part - sin(x_part, ans_ideal, n_approximation), '.g')

        X_new = np.concatenate((X_new, x_part))
        Y_new = np.concatenate((Y_new, y_part))
        Y_err_new = np.concatenate((Y_err_new, y_err_part))

    x = X_new.copy()
    y = Y_new.copy()
    y_err = Y_err_new.copy()

    sigma = np.sqrt(np.sum((y - sin(x, ans_ideal, n_approximation))**2) / len(x) )
    axs[2].axhline(y = 3*sigma)
    axs[2].axhline(y = -3*sigma)

    Condition = abs(y - sin1(x, ans_ideal, n_approximation)) < 3*sigma
    x, y, y_err = x[Condition], y[Condition], y_err[Condition]

    p0 = ans_ideal
    fun = lambda pp: (y - sin(x, pp, n_approximation))/y_err
    ans = spo.leastsq(fun, p0, full_output=1)
    sigma = np.sum((y - sin(x, ans[0], n_approximation))**2)/len(x)
    error = np.sqrt(np.diag(ans[1]*sigma))

    order_Error = -int(np.log10(error[1]))+1                  # evaluate order of Error

    Mean = np.mean(y)
    SS_res = np.sum((y -  sin(x, ans[0], n_approximation))**2)
    SS_tot = np.sum((y - Mean)**2)
    R_2 = 1 - SS_res/SS_tot

    chi_2 = np.sum(((y -  sin(x, ans[0], n_approximation))**2)/y_err**2)/( len(x) - (2*n_approximation + 1))

    def sin_chi(t):
        pp = ans[0]
        z = np.zeros(len(x))
        z += pp[0]
        for i in range(n_approximation):
            z += pp[2*i+2] * np.sin(2*np.pi*x*(i+1)/t + pp[2*i+3])
        chi_2_new = np.sum(((y -  z)**2)/y_err**2)/( len(x) - (2*n_approximation + 1))
        return (chi_2_new - chi_2 - 1)

    root = spo.fsolve(sin_chi, ans[0][1])

    xx = np.linspace(np.min(x), np.max(x), len(x))
    #axs[3].plot(xx, sin(xx, ans[0], n_approximation), '.r')
    axs[3].plot(x, y, '.g')

    plt.xlabel('BJD', fontsize = 20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.savefig(save_path + name + " light curve trends.png", dpi = 300)

    NName = name + "_detrended." + ftype                      #save data in the same file type
    completeName = os.path.join(save_path, NName)
    with open(completeName, 'w+') as f:
        for i in range(len(x)):
            f.write(str(x[i]) + ' ' + str(y[i]) + ' ' + str(y_err[i]) + '\n')

    return x, y, y_err, np.round(ans[0][1], order_Error), np.round(error[1], order_Error), ans[0][1]-root[0], R_2, chi_2, ans[0]

def phase_curve(T_ideal, answ, x, y, y_err, n_approximation, name, ftype, ratio, dpi_picture, dots_size, Start_phase, key_number, I):

    d = decimal.Decimal(str(T_ideal))
    if key_number == 1:
        order_Error = -d.as_tuple().exponent
    else:
        order_Error = -d.as_tuple().exponent-1
    Number_periods = (x - x[0])/T_ideal                         #To build phase curve
    Number_periods = Number_periods.astype(int)
    I_max = np.argmax(y)
    X_E = (x - x[0])/T_ideal - Number_periods
    X_E -= X_E[I_max]
    X_E[X_E < 0] += 1
    save_path = path_file + '/Results/' + name + '/'

    B = max(y) - min(y)
    hfont = {'fontname':'Helvetica'}
    fig = plt.figure(4 + I * 6)
    plt.gca().invert_yaxis()
    fig.set_size_inches(ratio*7, 7)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    strin = 'Phase (P = ' + str(np.round(Start_phase + x[I_max], order_Error)) + ' +' + str(np.round(T_ideal, order_Error)) + '*E)'
    plt.xlabel(strin, fontsize = 20, **hfont)
    plt.ylabel('$\Delta$T, mmag', fontsize = 20, **hfont)
    plt.plot(X_E, y,  color = 'green', linestyle = '', marker = '.', markersize = dots_size)
    plt.text(0, (np.min(y) + 1/30*B), name, fontsize = 20, **hfont)
    if key_number == 1:
        plt.savefig(save_path + name + "phase curve first.png", dpi = dpi_picture)
    else:
        plt.savefig(save_path + name + "phase curve.png", dpi = dpi_picture)
    plt.close()

    NName = name + " phase curve." + ftype                      #save data in the same file type
    completeName = os.path.join(save_path, NName)
    with open(completeName, 'w+') as f:
        for i in range(len(x)):
            f.write(str(X_E[i]) + ' ' + str(y[i]) + ' ' + str(y_err[i]) + '\n')

"""==========================================="""
"""COMPUTING APPROXIMATE VALUE OF PERIOD"""
"""==========================================="""

def Approximation_T(x, y, y_err, A, n_app_T, edge_appr_T, T_max, T_min, Presize_appr_T, name, dpi_picture, I):

    N_N = int(T_max/Presize_appr_T)         #number of dots in this area
    X_min = 0                               #just for fun(do not change)

    def sin2(t, T_Tt, pp, nn):              #approximation of function that take x data, period and parametrs and give the approximation function
        x = np.zeros(len(t))                #make array x lenth of x-data and full zero
        x += pp[0]
        for i in range(nn):                 #additions in Fourie series
            x += pp[2*i + 1]*np.sin(2*np.pi*t/T_Tt*(i+1)+pp[2*i+2])
        return x                            #return tha value of approximation function

    def sigma(xx, yy, yy_err, T_Tt, p00, nn):                              #function to find the sum of squares for each T
        fun = lambda pp: (yy - sin2(xx, T_Tt, pp, nn))/yy_err              #core of least squares
        ans = spo.leastsq(fun, p00, full_output=1)
        Sigma = np.sum((yy-sin2(xx, T_Tt, ans[0], nn))**2)/(len(x)*(len(x)-1))          #ans[0] - parametrs: amplitudes and phases
        return Sigma, ans[0]

    p0 = np.ones(2*n_app_T+1)
    p0[0], p0[1] = 0, A                   #main amplitude

    x_sigma = np.linspace(T_min, T_max, N_N)
    y_sigma = np.zeros(N_N)

    for i in range(len(x_sigma)):                                                 #for each dot
        if(x_sigma[i] == T_min):
            y_sigma[i], PP0 = sigma(x, y, y_err, x_sigma[i], p0, n_app_T)           #find y and ideal parametrs
        else:
            y_sigma[i], PP0 = sigma(x, y, y_err, x_sigma[i], PP0, n_app_T)          #start condition = ideal for previous

    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    fig = plt.figure(1 + I * 6)
    fig.set_size_inches(20, 6)
    save_path = path_file + '/Results/' + name + '/'
    plt.xlabel('Period', fontsize = 20)
    plt.ylabel('Sigma', fontsize = 20)
    plt.plot(x_sigma, y_sigma, color = '#FF0000', ls = '-', lw = 2)
    plt.savefig(save_path + name + "periodogram.png", dpi = dpi_picture)
    plt.close()

    value_error = False
    if ((np.min(y_sigma)/np.max(y_sigma)) < 0.3):
        value_error = True                                                  #there is no true minimum

    if value_error:
        Index = np.argmin(y_sigma)
        X_min = x_sigma[Index]
        PP0 = sigma(x, y, y_err, X_min, p0, n_app_T)[1]
        Error_program = 0

    else:
        X_min = 0
        order_ld = 0
        local_delta = 0
        PP0 = 0
        Error_program = 1

    return X_min, PP0, Error_program

"""==========================================="""
"""CREATING WINDOW AND GENERAL WIDJETS"""               #graphic interface
"""==========================================="""

def Automatic_work():

    window = tnk.Tk()
    bcg_cl = '#9999FF'
    window.title("Period D&P " + Version)
    w = 550
    h = 320
    window.geometry(str(w) + 'x' + str(h))
    window.config(bg=bcg_cl)
    window.resizable(width=False, height=False)

    lb_head = tnk.Label(window, font = ('Bookman Old Style', 18), text = 'Task file:', bg=bcg_cl)
    lb_head.place(x = 20, y = 30)
    lb_head = tnk.Label(window, font = ('Bookman Old Style', 15), text = 'Upper evaluation of T:', bg=bcg_cl)
    lb_head.place(x = 20, y = 100)
    lb_Par_file = tnk.Label(window, font = ('Bookman Old Style', 15), text = 'Name of file with parametrs:', bg=bcg_cl)
    lb_Par_file.place(x = 20, y = 170)

    ent_TaskFile = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_TaskFile.place(x = 110, y = 70)
    ent_Tmax = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_Tmax.place(x = 110, y = 135)
    ent_Par_file = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_Par_file.place(x = 110, y = 200)
    ent_Par_file.insert(0, 'Parametrs.txt')
    ent_TaskFile.insert(0, 'Task.txt')

    time_text = ['Time of calculations', 'min', 's']
    progress = ['Progress', 'from']
    lbprogress = [tnk.Label(window, font = ('Century', 14), text = progress[i], bg=bcg_cl) for i in range(2)]
    lbtime = [tnk.Label(window, font = ('Century', 14), text = time_text[i], bg=bcg_cl) for i in range(3)]
    enttime = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 4) for i in range(2)]
    ent_progress_1 = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 3)
    ent_progress_2 = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 3)
    lbprogress[0].place(x = 340, y = 50)
    lbprogress[1].place(x = 422, y = 80)
    ent_progress_1.place(x = 380, y = 80)
    ent_progress_2.place(x = 470, y = 80)
    lbtime[0].place(x = 340, y = 130)
    lbtime[1].place(x = 427, y = 160)
    lbtime[2].place(x = 530, y = 160)
    enttime[0].place(x = 370, y = 165)
    enttime[1].place(x = 475, y = 165)

    """==========================================="""
    """MAIN FUNCTION FOR AUTOMATIC MODE"""
    """==========================================="""

    def automatic_regime():
        start_time_0 = time.time()

        enttime[0].delete(0, len(enttime[0].get()))
        enttime[1].delete(0, len(enttime[1].get()))
        ent_progress_1.delete(0, len(ent_progress_1.get()))
        ent_progress_2.delete(0, len(ent_progress_2.get()))

        if (not os.path.exists('Results')):      # Create target Directory
            os.mkdir('Results')

        task_file = ent_TaskFile.get()
        TT_max = float(ent_Tmax.get())
        Parametrs_file = ent_Par_file.get()

        n_app_T, n_approximation, edge_appr_T, TT_min_par, Presize_appr_T, ratio, dpi_picture, dots_size, Start_phase, Error_program_par = read_parametrs(Parametrs_file)
        Task, Error_program_task = read_task(task_file)
        if Error_program_par:
            res = '    Error #2 -- Parametrs file\n'
            Error_2(Parametrs_file, Error_program_par)
        elif Error_program_task:
            res = '    Error #2 -- Task file\n'
            Error_2(task_file, Error_program_par)
        else:
            res = '  Name         T                   R^2                chi^2\n'
            N_stars = len(Task)
            ent_progress_1.insert(0, '0')
            ent_progress_2.insert(0, str(N_stars))

            for i in range(N_stars):
                try:
                    file = Task[i]
                    line = file.split('.')
                    name = line[0]
                    ftype = line[1]
                    res += name + '  '

                    sub_name = path_file + '/Results/' + name
                    if (not os.path.exists(sub_name)):
                        os.mkdir(sub_name)

                    x, y, y_err, Error_program_data = read_data(file)
                    if Error_program_data:
                        res += '     Error #2 -- data file\n'
                        Error_2(file, Error_program_data)
                    else:
                        print(name)
                        A0 = (np.max(y)-np.min(y)) / 2
                        Tappr, ans_start, Error_program_app_T = Approximation_T(x, y, y_err, A0, n_app_T, edge_appr_T, TT_max, TT_min_par, Presize_appr_T, name, dpi_picture, i)
                        if Error_program_app_T:
                            res += '   Error #1\n'
                            Error_1()
                        else:
                            ans_ideal, T = first_approximation(Tappr, A0, x, y, y_err, n_approximation, name, n_app_T, ans_start, dpi_picture, dots_size, ratio, i)
                            phase_curve(T,ans_ideal, x, y, y_err, n_approximation, name, ftype, ratio, dpi_picture, dots_size, Start_phase, 1, i)
                            x, y, y_err,  T, err_T, root, R_2, chi_2, ans_ideal = remove_trends_2(x, y, y_err, ans_ideal, name, ftype, n_approximation, dpi_picture, dots_size, ratio, i)
                            phase_curve(T,ans_ideal, x, y, y_err, n_approximation, name, ftype, ratio, dpi_picture, dots_size, Start_phase, 2, i)

                            # res += str(T) + ' +- ' + str(err_T)+ '  ' + str(root) + '   ' + str(R_2)+ '  ' + str(chi_2) + '\n'
                            res += str(T) + ' +- ' + str(err_T) + '   ' + str(R_2)+ '  ' + str(chi_2) + '\n'

                            k = int(ent_progress_1.get())
                            ent_progress_1.delete(0, len(ent_progress_1.get()))
                            ent_progress_1.insert(0, str(k+1))

                except:
                    print("Problem with " + str(i+1) + " star (" + name + "). Please check in manual mode")
                    res += 'Problem. Please check manually'
                    res += '\n'

            task_file = task_file.split('.')[0]

            results_path =  path_file + '/Results/' + 'results_' + task_file + '.dat'
            with open(results_path, 'w') as f:
                f.writelines(res)

        t_0 = time.time() - start_time_0
        enttime[1].insert(0, str(round(t_0)-60*int(t_0/60)))
        enttime[0].insert(0, str(int(t_0/60)))

    btn = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Calculate', bg = "blue", fg = "white", height = 1, width = 13, command = automatic_regime)
    btn.place(x = 50, y = 255)
    window.mainloop()

Automatic_work()
