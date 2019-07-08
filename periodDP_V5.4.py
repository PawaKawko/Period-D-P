#The main idea here that we try to approximate the light curve by Fourier series with different periods 
#and choose that one, for which the sum of square deviations dots from the approximation is the smallest.
#Then programm build a light curve and phase curve. All dots that are stands out from the approximation 
#is cutted off. Program writes in the file the pictures of phase curves and data with cutted points

#try
"""==========================================="""
"""IMPORTING LIBRUARIES"""
"""==========================================="""

import scipy.optimize as spo     #for the method of LS
import numpy as np               #for math stuff
import matplotlib.pyplot as plt  #for plotting
import time                      #to know time of calculations
import tkinter as tnk            #graphic interface
import os                        #to work with directories

"""==========================================="""
"""Path to files"""
"""==========================================="""

path_file = os.getcwd()

"""==========================================="""
"""Errors"""
"""==========================================="""

def Error_1():
    window_error = tnk.Tk()
    bcg_cl = '#ffff00'
    window_error.title("Period D&P V5.2")
    w = 550
    h = 180
    window_error.geometry(str(w) + 'x' + str(h))
    window_error.config(bg=bcg_cl)
    window_error.resizable(width=False, height=False)
    
    lb_error = tnk.Label(window_error, font = ('Algerian', 19), text = 'Error #1', bg=bcg_cl)                             
    lb_describtion_1 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = 'The program has not found minimum in periodogram', bg=bcg_cl)                 #words and labels
    lb_describtion_2 = tnk.Label(window_error, font = ('Bookman Old Style', 14), text = 'Please try another period or its error', bg=bcg_cl)
    
    lb_error.place(x = 200, y = 30)               #their place on the window
    lb_describtion_1.place(x = 20, y = 80)       
    lb_describtion_2.place(x = 90, y = 110) 
    window_error.mainloop()

"""==========================================="""
"""TRIGONOMETRIC POLYNOMIAL FUNCTIONS"""
"""==========================================="""

def sin(t, pp, n):                  #approximation of function by Fourie series (t -> x_data, pp - parameters)
    x = np.zeros(len(t))            #creating array with the size of data
    x += pp[0]                      #constant
    for i in range(n):
        x += pp[2*i+2]*np.sin(2*np.pi*t*(i+1)/pp[1]+pp[2*i+3])      # x = SUM( A*sin(t + φ))
    return x

def sin1(t, pp, n):                 #the same as sin but give you not array, but a value
    y = pp[0]
    for i in range(n):
        y += pp[0] + pp[2*i+2]*np.sin(2*np.pi*t/pp[1]*(i+1)+pp[2*i+3])
    return y

"""==========================================="""
"""READING DATA FROM FILE"""
"""==========================================="""

def read_data(name, ftype):
    Name = path_file + '\\data\\' + name + ftype    #data is stored in the same sirectory in the folder "data"
    with open(Name, 'r') as file:                   #each file should be named "Name_star.type"
        x, y, y_err = [], [], []                    #set arrays
        lines = file.readlines()                    #lines - list; read() - not work
        for i in lines: 
            if (not i.startswith("#")):             #to avoid comments
                data = i.split()                    #split into words because of spaces
                x.append(float(data[0]))
                y.append(float(data[1]))
                y_err.append(float(data[2]))
        x, y, y_err = np.array(x), np.array(y), np.array(y_err)      #to make arrays more cool and suitable with method of LS
    Number_of_elements_0 = len(x)  
    return x, y, y_err, Number_of_elements_0

"""==========================================="""
"""READING PARAMENTRS FROM FILE"""
"""==========================================="""
def read_parametrs(Parametrs_file):
    parametrs = []
    with open(Parametrs_file, 'r') as file:  
        for line in file:           
            if (not line.startswith("#")):         #to avoid comments
                parametrs.append(float(line)) 
    n_app_T = int(parametrs[0])                    #number of additions in Fourie series in function Approximation T
    n_becoming_perfect = int(parametrs[1])         #number of additions in Fourie series in function becoming perfect
    n_bec_per_sec = int(parametrs[2])              #number of additions in Fourie series in function becoming perfect second
    edge_appr_T = float(parametrs[3])              #to cut minimum
    Parametr_sigma = float(parametrs[4])           #to cut phase diagramm
    TT_min_par = float(parametrs[5])               #the minimum value of period in Periodogram
    Presize_appr_T = float(parametrs[6])           #the distance between points in the Periodogram
    ratio = float(parametrs[7])                    #size of Phase picture (x:y)
    max_width = float(parametrs[8])                #to cut Phase diagramm
    N_cutting = int(parametrs[9])                  #Number of pictures of Phase diagram and also detalization of cutting it
    N_fragmentation = int(parametrs[10])           #detalization of cuttind phase diagram
    dpi_picture = int (parametrs[11])              #quality of picture
    dots_size = int(parametrs[12])                 #size of dots ob phase curves
    return n_app_T, n_becoming_perfect, edge_appr_T, Parametr_sigma, TT_min_par, Presize_appr_T, ratio, N_cutting, n_bec_per_sec, max_width, N_fragmentation, dpi_picture, dots_size
"""==========================================="""
"""CALCULATING PRESIZE VALUE OF PERIOD"""
"""==========================================="""
def becoming_perfect(Tappr, A0, x, y, y_err, n_becoming_perfect, name, n_app_T, ans_start, dpi_picture, dots_size, I = 0, Repeats = 0):
    p0 = np.zeros(2*n_becoming_perfect + 2)             #start conditions 
    p0[0] = ans_start[0]                                #first = ideal from periodogram
    p0[1] = Tappr
    if(n_becoming_perfect > n_app_T):
        for i in range(2*n_app_T):
             p0[i+2] = ans_start[i+1]      
        for i in range(2*n_app_T + 2, 2*n_becoming_perfect + 2):
             p0[i] = 1   
    else:
        for i in range(2*n_becoming_perfect + 2):
             p0[i] = ans_start[i]     
             
    fun = lambda pp: (y - sin(x, pp, n_becoming_perfect))/y_err       #core of least squares
    ans = spo.leastsq(fun, p0, full_output=1)
    sigma = np.sum((y - sin(x, ans[0], n_becoming_perfect))**2)/len(x)
    error = np.sqrt(np.diag(ans[1]*sigma))
    
    T_ideal = ans[0][1]
    error_T = error[1]
    ans_ideal  = ans[0]         #ideal parametrs
                             
    order_Error = -int(np.log10(error_T))+1   #evaluate order of Error
    save_path = path_file + '\\Results\\' + name + '\\'         #save results in the folder "Results"

    fig = plt.figure(I*(Repeats+2) + 2)             #plot dots and curve
    plt.gca().invert_yaxis()                        #to invert y axis
    fig.set_size_inches(15, 6)                      #size
    plt.rc('xtick', labelsize=20)                   #size of tics
    plt.rc('ytick', labelsize=20) 
    plt.plot(x, y, '.b')                            #blue dots
    plt.xlabel('BJD', fontsize = 20)                #name of axis
    plt.ylabel('V mmag', fontsize = 20)
    plt.title('Light curve', fontsize = 20)
    plt.savefig(save_path + name + " light curve.png", dpi = 300)                       #without approximation
    xx = np.linspace(min(x), max(x), len(x))                                            #to plot approximation on the parts, where are not data
    plt.plot(xx, sin(xx, ans_ideal, n_becoming_perfect), '-r')
    plt.savefig(save_path + name + " light curve with approximation.png", dpi = dpi_picture)    #with approximation
    plt.show()                                                                          #to show plot during process

    return ans_ideal, np.round(T_ideal, order_Error), np.round(error_T, order_Error)

def becoming_perfect_second(I, answ, x, y, y_err, n_becoming_perfect, name, ftype, Parametr, n, answ_2, ratio, max_width, N_cutting, N_fragmentation, dpi_picture, dots_size, I_star = 0):  
    
    p0 = np.zeros(2*n+2)                                        #take ideal start conditions either from becoming_perfect or becoming_perfect_second
    if(I == 0):
        if(n > n_becoming_perfect):
            for i in range(2*n_becoming_perfect + 2):  
                p0[i] = answ[i]
            for i in range(2*n_becoming_perfect + 2, 2*n+2):  
                p0[i] = 1
        else:
             for i in range (2*n+2):  
                p0[i] = answ[i]
    else:
        for i in range(2*n+2):  
             p0[i] = answ_2[i]
        
    fun = lambda pp: (y - sin(x, pp, n))/y_err       #core of least squares
    ans = spo.leastsq(fun, p0, full_output=1)
    sigma = np.sqrt(np.sum((y - sin(x, ans[0], n))**2)/len(x))
    error = np.sqrt(np.diag(ans[1]*sigma))
    
    order_Error = -int(np.log10(error[1]))+1    
    T_ideal = ans[0][1]  
    Number_of_elements = len(x)
    ans_id = ans[0]
 
    Number_periods = (x - x[0])/T_ideal                         #To cuild phase curve (need to rearrange x array)
    X_E = np.zeros(Number_of_elements)
    y_max = y[0]
    for i in range(Number_of_elements): 	                    #find y max
        if (y[i] > y_max):
            y_max = y[i]
            I_max = i
    for i in range(Number_of_elements):
        X_E[i] = (x[i] - x[0]) - int(Number_periods[i])*T_ideal
    delta = X_E[I_max]
    for i in range(Number_of_elements):                         #shift plot so it starts from the minimum
        X_E[i] -= delta
        if (X_E[i] < 0):
            X_E[i] += T_ideal   
       
    save_path = path_file + '\\Results\\' + name + '\\'
    
    A = max(x) - min(x)
    B = max(y) - min(y)
    N_periods = np.round(((max(x) - min(x))/T_ideal), 1)
            
    hfont = {'fontname':'Helvetica'}
    fig = plt.figure(3 + I + I_star*(N_cutting + 2))  
    plt.gca().invert_yaxis()            
    fig.set_size_inches(ratio*7, 7)                 
    plt.xlim(-0.02, 1.02)
    plt.ylim((max(y) + 0.05*B), (min(y) - 0.1*B))
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    strin = 'Phase (T = ' + str(np.round(2457000 + x[I_max], 5)) + ' +' + str(np.round(T_ideal, 5)) + '*E)'
    plt.xlabel(strin, fontsize = 20, **hfont)
    plt.ylabel('V, mmag', fontsize = 20, **hfont)
    plt.plot(X_E/T_ideal, y,  color = 'green', linestyle = '', marker = '.', markersize = dots_size)
    plt.text(0, (np.min(y) - 1/30*B), name, fontsize = 20, **hfont)
    plt.savefig(save_path + name + "phase curve " + str(I) + ".png", dpi = dpi_picture)
    plt.show()
    
    NName = name + " phase curve" + str(I) + ftype                      #save data in the same file type
    completeName = os.path.join(save_path, NName) 
    with open(completeName, 'w+') as f:
        for i in range(Number_of_elements):
            f.write(str(X_E[i]/T_ideal) + ' ' + str(y[i]) + ' ' + str(y_err[i]) + '\n')

    X_0 = list(x)                   #copy arrays
    y_0 = list(y)
    y_0_error = list(y_err)
 
    k=0
    #d = max_width  / (1 + (max_width - Parametr)* I/((N_cutting-1)*Parametr))    #width of cutting     
    if (I < (int(N_cutting/2)-1)):
        d = max_width - (max_width - Parametr)*I/(int(N_cutting/2))
    else:
        N_times = N_cutting - int(N_cutting/2) + 1
        d_0 = Parametr + (max_width - Parametr)*2/(int(N_cutting/2))
        d = d_0 - (d_0 - Parametr) * (I + 2 - int(N_cutting/2))/N_times
    
    if not (I == N_cutting):
        for i in range(Number_of_elements):                             #for each dot  
            left = -d                                                   #find in the range from -d to d
            right = d
            for m in range(3):                                          #3 iterations
                x_delta = np.linspace((left), (right), N_fragmentation)
                x_delta_distance = (x_delta)**2 + ( ( y[i]  -  sin((x[i] + x_delta*A/(ratio*N_periods) ), ans_id, n)  ) / B  )**2           #distance between dot and approximation that correspond to the size of phase diagramm. It is calculating for N_fragmentation dots m time
                Minimum = x_delta_distance[0]
                J = 0
                for j in range(1, N_fragmentation):             #find a minimum through distances
                    if(x_delta_distance[j] < Minimum):
                        Minimum = x_delta_distance[j]
                        J = j
                if (J == 0):
                     left = x_delta[0]
                     right = x_delta[1]
                elif (J == (N_fragmentation-1)):
                     left = x_delta[N_fragmentation-2]
                     right = x_delta[N_fragmentation-1]
                else:
                    left = x_delta[J-1]
                    right = x_delta[J+1]
                    
            if (np.sqrt(Minimum) > d):                          #if dot's deviation is bigger than d => cut off
                del X_0[i-k]
                del y_0[i-k]
                del y_0_error[i-k] 
                k+=1       
    
    x = np.array(X_0)
    y = np.array(y_0)
    y_err = np.array(y_0_error)
    return np.round(T_ideal, order_Error), np.round(error[1], order_Error), x, y, y_err, ans_id

"""==========================================="""
"""COMPUTING APPROXIMATE VALUE OF PERIOD"""
"""==========================================="""

def Approximation_T(XxX, YyY, YyY_err, A, n_app_T, edge_appr_T, TTT_max, TTT_min, Presize_appr_T, name, dpi_picture, I = 0, N_cutting = 0):
    
    T_minnn = TTT_min
    T_maxxx = TTT_max                       #range of periods that can be (not take T_min=0)
    N_N = int(T_maxxx/Presize_appr_T)       #number of dots in this area
    X_minn = 0                              #just for fun(do not change)

    def sin2(t, T_Tt, pp, nn):              #approximation of function that take x data, period and parametrs and give the approximation function
        x = np.zeros(len(t))                #make array x lenth of x-data and full zero  
        x += pp[0] 
        for i in range(nn):                 #additions in Fourie series
            x += pp[2*i + 1]*np.sin(2*np.pi*t/T_Tt*(i+1)+pp[2*i+2])
        return x                            #return tha value of approximation function

    def sigma(XxXx, YyYy, YyYy_err, T_Tt, p00, nn):                              #function to find the sum of squares for each T
        fun = lambda pp: (YyYy - sin2(XxXx, T_Tt, pp, nn))/YyYy_err              #core of least squares
        ans = spo.leastsq(fun, p00, full_output=1)
        Sigma = np.sum((YyYy-sin2(XxXx, T_Tt, ans[0], nn))**2)/len(XxX)          #ans[0] - parametrs: amplitudes and phases
        return Sigma, ans[0]
    
    p0 = []                     
    p0.append(0)
    p0.append(A)
    for i in range(2*n_app_T-1):
        p0.append(1)
    YyY_sigma = []                #arrays for finding minimum
    XxX_sigma = []
    
    fig = plt.figure(1 + I * (N_cutting + 2))
    fig.set_size_inches(20, 6)
    T_T = np.linspace(T_minnn, T_maxxx, N_N)   
    save_path = path_file + '\\Results\\' + name + '\\'
    
    for i in T_T:                                                           #for each dot
        XxX_sigma.append(i)                                                 #fill x-array
        if(i == T_minnn):
            VALUE, PP0 = sigma(XxX, YyY, YyY_err, i, p0, n_app_T)           #find y and ideal parametrs
        else:
            VALUE, PP0 = sigma(XxX, YyY, YyY_err, i, PP0, n_app_T)          #start condition = ideal for previous
        YyY_sigma.append(VALUE)
    plt.plot(XxX_sigma, YyY_sigma, color = '#FF0000', ls = '-', lw = 2) 
    plt.savefig(save_path + name + "periodogram.png", dpi = dpi_picture)        
    
    
    value_error = False
    if ((np.min(YyY_sigma)/np.max(YyY_sigma)) < 0.3):
        value_error = True
    if value_error:
        YYY_min = min(YyY_sigma)                                                #find absolut minimum and normal periodogram
        for i in range(N_N):
            YyY_sigma[i] = (YyY_sigma[i] - YYY_min)/np.max(YyY_sigma)
        
        value = False                                                           #find first minimum (because others one can be 2T , 3T  etc.)
        for i in range(N_N):
            if (YyY_sigma[i] < edge_appr_T) and (not value):
                Index_1 = i
                value = True
            if value:
                if (YyY_sigma[i] > edge_appr_T):
                    Index_2 = i
                    break
            if (i == (N_N - 1)):            #for manual
                Index_2 = i       
    
        Y_minn = YyY_sigma[Index_1]                                             #find exact minimum
        for i in range(Index_1, Index_2):
            if (YyY_sigma[i] < Y_minn):
                Y_minn = YyY_sigma[i]
                X_minn = XxX_sigma[i]
                VALUE, PP0 = sigma(XxX, YyY, YyY_err, XxX_sigma[i], p0, n_app_T)
    
        local_delta = (T_maxxx-T_minnn)/N_N
        order_ld = -int(np.round(np.log10(local_delta)))+1
        Error_program = 0
        
    else:
        X_minn = 0
        order_ld = 0
        local_delta = 0
        PP0 = 0
        Error_program = 1
    #print(np.round(X_minn, order_ld), '+-', np. round(local_delta, order_ld))  
    return np.round(X_minn, order_ld), np. round(local_delta, order_ld), PP0, Error_program

"""==========================================="""
"""CREATING WINDOW AND GENERAL WIDJETS"""               #graph interphase
"""==========================================="""
def Manual_work():
    
    """==========================================="""
    """MAIN FUNCTION FOR MANUAL MODE"""
    """==========================================="""
    def do_for_single_star():
        
        start_time = time.time()                        #time of begining
        enttime[1].delete(0, len(enttime[1].get()))     #clear windows     
        enttime[0].delete(0, len(enttime[0].get()))
        entT.delete(0, len(entT.get()))
        entdT.delete(0, len(entdT.get()))
    
        name = ent_StarName.get()                       #get values that input user
        ftype = ent_TypeFile.get()
        T = float(init_param[0].get())                  #they strings
        dT = float(init_param[1].get())
        Parametrs_file = init_param[2].get()
                                                        #read parametrs
        n_app_T, n_becoming_perfect, edge_appr_T, Parametr, TT_min_par, Presize_appr_T, ratio, N_cutting, n_bec_per_sec, max_width, N_fragmentation, dpi_picture, dots_size = read_parametrs(Parametrs_file) 
        res = 'Iteration  Period\n'             #for writing results
        cut_res = 'Iteration   % of cutted dots \n'
        
        if (not os.path.exists('Results')):      # Create target Directory
            os.mkdir('Results')       
        sub_name = path_file + '\\Results\\' + name
        if (not os.path.exists(sub_name)):
            os.mkdir(sub_name)
        
        x, y, y_err, Number_of_elements0 = read_data(name, ftype)       #getting data from file
        A0 = (max(y)-min(y)) / 2                                                            #approximate amplitude
        Period, Period_error, ans_start, Error_program = Approximation_T(x, y, y_err, A0, n_app_T, edge_appr_T, (T + dT), (T - dT), Presize_appr_T, name, dpi_picture) #approximate period
        if not (Error_program):
            ans_ideal, T, ΔT = becoming_perfect(Period, A0, x, y, y_err, n_becoming_perfect, name, n_app_T, ans_start, dpi_picture, dots_size)                         #more presize period
            res += '    ' + str(0) + '      ' + str(T) + '   ' + str(ΔT) + '\n'                 #write basic iteration
            cut_res += '   ' + str(0) + '     ' + str(0) + '\n'
            ans_ideal_2 = 1
            T_true = 0
            T_array = np.zeros(N_cutting + 1)
            
            for indicator in range(N_cutting + 1):          #N_cutting times cut phase diagram
                T, ΔT, x, y, y_err, ans_ideal_2 = becoming_perfect_second(indicator, ans_ideal, x, y, y_err, n_becoming_perfect, name, ftype, Parametr, n_bec_per_sec, ans_ideal_2, ratio, max_width, N_cutting, N_fragmentation, dpi_picture, dots_size)
                if not (indicator == N_cutting):
                    res += '    ' + str(indicator + 1) + '      ' + str(T) + '     ' + str(ΔT) + '\n'
                    Number_of_elements = np.round((1 - len(x)/Number_of_elements0)*100, 1)
                    cut_res += '   ' + str(indicator + 1) + '            ' + str(Number_of_elements) + '\n'
                    T_true += T
                    T_array[indicator] = T
                
            T_true = T_true/(N_cutting)
            Ssigma = 0
            for indicator in range(N_cutting):
                Ssigma += (T_array[indicator] - T_true)**2
            Ssigma = 3*np.sqrt(Ssigma/(N_cutting*(N_cutting-1)))
            order_Error = -int(np.log10(Ssigma))+1    
            res += '    Period: ' +  str(np.round(T_true, order_Error)) + ' +- ' + str(np.round(Ssigma, order_Error)) + '\n'
            res += '\n'
            
            results_path =  path_file + '\\Results\\' + 'results_' + name + '.dat'
            cutted_dots_path = path_file + '\\Results\\' + 'Number_of_cutted_dots_' + name + '.dat'
            with open(results_path, 'w') as f:
                f.writelines(res)
            with open(cutted_dots_path, 'w') as f:
                f.writelines(cut_res)
            entT.insert(0, str(np.round(T_true, order_Error)))
            entdT.insert(0, str(np.round(Ssigma, order_Error)))
            t_0 = time.time() - start_time
            enttime[1].insert(0, str(round(t_0)-60*int(t_0/60)))    #fill windows with the values
            enttime[0].insert(0, str(int(t_0/60)))
        else:
            Error_1()        
        
    """==========================================="""
    """CLEARING WINDOW"""
    """==========================================="""
    
    def clear_win():                                            #realising the button "Clear"
        for i in range(3):
            init_param[i].delete(0, len(init_param[i].get()))
        ent_TypeFile.delete(0, len(ent_TypeFile.get()))
        ent_StarName.delete(0, len(ent_StarName.get()))
        entT.delete(0, len(entT.get()))
        entdT.delete(0, len(entdT.get()))
        for i in range(2):
            enttime[i].delete(0, len(enttime[i].get()))    
        
    window = tnk.Tk()                                           #start window
    bcg_cl = '#9999FF'                                          #background color
    window.title("Period D&P V5.2")                             #title of the window
    w = 900                                                     #width and height
    h = 350
    window.geometry(str(w) + 'x' + str(h))                      #set size
    window.config(bg=bcg_cl)    
    window.resizable(width=False, height=False)                 #not to change size
    
    lb_head = tnk.Label(window, font = ('Algerian', 19), text = 'Name of star:', bg=bcg_cl)                 #words and labels
    
    lb_TypeFile = tnk.Label(window, font = ('Bookman Old Style', 14), text = 'Type of file:', bg=bcg_cl)
    lb_par = tnk.Label(window, font = ('Algerian', 19), text = 'Parameters:', bg=bcg_cl)
    ent_StarName = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)                          #Entry_window
    ent_TypeFile = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    
    lb_head.place(x = 20, y = 30)               #their place on the window
    lb_TypeFile.place(x = 40, y = 70)
    ent_StarName.place(x = 250, y = 35)
    ent_TypeFile.place(x = 250, y = 72)
    lb_par.place(x = 20, y = 135)
    
    text_init = ['T approximately', 'Error of Τ', 'File with parameters']
    init_param = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12) for i in range(3)]
    init_param_lables = [tnk.Label(window, font = ('Century', 14), text = text_init[i], bg=bcg_cl) for i in range(3)]
    for i in range(3):
        init_param_lables[i].place(x = 40, y = 170 + i * 35)
        init_param[i].place(x = 250, y = 175 + i * 35)
    
    init_param[2].insert(0, 'Parametrs.txt')
        
    lb_Results = tnk.Label(window, font = ('Algerian', 19), text = 'Results:', bg=bcg_cl)
    lbT = tnk.Label(window, font = ('Century', 14), text = 'Period', bg=bcg_cl)
    lbpm = tnk.Label(window, font = ('Century', 14), text = '+-', bg=bcg_cl)
    entT = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 8)
    entdT = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 8)
    time_text = ['Time of calculations', 'min', 's']
    lbtime = [tnk.Label(window, font = ('Century', 14), text = time_text[i], bg=bcg_cl) for i in range(3)]
    enttime = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 7) for i in range(2)]
    
    lb_Results.place(x= 500, y = 70)
    lbT.place(x = 520, y = 110)
    lbpm.place(x = 700, y = 140)
    entT.place(x = 600, y = 140)
    entdT.place(x = 730, y = 140)
    lbtime[0].place(x = 520, y = 200)
    lbtime[1].place(x = 690, y = 230)
    lbtime[2].place(x = 825, y = 230)
    enttime[0].place(x = 600, y = 230)
    enttime[1].place(x = 735, y = 230)
    
    btn1 = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Calculate', bg = "blue", fg = "white", height = 1, width = 13, command = do_for_single_star)  #realising buttons
    btn2 = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Clear', bg = "red", fg = "white", height = 1, width = 13, command = clear_win)                
    btn1.place(x = 200, y = 300)
    btn2.place(x = 20, y = 300)
    window.mainloop()
    
def Automatic_work():
        
    window = tnk.Tk()
    bcg_cl = '#9999FF'
    window.title("Period D&P V5.2")
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
    ent_TaskFile.insert(0, 'Task1.txt')
    ent_Tmax.insert(0, '10')
    
    """==========================================="""
    """MAIN FUNCTION FOR AUTOMATIC MODE"""
    """==========================================="""
    
    def automatic_regime():
        enttime[0].delete(0, len(enttime[0].get()))
        enttime[1].delete(0, len(enttime[1].get()))
        ent_progress_1.delete(0, len(ent_progress_1.get()))
        ent_progress_2.delete(0, len(ent_progress_2.get()))

        task_file = ent_TaskFile.get()
        TT_max = float(ent_Tmax.get())
        Parametrs_file = ent_Par_file.get()
        n_app_T, n_becoming_perfect, edge_appr_T, Parametr, TT_min_par, Presize_appr_T, ratio, N_cutting, n_bec_per_sec, max_width, N_fragmentation, dpi_picture, dots_size = read_parametrs(Parametrs_file)
        
        with open(task_file, 'r') as f:
            s = f.readlines()
        N_stars = len(s)       
        i = 0
        while (i < len(s)):
            if (s[i][0] == "#"):
                del s[i]
            i += 1       
        res = '  Name     T        time'             
        start_time_0 = time.time()
        ent_progress_1.insert(0, '0')
        ent_progress_2.insert(0, str(len(s)))
        if (not os.path.exists('Results')):      # Create target Directory
                os.mkdir('Results')
        
        for i in range(N_stars):
            #try:
            start_time = time.time()
            a = s[i].split()
            name = a[0]
            ftype = a[1]
            res += name + '\n'
            
            sub_name = path_file + '\\Results\\' + name
            if (not os.path.exists(sub_name)):
                os.mkdir(sub_name)

            x, y, y_err, Number_of_elements0 = read_data(name, ftype)
            A0 = (max(y)-min(y)) / 2
            Tappr, Terr, ans_start, Error_program = Approximation_T(x, y, y_err, A0, n_app_T, edge_appr_T, TT_max, TT_min_par, Presize_appr_T, name, dpi_picture, i, N_cutting)
            if not Error_program:
                ans_ideal, T, ΔT = becoming_perfect(Tappr, A0, x, y, y_err, n_becoming_perfect, name, n_app_T, ans_start, dpi_picture, dots_size, i, N_cutting)
                arrT = []
                arrT.append(T)
                ans_ideal_2=[]
                
                T_true = 0
                K_index = 0
                for indicator in range(N_cutting + 1):
                    T, ΔT, x, y, y_err, ans_ideal_2 = becoming_perfect_second(indicator, ans_ideal, x, y, y_err, n_becoming_perfect, name, ftype, Parametr, n_bec_per_sec, ans_ideal_2, ratio, max_width, N_cutting, N_fragmentation, dpi_picture, dots_size, i)               
                    if not (T == arrT[K_index]):
                        arrT.append(T)
                        K_index += 1
                        T_true += T
                    
                T_true = T_true/(K_index+1)
                Ssigma = 0
                for indicator in range(K_index+1):
                    Ssigma += (arrT[indicator] - T_true)**2
                Ssigma = 3*np.sqrt(Ssigma/(K_index*(K_index+1)))
                  
                t0 = time.time() - start_time   
                
                order_Error = -int(np.log10(Ssigma))+1    
                res += str(np.round(T_true, order_Error)) + ' +- ' + str(np.round(Ssigma, order_Error)) + '    '+ str(int(t0/60)) + 'min ' + str(round(t0)-60*int(t0/60)) + 's' + '\n'                
                start_time = t0
                k = int(ent_progress_1.get())
                ent_progress_1.delete(0, len(ent_progress_1.get()))
                ent_progress_1.insert(0, str(k+1))
                
            else:
                res += name + '   Error #1\n'
                t0 = time.time() - start_time   
                start_time = t0
            #except: 
             #   print("Problem with " + str(i+1) + " star. Please check in manual mode")
              #  res += 'Problem. Check ' + name + 'manually'
               # res += '\n'
                
        task_file = str(task_file)
        for j in range(len(task_file)):
            if task_file[len(task_file) - j -1] == '.':
                task_file = task_file[:(len(task_file) - j - 1)]
                break
        
        results_path =  path_file + '\\Results\\' + 'results_' + task_file + '.dat'
        with open(results_path, 'w') as f:
            f.writelines(res)
        t_0 = time.time() - start_time_0  
        enttime[1].insert(0, str(round(t_0)-60*int(t_0/60)))
        enttime[0].insert(0, str(int(t_0/60)))
        
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
    
    btn = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Calculate', bg = "blue", fg = "white", height = 1, width = 13, command = automatic_regime)
    btn.place(x = 50, y = 255) 
    window.mainloop()

window_0 = tnk.Tk()
bcg_cl = '#9999FF'
window_0.title("Period D&P V5.2")
w = 390
h = 100
window_0.geometry(str(w) + 'x' + str(h))
window_0.config(bg=bcg_cl)
window_0.resizable(width=False, height=False)

btn1 = tnk.Button(window_0, font = ('Bookman Old Style', 14), text = 'Manual Work', bg = "blue", fg = "white", height = 1, width = 13, command = Manual_work)
btn2 = tnk.Button(window_0, font = ('Bookman Old Style', 14), text = 'Automatic Work', bg = "green", fg = "white", height = 1, width = 13, command = Automatic_work)
btn1.place(x = 200, y = 30)
btn2.place(x = 20, y = 30)
window_0.mainloop()