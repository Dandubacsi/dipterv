import numpy as np

def TableLookup1D(xdata,ydata,xvalues):
    out = np.zeros(len(xvalues))
    for i, x in enumerate(xvalues):
        if x>max(xdata) or x<min(xdata):
            print('Interpolation point (index = '+str(i)+') is not inside the given datapoints')
        for j in range(len(xdata)):
            if xdata[j]>x:
                index = j
                break
        out[i] = (ydata[index]-ydata[index-1])/(xdata[index]-xdata[index-1])*(x-xdata[index-1])+ydata[index-1]
    return out

def PrefixCorrection(value):
    if value[-1] == 'f':
        return float(value[0:-1])*1e-15
    if value[-1] == 'p':
        return float(value[0:-1])*1e-12
    if value[-1] == 'n':
        return float(value[0:-1])*1e-9
    if value[-1] == 'u':
        return float(value[0:-1])*1e-6
    if value[-1] == 'm':
        return float(value[0:-1])*1e-3
    if value[-1] == 'k':
        return float(value[0:-1])*1e+3
    if value[-1] == 'M':
        return float(value[0:-1])*1e+6
    if value[-1] == 'G':
        return float(value[0:-1])*1e+9
    if value[-1] == 'G':
        return float(value[0:-1])*1e+12
    if value[-1] == 'T':
        return float(value[0:-1])*1e+15
    return float(value)
