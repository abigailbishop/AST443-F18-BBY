# Uses error to determine significant digits to display data in
# Usage:  r'$(%.0f \pm %.0f)x10^{%i}$' % sigfigs(value, error)
def sigfigs(value, error):
    sigfigs = 0
    if error >= 1.:
        e = str(int(error))
        for i in range(len(e),0,-1):
            if int(e[i-1]) != 0:
                sigfigs = len(e) - i
                break
        value = value / 10.**sigfigs
        error = error / 10.**sigfigs
        return(value, error, sigfigs)
    elif error == 0.:
        return(value, 0, 0)
    else:
        e = str(float(error))
        for i in range(2,len(e),1):
            if int(e[i]) != 0:
                sigfigs = 1 - i
                break
        value = value / 10.**sigfigs
        error = error / 10.**sigfigs
        return(value, error, sigfigs)
