# Returns index of calculated baseline in the np.loadtxt file based on 
#    ladder baseline
def get_baseline(baseline_ladder):
    baseline_ladder = float(baseline_ladder)
    if baseline_ladder == 30. :
        return 0
    elif baseline_ladder == 35.:
        return 1
    elif baseline_ladder == 40.:
        return 2
    elif baseline_ladder == 45.:
        return 3
    elif baseline_ladder == 50.:
        return 4
    elif baseline_ladder == 55.:
        return 5
    elif baseline_ladder == 60.:
        return 6
    elif baseline_ladder == 65.:
        return 7
    elif baseline_ladder == 70.:
        return 8
    elif baseline_ladder == 77.:
        return 9
    else:
        print("Baseline value doesn't exist. Please check your file")
        print("Attempted baseline: {}".format(baseline_ladder))
        return None
